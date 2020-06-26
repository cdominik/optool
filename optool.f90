subroutine usage()
  write(*,'("")')
  write(*,'("===============================================================================")')
  write(*,'("optool uses DHS (Min+2005) and a new effective medium solver (Lefèvre,Min+2020)")')
  write(*,'("to compute dust opacities and scattering properties for a size distribution of")')
  write(*,'("grains in a wavelength interval.  Without arguments, it produces the DIANA")')
  write(*,'("standard opacities (Woitke,Min+ 2016). With arguments, set your composition and")')
  write(*,'("other parameters.")')
  write(*,'("")')
  write(*,'("-c ?                      List available materials")')
  write(*,'("-c KEY-or-FILE [Mfrac]    Add material with mass fraction to core")')
  write(*,'("-m KEY-or-FILE [Mfrac]    Set material and mass fraction of mantle")')
  write(*,'("-p POROSITY [P_MANTLE]    Set porosity, possibly different for core and mantle")')
  write(*,'("-fmax VHMAX               Maximum volume fraction of vacuum in DHS computation")')
  write(*,'("-a AMIN AMAX [APOW [NA]]  Set up grain radius radius grid (unit: micron)")')
  write(*,'("-l LMIN LMAX [NLAM]       Set up wavelength grid          (unit: micron)")')
  write(*,'("-l FILE                   Read wavelength grid from file, e.g. some/file.lnk")')
  write(*,'("-t [TMIN [TMAX [NT]]]     Compute mean opacities on this temperature grid")')
  write(*,'("-d [NSUB]                 Write NA files for specific grain sizes")')
  write(*,'("-s                        Add the scattering matrix to the produced output")')
  write(*,'("-radmc [LABEL]; -fits     Special output options")')
  write(*,'("-o [DIRECTORY]            Output to DIRECTORY instead of current working dir")')
  write(*,'("-h                        Show this message")')
  write(*,'("===============================================================================")')
end subroutine usage

module Defs
  ! Module with constants and the basic data structures used in the program
  implicit none
  integer, public, parameter     :: dp = selected_real_kind(P=15)
  ! ----------------------------------------------------------------------
  ! Physics and math constants
  ! ----------------------------------------------------------------------
  real (kind=dp), public, parameter :: pi = 3.1415926535897932384_dp
  ! ----------------------------------------------------------------------
  ! Some global switches
  ! ----------------------------------------------------------------------
  logical, public                :: verbose   = .false. ! additional output
  logical, public                :: blendonly = .false. ! only blend materials and write result out
  logical, public                :: CLBlend   = .true.  ! use Charléne Lefévre's new blender
  logical, public                :: split     = .false. ! split to many files
  ! ----------------------------------------------------------------------
  ! Lambda is shared, because multiple routines need it
  ! ----------------------------------------------------------------------
  real (kind=dp), allocatable    :: lam(:)     ! wavelength array
  integer                        :: nlam       ! nr of wavelength points
  ! ----------------------------------------------------------------------
  ! Material properties
  ! ----------------------------------------------------------------------
  integer                      :: mat_nm       ! number of materials specified
  integer,        allocatable  :: mat_num(:)   ! number associated to a component
  character*500,  allocatable  :: mat_loc(:)   ! either 'core' or 'mantle'
  character*500,  allocatable  :: mat_lnk(:)   ! the lnk key of file path
  real (kind=dp), allocatable  :: mat_rho(:)   ! specific mass density of material
  real (kind=dp), allocatable  :: mat_mfr(:)   ! mass fraction of each component
  real (kind=dp), allocatable  :: mat_e1(:,:)  ! Real      part of refractive index
  real (kind=dp), allocatable  :: mat_e2(:,:)  ! Imaginary part of refractive index
  ! ----------------------------------------------------------------------
  ! Mueller matrix structure, records only non-zero elements of the matrix
  ! ----------------------------------------------------------------------
  type mueller
     real (kind=dp) :: F11(180),F12(180),F22(180),F33(180),F44(180),F34(180)
  end type mueller
  ! ----------------------------------------------------------------------
  ! The particle structure, contains particle and scatteting properties
  ! ----------------------------------------------------------------------
  type particle
     real (kind=dp)              :: rv,rvmin,rvmax ! grain radius with min and max
     real (kind=dp)              :: rho            ! mass density in g.cm^-3
     real (kind=dp), allocatable :: Kabs(:),Ksca(:),Kext(:) ! Opacities
     real (kind=dp), allocatable :: g(:)           ! asymmetry, Henyey Greenstein
     TYPE(MUELLER),  allocatable :: F(:)           ! Mueller matrix elements
  end type particle
  ! ----------------------------------------------------------------------
  ! The output directory
  ! ----------------------------------------------------------------------
  character*500                  :: outdir = ''    ! Output directory 
  ! ----------------------------------------------------------------------
  ! Interfaces
  ! ----------------------------------------------------------------------
end module Defs

! ----------------------------------------------------------------------
! ----------------------------------------------------------------------
! ----------------------------------------------------------------------
!
!                         optool, the main program
!
! ----------------------------------------------------------------------
! ----------------------------------------------------------------------
! ----------------------------------------------------------------------

program optool
  use Defs
  use omp_lib
  implicit none
  integer         :: na              ! nr of sizes for size distribution
  real (kind=dp)  :: amin,amax       ! min and max size of grains
  real (kind=dp)  :: apow            ! power law index f(a) ~ a^(-apow)
  real (kind=dp)  :: fmax            ! maximum fraction of vaccum for DHS
  real (kind=dp)  :: p_core, p_mantle! porosity for core and mantle

  real (kind=dp)  :: lmin,lmax       ! min and max wavelength

  logical         :: write_mean_kap  ! Should mean kappas be computed?
  logical         :: write_scatter   ! Should a fits file be written?
  logical         :: write_fits      ! Should a fits file be written?
  logical         :: for_radmc       ! Should the scattering matrix use RADME-3D convention?
  real (kind=dp)  :: tmin,tmax       ! min and max temperature for mean opacities
  integer         :: nt              ! number of temperature steps

  integer         :: nm              ! nr of grain materials

  type(particle)  :: p
  integer         :: i,ndone         ! counter
  integer         :: im,ia           ! for material, radius
  character*100   :: tmp
  character*100   :: value

  logical         :: have_mantle
  logical         :: arg_is_value, arg_is_number  ! functions to test arguments
  character*500   :: fitsfile,meanfile            ! file names for output
  character*500   :: make_file_path               ! Function
  character*50    :: radmclbl   = ""              ! file label for RADMC-3D compatible
  character*50    :: label                        ! for use in file names

  real (kind=dp)  :: asplit,afact,afsub,amaxsplit,aminsplit
  integer         :: nsubgrains = 5,nsub

  real (kind=dp), allocatable :: e1d(:),e2d(:)

  ! ----------------------------------------------------------------------
  ! Defaults values for parameters and switches
  ! ----------------------------------------------------------------------
  amin           = 0.05       ! micrometer
  amax           = 3000.      ! micrometer
  apow           = 3.50_dp    ! a minus sign will be added internally
  na             = 50

  lmin           = 0.05_dp    ! micrometer
  lmax           = 10000.0_dp ! micrometer
  nlam           = 300

  tmin           = 10d0       ! K
  tmax           = 1d4        ! K
  nt             = 200

  fmax           = 0.8_dp     ! volume fraction DHS
  p_core         = 0.0_dp     ! porosity core
  p_mantle       = 0.0_dp     ! porosity mantle

  nm             = 0          ! number of materials - zero to start with

  write_fits     = .false.    ! Default is to write ASCII output
  write_scatter  = .false.    ! Default is to not write scattering matrix
  write_mean_kap = .false.    ! Default is to not compute mean kappas
  for_radmc      = .false.    ! Default is to use optool conventions.

  ! ----------------------------------------------------------------------
  ! Allocate space for up to 12 different materials
  ! ----------------------------------------------------------------------
  allocate(mat_num(12))
  allocate(mat_loc(12))
  allocate(mat_lnk(12))
  allocate(mat_mfr(12))
  allocate(mat_rho(12))
  have_mantle = .false.

  ! ----------------------------------------------------------------------
  ! Initialize rho, because we need the fact that it has not been set
  ! to decide what to do with lnk files where it is missing
  ! ----------------------------------------------------------------------
  do im=1,12
     mat_rho(im) = 0.d0
  enddo

  ! ----------------------------------------------------------------------
  ! Process the command line arguments
  ! ----------------------------------------------------------------------

  ! Loop over all command line arguments
  call getarg(1,tmp)
  i = 1
  do while(tmp.ne.' ')
     select case(tmp)

        ! ----------------------------------------------------------------------
        ! Definition of a material for the mix
        ! ----------------------------------------------------------------------
     case('-c','-m')
        i  = i+1
        nm = nm+1; mat_num(nm) = nm

        ! First value is the material key or refindex file path
        call getarg(i,value); read(value,'(A)') mat_lnk(nm)

        if (value .eq. '?') then
           call ListBuiltinMaterials()
           stop
        endif
        ! Second value is the volume fraction
        if (.not. arg_is_value(i+1)) then
           print *, "WARNING: 1.0 used for missing mass fraction of material: ",trim(mat_lnk(nm))
           mat_mfr(nm) = 1.0
        else
           i = i+1; call getarg(i,value); read(value,*) mat_mfr(nm)
        endif

        if (mat_mfr(nm).eq.0) then
           print *, "WARNING: Ignoring material with zero mass fraction: ",trim(mat_lnk(nm))
           nm = nm-1
        else
           ! Set the type, and make sure we have at most one mantle material
           if (tmp.eq.'-m') then
              ! This is the mantle material
              if (have_mantle) then
                 print *, "ERROR: Only one mantle material is allowed"; stop
              else
                 if (nm.eq.1) then
                    print *,"ERROR: at least one core material must be specified"; stop
                 endif
                 mat_loc(nm) = 'mantle'
                 have_mantle = .true.
              endif
           else
              ! This is a core material.
              if (have_mantle) then
                 print *,"ERROR: Mantle material must be specified last"; stop
              endif
              mat_loc(nm) = 'core'
           endif
        endif

        ! There might be a density, as a third argument
        if (arg_is_value(i+1)) then
           i = i+1; call getarg(i,value); read(value,*) mat_rho(nm)
        endif

        ! ----------------------------------------------------------------------
        ! Grain size setup
        ! ----------------------------------------------------------------------
     case('-a')
        ! ----------------------------------------------------------------------
        ! -a expects 2, 3, ro 4 values:  amin amax [na [apow]]
        ! ----------------------------------------------------------------------
        if (.not. arg_is_number(i+1) .or. .not. arg_is_number(i+2)) then
           print *,"ERROR: -a needs 2-4 values: amin amax [na [apow]]"
           stop
        endif
        i=i+1; call getarg(i,value); read(value,*) amin
        i=i+1; call getarg(i,value); read(value,*) amax
        ! Lets see if there is more, i.e. na
        if (arg_is_value(i+1)) then
           i=i+1; call getarg(i,value); read(value,*) na
           ! Lets see if there is more, i.e. apow
           if (arg_is_value(i+1)) then
              i=i+1; call getarg(i,value); read(value,*) apow
           endif
        endif
     case('-amin','--amin')
        i = i+1; call getarg(i,value); read(value,*) amin
     case('-amax','--amax')
        i = i+1; call getarg(i,value); read(value,*) amax
     case('-na')
        i = i+1; call getarg(i,value); read(value,*) na
     case('-apow','--apow')
        i = i+1; call getarg(i,value); read(value,*) apow

        ! ----------------------------------------------------------------------
        ! Wavelength setup
        ! ----------------------------------------------------------------------
     case('-l')
        ! ----------------------------------------------------------------------
        ! -l expects a file name, or 2-3 numbers: lmin lmax [nlam]
        ! ----------------------------------------------------------------------
        if (.not. arg_is_value(i+1)) then
           print *,"ERROR: -l needs a file or numbers as values"; stop
        else if (.not. arg_is_number(i+1)) then
           ! Could be a file name.  If yes, read the lambda grid from it
           call getarg(i+1,value)
           call check_for_file(trim(value))
           i=i+1
           call read_lambda_grid(trim(value))
        else if (.not. arg_is_value(i+2)) then
           ! First arg was a number, this one is not.  Bad.
           print *,"ERROR: -l needs 2-3 values: lmin lmax [nlam]"; stop
        else
           ! We have 2 numbers, these are lmin and lmax
           i = i+1; call getarg(i,value); read(value,*) lmin
           i = i+1; call getarg(i,value); read(value,*) lmax
           if (arg_is_value(i+1)) then
              i=i+1; call getarg(i,value); if (value.ne.'') read(value,*) nlam
           endif
        endif
     case('-lmin','--lmin')
        i = i+1; call getarg(i,value); read(value,*) lmin
     case('-lmax','--lmax')
        i = i+1; call getarg(i,value); read(value,*) lmax
     case('-nlam','--nlam','-nl','--nl')
        i = i+1; call getarg(i,value); read(value,*) nlam

        ! ----------------------------------------------------------------------
        ! Grain geometry
        ! ----------------------------------------------------------------------
     case('-p','-porosity','--porosity')
        i = i+1;  call getarg(i,value); read(value,*) p_core
        if (arg_is_value(i+1)) then
           i=i+1; call getarg(i,value); read(value,*) p_mantle
        else
           p_mantle = p_core
        endif
     case('-fmax')
        i = i+1; call getarg(i,value); read(value,*) fmax

        ! ----------------------------------------------------------------------
        ! Temperature setup for mean kappas
        ! ----------------------------------------------------------------------
     case('-t')
        write_mean_kap = .true.
        if (arg_is_number(i+1)) then
           i=i+1;        call getarg(i,value); read(value,*) tmin
           if (arg_is_number(i+1)) then
              i=i+1;     call getarg(i,value); read(value,*) tmax
              if (arg_is_number(i+1)) then
                 i=i+1;  call getarg(i,value); read(value,*) nt
              endif
           endif
        endif
     case('-fits')
        write_fits = .true.

        ! ----------------------------------------------------------------------
        ! Various other switches
        ! ----------------------------------------------------------------------
     case('-o')
        if (arg_is_value(i+1)) then
           i=i+1; call getarg(i,outdir)
        else
           outdir = 'output'
        endif
     case('-s','-scatter','-scat')
        write_scatter=.true.
     case('-radmc','-radmc3d')
        for_radmc = .true.
        if (arg_is_value(i+1)) then
           i=i+1
           call getarg(i,radmclbl)
        endif
     case('-d')
        split = .true.
        if (arg_is_value(i+1)) then
           i=i+1; call getarg(i,value); read(value,*) nsubgrains
        endif
     case('-b','-blendonly','--blendonly')
        ! Write blended refractive index to file and exit
        blendonly = .true.
     case('-B')
        ! Use the old blender
        CLBlend = .false.
     case('-v','-verbose','--verbose')
        verbose = .true.
     case('?','-h','--help','-help','help')
        call usage()
        stop
     case default
        write(*,*) "ERROR: Option or Arg: >",trim(tmp),'> not recognized'
        write(*,*) "For help, try: optool -h     ... or find the user guide OpTool.pdf"
        stop
     end select
     i = i+1
     call getarg(i,tmp)   ! Get the next argument for the loop.
  enddo

  ! ----------------------------------------------------------------------
  ! Sanity checks and preparations
  ! ----------------------------------------------------------------------
  if (split .and. write_mean_kap) then
     write(*,*) 'ERROR: With both -d and -t options, the dustkapmean.dat file'
     write(*,*) '       would only reflect the final size bin.'
     stop
  endif
  if (split .and. blendonly) then
     write(*,*) 'WARNING: Turning off -s for -blendonly'
     split = .false.
  endif
#ifndef USE_FITSIO
  if (write_fits) then
     write(*,*) 'ERROR: Support for writing FITS files needs to be compiled in'
     write(*,*) '       Try: "make clean", and then "make fitsio=true"'
     stop
  endif
#endif
  if (trim(outdir) .ne. '') then
     call make_directory(outdir)
  endif
  meanfile       = make_file_path(outdir,"dustkapmean.dat")
  fitsfile       = make_file_path(outdir,"dustkappa.fits")
  
  ! ----------------------------------------------------------------------
  ! Default grain composition if nothing is specified
  ! ----------------------------------------------------------------------
  if (nm.eq.0) then
     ! Set default composition will set a default composition here, the DIANA opacities
     write(*,'("No materials specified, using DIANA standard")')
     nm = 2
     mat_lnk(1) = 'pyr-mg70' ; mat_loc(1)  = 'core' ; mat_mfr(1)     = 0.87d0
     mat_lnk(2) = 'c-z'      ; mat_loc(2)  = 'core' ; mat_mfr(2)     = 0.1301d0
     p_core     = 0.25d0
  endif
  mat_nm   = nm

  if (nm .ge. 10) then
     write(*,*) 'ERROR: Too many materials'
     stop
  endif
  do im = 1, nm
     mat_loc(im) = trim(mat_loc(im))
     mat_lnk(im) = trim(mat_lnk(im))
  enddo

  ! ----------------------------------------------------------------------
  ! Write a setup summary to the screen
  ! ----------------------------------------------------------------------
  call write_header(6,'',amin,amax,apow,na,lmin,lmax, &
       p_core,p_mantle,fmax,mat_mfr,mat_nm)
  
  ! ----------------------------------------------------------------------
  ! Make a logarithmic lambda grid, unless read in from file
  ! ----------------------------------------------------------------------
  if (allocated(lam)) then
     ! Lam was allocated by reading from a file
     continue
  else
     allocate(lam(nlam))
     do i = 1,nlam
        lam(i)=10.0_dp**(log10(lmin)+log10(lmax/lmin)*(i-1)/(nlam-1))
     enddo
  endif
  ! Allocate all lambda-dependant arrays
  allocate(p%Kabs(nlam))
  allocate(p%Kext(nlam))
  allocate(p%Ksca(nlam))
  allocate(p%g(nlam))
  allocate(p%F(nlam))

  ! Allocate space for the refractive indices
  allocate(mat_e1(nm+1,nlam),mat_e2(nm+1,nlam))

  ! ----------------------------------------------------------------------
  ! Get the refractory index data for all materials
  ! ----------------------------------------------------------------------
  allocate(e1d(nlam),e2d(nlam))
  do im=1,nm
     call GetAndRegridLNK(mat_lnk(im),lam(1:nlam),e1d(1:nlam),e2d(1:nlam), &
          nlam,.true.,mat_rho(im))
     mat_e1(im,1:nlam)    = e1d(1:nlam)
     mat_e2(im,1:nlam)    = e2d(1:nlam)
  enddo
  deallocate(e1d,e2d)

  ! ----------------------------------------------------------------------
  ! Loop for splitting the output into files by grain size
  ! ----------------------------------------------------------------------
  if (split) then
     write(*,'("Computing opacities for ",I3," different grain size bins")') na
     nsub = nsubgrains
     if (mod(nsub,2).eq.0) nsub = nsub+1
     afact = (amax/amin)**(1.d0/real(na))   ! factor to next grain size
     afsub = afact**(1.d0/real(nsub-1))     ! Factor to next subgrain size
     ndone = 0; call tellertje(ndone,na)
     
     !$OMP parallel do if (split)                                                 &
     !$OMP default(none)                                                          &
     !$OMP shared(amin,afact,afsub,nsub,apow,fmax,p_core,p_mantle,mat_mfr,mat_nm) &
     !$OMP shared(lmin,lmax,write_scatter,for_radmc,write_fits,radmclbl,ndone,na) &
     !$OMP private(ia,asplit,aminsplit,amaxsplit,label,fitsfile,p)
     do ia=1,na
        asplit    = amin  *afact**real(ia-1+0.5)
        aminsplit = asplit*afsub**real(-nsub/2)
        amaxsplit = asplit*afsub**real(+nsub/2)
        call ComputePart(p,aminsplit,amaxsplit,apow,nsub,fmax,p_core,p_mantle, &
             mat_mfr,mat_nm,.false.)

        ! We do not allow the output to be parallel
        !$OMP critical
        ndone = ndone + 1
        call tellertje(ndone,na)
        write(label,'(I3.3)') ia
        if (write_fits) then
#ifdef USE_FITSIO
           write(fitsfile,'(A,"_",A,".fits")') "dustkappa",trim(label)
           fitsfile = make_file_path(outdir,fitsfile)
           call write_fits_file(p,aminsplit,amaxsplit,apow,nsub, &
                fmax,p_core,p_mantle,mat_nm,mat_mfr,fitsfile)
#endif
        else
           if (radmclbl .ne. ' ') label = trim(radmclbl) // "_" // label
           call write_ascii_file(p,aminsplit,amaxsplit,apow,nsub,lmin,lmax, &
                fmax,p_core,p_mantle,mat_mfr,mat_nm, &
                label,write_scatter,for_radmc,.false.)
        endif
        !$OMP end critical
     enddo
     !$OMP end parallel DO
     stop
  endif

  ! ----------------------------------------------------------------------
  ! Call the main routine to compute the opacities and scattering matrix
  ! ----------------------------------------------------------------------
  call ComputePart(p,amin,amax,apow,na,fmax,p_core,p_mantle,mat_mfr,nm,.true.)

  ! ----------------------------------------------------------------------
  ! Write the output
  ! ----------------------------------------------------------------------
  if (write_fits) then
#ifdef USE_FITSIO
     call write_fits_file(p,amin,amax,apow,nsub, &
          fmax,p_core,p_mantle,mat_nm,mat_mfr,fitsfile)
#endif
     continue
  else
     call write_ascii_file(p,amin,amax,apow,na,lmin,lmax, &
          fmax,p_core,p_mantle,mat_mfr,mat_nm, &
          radmclbl,write_scatter,for_radmc,.true.)
  endif

  ! ----------------------------------------------------------------------
  ! Produce and write mean opacities
  ! ----------------------------------------------------------------------
  if (write_mean_kap) then
     call mean_opacities(lam,nlam,p%kabs,p%ksca,p%g,tmin,tmax,nt,meanfile)
  endif

  deallocate(mat_num)
  deallocate(mat_loc)
  deallocate(mat_lnk)
  deallocate(mat_mfr)
  deallocate(mat_rho)
  deallocate(p%Kabs)
  deallocate(p%Kext)
  deallocate(p%Ksca)
  deallocate(p%g)
  deallocate(p%F)

end program optool

! ----------------------------------------------------------------------
! ----------------------------------------------------------------------
! ----------------------------------------------------------------------
!
!     ComputePart, the central routine avaraging properties over sizes
!
! ----------------------------------------------------------------------
! ----------------------------------------------------------------------
! ----------------------------------------------------------------------

subroutine ComputePart(p,amin,amax,apow,na,fmax,p_c,p_m,mfrac0,nm0,progress)
  ! ----------------------------------------------------------------------
  ! Main routine to compute absorption cross sections and the scattering matrix.
  !
  ! INPUT
  !   amin       Minimum grain size to consider, in microns
  !   amax       Maximum grain size to consider, in microns
  !   apow       The powerlaw exponent for the size distribution, f(a) ~ a^{-apow}
  !   na         The number of grains to consider between amin and amax
  !   fmax       The maximum volume fraction of vacuum for the DHS computations
  !   p_c        The porosity of the core, volume fraction of vacuum
  !   p_m        The porosity of the mantle, if there is one
  !   mfrac0     An array of nm mass fraction for the various materials.
  !   nm0        The number of materials, also the size of the array mfrac
  !   progress   When .true. give information about progress
  !
  ! OUTPUT
  !   p          A "particle" structure to return all the information in.
  !              See the module Defs for the definition.
  ! ----------------------------------------------------------------------
  use Defs
  use omp_lib
  implicit none

  real (kind=dp)                 :: amin,amax,apow   ! min and max grain size, and power law exp
  real (kind=dp)                 :: fmax             ! maximum fraction of vaccum for DHS
  real (kind=dp)                 :: p_c,p_m          ! porosity, core and mantle

  integer                        :: nm0,nm,im        ! nr of grain materials in a composite grain
  real (kind=dp),allocatable     :: rho(:)           ! specific material dnesities
  real (kind=dp)                 :: mfrac0(nm0)      ! mass fractions, input
  real (kind=dp),allocatable     :: mfrac(:)         ! mass fractions, modified
  real (kind=dp)                 :: vfrac_mantle     ! volume fraction of mantle
  real (kind=dp)                 :: mfrac_mantle     ! mass   fraction of mantle

  TYPE (PARTICLE)                :: p

  logical                        :: progress

  integer                        :: i,j              ! counters for various loops
  integer                        :: mn_max           ! maximum number of grain material
  integer, parameter             :: n_ang = 180      ! number of angles FIXME let the user set this?
  integer                        :: na               ! Number of grains sizes between amin and amax
  integer                        :: nf,if            ! Number of DHS volume fractions
  integer                        :: ns,is            ! Number of grains sizes
  integer                        :: ilam,il,ndone    ! Counter for wavelengths
  integer                        :: i_mantle         ! Index of mantle material, if any
  integer                        :: err,spheres,toolarge ! Error control for Mie routines

  real (kind=dp)                 :: cext, csca, cabs
  real (kind=dp)                 :: cext_ff, csca_ff, cabs_ff
  real (kind=dp)                 :: qext, qsca, qabs, gqsc

  real (kind=dp), allocatable    :: f11(:,:),  f12(:,:),  f22(:,:),  f33(:,:),  f34(:,:),  f44(:,:)
  real (kind=dp), allocatable    :: Mief11(:), Mief12(:), Mief22(:), Mief33(:), Mief34(:), Mief44(:)

  real (kind=dp), allocatable    :: mu(:)
  real (kind=dp), allocatable    :: M1(:,:), M2(:,:), S21(:,:), D21(:,:)
  real (kind=dp), allocatable    :: r(:),nr(:)       ! Grain radii and size distribution
  real (kind=dp), allocatable    :: f(:),wf(:)       ! Needed for Gauss-Legendre integration

  real (kind=dp)                 :: rmie, lmie
  real (kind=dp)                 :: e1mie, e2mie
  real (kind=dp)                 :: csmie, cemie
  real (kind=dp)                 :: theta

  real (kind=dp)                 :: aminlog,amaxlog,pow
  real (kind=dp)                 :: rad,r1,rcore
  real (kind=dp)                 :: tot,tot2,mtot,vtot
  real (kind=dp)                 :: mass,vol
  real (kind=dp)                 :: rho_av,rho_core,rho_mantle
  real (kind=dp)                 :: wvno

  ! Effective refractive index
  real (kind=dp)                 :: e1mg,       e2mg
  real (kind=dp), allocatable    :: e1(:,:),    e2(:,:)
  real (kind=dp), allocatable    :: e1blend(:), e2blend(:)
  real (kind=dp), allocatable    :: e1mantle(:),e2mantle(:)
  real (kind=dp), allocatable    :: vfrac(:)
  complex (kind=dp), allocatable :: epsj(:)
  complex (kind=dp)              :: m,min,eps_eff

  character (len=3)              :: meth   ! Method for computing opacities

  ns     = na    ! number of subgrains to compute
  meth   = 'DHS' !

  ! ----------------------------------------------------------------------
  ! Allocate the necessary arrays
  ! ----------------------------------------------------------------------
  mn_max = nm0+1 ! Allocate one more, because vacuum will also be a material
  allocate(Mief11(n_ang),Mief12(n_ang),Mief22(n_ang),Mief33(n_ang))
  allocate(Mief34(n_ang),Mief44(n_ang))
  allocate(mu(n_ang),M1(n_ang,2),M2(n_ang,2),S21(n_ang,2),D21(n_ang,2))

  allocate(vfrac(mn_max),mfrac(mn_max),rho(mn_max))
  allocate(epsj(mn_max))
  allocate(f11(nlam,n_ang),f12(nlam,n_ang),f22(nlam,n_ang))
  allocate(f33(nlam,n_ang),f34(nlam,n_ang),f44(nlam,n_ang))

  allocate(e1(mn_max,nlam),e2(mn_max,nlam))
  allocate(e1blend(nlam),e2blend(nlam))

  ! Set the number of f values between 0 and fmax, for DHS
  if (fmax .eq. 0e0) then
     ! DHS is turned off by setting fmax to 0.
     nf = 1
  else
     nf = 20
  endif
  allocate(r(ns),nr(ns))
  allocate(f(nf),wf(nf))
  allocate(e1mantle(nlam),e2mantle(nlam))

  ! ----------------------------------------------------------------------
  ! Make copies of arrays and numbers we want to change in this routine
  ! without affecting values in the calling routine
  ! ----------------------------------------------------------------------
  nm          = nm0
  mfrac(1:nm) = mfrac0(1:nm)
  rho(1:nm)   = mat_rho(1:nm)
  
  ! ----------------------------------------------------------------------
  ! Normalize the mass fractions
  ! ----------------------------------------------------------------------
  tot = 0.0_dp
  do im=1,nm
     tot=tot+mfrac(im)
  enddo
  mfrac = mfrac/tot

  ! ----------------------------------------------------------------------
  ! Identify the mantle material and take it out of the main list
  ! ----------------------------------------------------------------------
  i_mantle = 0
  if (trim(mat_loc(nm)).EQ."mantle") then
     ! Remember where the mantle information is
     i_mantle     = nm
     mfrac_mantle = mfrac(nm)
     ! Reduce the number of materials, so that the Bruggeman
     ! rule is only applied to the core materials
     nm = nm-1
  else
     ! These are just to keep the compiler happy
     vfrac_mantle = 0.d0
     mfrac_mantle = 0.d0
     rho_mantle   = 0.d0
  endif

  ! ----------------------------------------------------------------------
  ! Create the size distribution
  ! ----------------------------------------------------------------------
  aminlog = log10(amin)
  amaxlog = log10(amax)
  pow     = -apow
  if (ns.eq.1) then
     ! Just one size
     r(1)  = 10d0**((aminlog+amaxlog)/2d0)
     nr(1) = r(1)**(pow+1d0)
  else
     tot = 0d0
     ! Power-law size distribution
     do is=1,ns
        r(is)=10d0**(aminlog + (amaxlog-aminlog)*real(is-1)/real(ns-1))
        nr(is) = r(is)**(pow+1d0)
        if (r(is).lt.amin .or. r(is).gt.amax) nr(is) = 0.0_dp
        tot=tot+nr(is)*r(is)**3 ! for volume normalization
     enddo
     ! normalize the grain numbers so that the total volume is 1
     do is=1,ns
        nr(is)=1.*nr(is)/tot
     enddo
  endif

  ! ----------------------------------------------------------------------
  ! Copy the refractory index data for all materials into local arrays
  ! ----------------------------------------------------------------------
  ! We do this, because during mixing we overwrite some of the values
  do im=1,nm
     e1(im,1:nlam)    = mat_e1(im,1:nlam)
     e2(im,1:nlam)    = mat_e2(im,1:nlam)
  enddo
  if (i_mantle.gt.0) then
     e1mantle(1:nlam) = mat_e1(i_mantle,1:nlam)
     e2mantle(1:nlam) = mat_e2(i_mantle,1:nlam)
     rho_mantle       = rho(i_mantle)
  endif

  ! ----------------------------------------------------------------------
  ! Turn the mass fractions into volume fractions and compute rho_core
  ! ----------------------------------------------------------------------
  mtot = 0.d0
  vtot = 0.d0
  do im = 1,nm
     mtot = mtot + mfrac(im)
     vfrac(im) = mfrac(im)/rho(im)
     vtot      = vtot+vfrac(im)
  enddo
  rho_core  = mtot/vtot  ! No porosity included yet, will be done below
  ! Normalize the volume fractions of the core
  vfrac(1:nm) = vfrac(1:nm)/vtot

  ! ----------------------------------------------------------------------
  ! Add vacuum as an extra material, to model porosity
  ! ----------------------------------------------------------------------
  nm = nm+1
  e1(nm,1:nlam)     = 1.0_dp
  e2(nm,1:nlam)     = 0.0_dp
  rho(nm)           = 0.0_dp
  vfrac(nm)         = p_c
  vfrac(1:nm-1)     = vfrac(1:nm-1)*(1.0_dp-p_c) ! renormalize to 1
  rho_core          = rho_core * (1.d0-p_c)

  ! ----------------------------------------------------------------------
  ! A this point, the core is properly normalized, including the
  ! porosity. Note that the mantle is no longer in the normalization
  ! mix at this point.
  ! We now compute the average density.
  ! ----------------------------------------------------------------------

  if (i_mantle.eq.0) then
     ! No mantle, rho is already correct
     rho_av = rho_core
  else
     ! Compute the valume fraction of the mantle, from the mass fractions.
     !
     ! The following is the solution of these two equations:
     !
     ! 1. rho_av  = rho_c*(1-vfrac_m) + rho_m*vfrac_m
     ! 2. mfrac_m = rho_m*vfrac_m/rho
     !
     ! Given mfrac_m, rho_c, and rho_m, we solve for rho and vfrac_m
     !
     if (p_m.gt.0d0) then
        ! reduce the density of the mantle material because it is porous as well.
        rho_mantle = rho_mantle * (1.d0-p_m)
     endif
     rho_av        = rho_core / (1.d0 + mfrac_mantle*(rho_core/rho_mantle-1.d0))
     vfrac_mantle  = mfrac_mantle * rho_av / rho_mantle
     if (verbose) then
        write(*,*) 'Volume fractions: vfrac(1:i_mantle-1)*(1.d0-vfrac_mantle),vfrac_mantle'
        write(*,*) 'Total: ',sum(vfrac(1:i_mantle-1)*(1.d0-vfrac_mantle))+vfrac_mantle
        write(*,*) 'Core and Mantle porosities ',p_c,p_m
     endif
  endif

  if (verbose) then
     write(*,'("Average bulk density = ",f6.3, " g/cm3")') rho_av
  endif

  ! ----------------------------------------------------------------------
  ! Now do the mixing, for all wavelengths
  ! ----------------------------------------------------------------------

  ! Loop over wavelengths
  do il=1,nlam
     if (nm0.eq.1 .and. p_c.eq.0) then
        ! Solid core, single material, nothing to blend for the core
     else
        ! Blend the core materials
        do im=1,nm
           epsj(im) = ((dcmplx(e1(im,il),e2(im,il))))**2
        enddo
        if (CLBlend) then
           ! Charlène Lefévre's new blender
           call brugg(vfrac,nm,epsj,eps_eff)
        else
           ! The old Blender from OpacityTool.
           ! Never fails, but might not be converged
           ! ASKMICHIEL
           call Blender(vfrac,nm,epsj,eps_eff)
        endif
        e1blend(il) = dreal(cdsqrt(eps_eff))
        e2blend(il) = dimag(cdsqrt(eps_eff))
     endif
     if (i_mantle.gt.0) then
        ! We do have a mantle to add
        if (p_m.gt.0d0) then
           ! The mantle is porous
           vfrac(1) = 1.d0-p_m
           vfrac(2) = p_m
           epsj(1)  = (dcmplx(e1mantle(il),e2mantle(il)))**2
           epsj(2)  = (dcmplx(1.d0,0.d0))**2
           if (CLBlend) then
              call brugg(vfrac,2,epsj,eps_eff)
           else
              call Blender(vfrac,2,epsj,eps_eff)
           endif
           e1mantle(il) = dreal(cdsqrt(eps_eff))
           e2mantle(il) = dimag(cdsqrt(eps_eff))
        endif
        ! Now we have the mantle material ready - put it on the core
        call maxgarn_2compo(e1blend(il),e2blend(il),e1mantle(il),e2mantle(il),vfrac_mantle, &
             e1mg,e2mg)
        e1blend(il) = e1mg
        e2blend(il) = e2mg
     endif
  enddo ! end of wavelength loop over il

  ! ----------------------------------------------------------------------
  ! We are done with all the blending, so from now on we have ony one material
  ! ----------------------------------------------------------------------
  nm = 1

  ! ----------------------------------------------------------------------
  ! Write the derived n and k to a file
  ! ----------------------------------------------------------------------
  if (blendonly) then
     write(*,'("Writing the blended n and k to blended.lnk, and exiting")')
     call remove_file_if_exists('blended.lnk')
     open(unit=20,file='blended.lnk')
     write(20,'(i5 f5.2)') nlam,rho_av
     do ilam=1,nlam
        write(20,'(3e15.4)') lam(ilam),e1blend(ilam),e2blend(ilam)
     enddo
     close(unit=20)
     stop
  endif

  ! ----------------------------------------------------------------------
  ! Initialize the scattering matrix elements
  ! ----------------------------------------------------------------------
  do il=1,nlam
     do j=1,n_ang
        f11(il,j) = 0d0
        f12(il,j) = 0d0
        f22(il,j) = 0d0
        f33(il,j) = 0d0
        f34(il,j) = 0d0
        f44(il,j) = 0d0
     enddo
  enddo

  ! ----------------------------------------------------------------------
  ! Check how we are going to average over hollow sphere components
  ! ----------------------------------------------------------------------
  if (nf.gt.1 .and. fmax.gt.0.01e0) then
     ! Get the weights for Gauss-Legendre integration
     call gauleg2(0.01e0,fmax,f(1:nf),wf(1:nf),nf)
  else if (fmax.eq.0e0) then
     ! Just a compact sphere, weight is 1
     f(1:nf)  = 0d0
     wf(1:nf) = 1d0/real(nf)
  else
     ! Just one fixed volume fraction of hollow
     f(1)  = fmax
     wf(1) = 1d0
  endif

  ndone = 0
  if (progress) call tellertje(ndone,nlam)
  ! ----------------------------------------------------------------------
  ! Start the main loop over all wavelengths
  ! ----------------------------------------------------------------------
  !$OMP parallel do if (.not. split)                                   &
  !$OMP default(none)                                                  &
  !$OMP shared(f11,f12,f22,f33,f34,f44)                                &
  !$OMP shared(r,lam,nlam,e1blend,e2blend,p,nr,meth,nf,ns,rho_av,wf,f) &
  !$OMP shared(split,progress,ndone)                                   &
  !$OMP private(r1,is,if,rcore,rad)                                    &
  !$OMP private(csca,cabs,cext,mass,vol,theta,mu)                      &
  !$OMP private(cemie,csmie,e1mie,e2mie,rmie,lmie,qabs,qsca,qext,gqsc) &
  !$OMP private(cext_ff,cabs_ff,csca_ff,err,spheres,toolarge)          &
  !$OMP private(m1,m2,d21,s21,m,wvno,min)                              &
  !$OMP private(Mief11,Mief12,Mief22,Mief33,Mief34,Mief44)             &
  !$OMP private(tot,tot2)
  do ilam = 1,nlam
     !$OMP critical
     ndone = ndone+1
     if (progress) call tellertje(ndone,nlam)
     !$OMP end critical
     csca     = 0d0
     cabs     = 0d0
     cext     = 0d0
     mass     = 0d0
     vol      = 0d0

     do j=1,n_ang/2  ! FIXME could be outside of the loop!
        theta=(real(j)-0.5)/real(n_ang/2)*pi/2d0
        mu(j)=cos(theta)
     enddo

     ! ----------------------------------------------------------------------
     ! Start the main loop over all particle sizes
     ! ----------------------------------------------------------------------
     do is=1,ns
        r1       = r(is)
        err      = 0
        spheres  = 0
        toolarge = 0
        ! ----------------------------------------------------------------------
        ! Start the loop over the DHS f factors
        ! ----------------------------------------------------------------------
        cext_ff = 0.d0; cabs_ff = 0.d0; csca_ff = 0.d0
        min = dcmplx(1d0,0d0)
        do if=1,nf
           rad  = r1 / (1d0-f(if))**(1d0/3d0)
           m    = dcmplx(e1blend(ilam),-e2blend(ilam))
           wvno = 2d0*pi / lam(ilam)

           if (f(if) .eq. 0d0) then
              ! solid sphere
              spheres = 1
           else if (r1*wvno.gt.10000d0) then
              ! Sphere is too large
              toolarge = 1
           else if (meth(1:3) .eq. 'DHS') then
              rcore = rad*f(if)**(1d0/3d0)
              call DMiLay(rcore, rad, wvno, m, min, mu, &
                   n_ang/2, qext, qsca, qabs, gqsc, &
                   m1, m2, s21, d21, n_ang ,err)
           else
              rcore = rad*0.999
              call DMiLay(rcore, rad, wvno, min, m, mu, &
                   n_ang/2, qext, qsca, qabs, gqsc, &
                   m1, m2, s21, d21, n_ang ,err)
           endif
           if (err.eq.1 .or. spheres.eq.1 .or. toolarge.eq.1) then
              rad   = r1
              rcore = rad
              rmie  = rad
              lmie  = lam(ilam)
              e1mie = e1blend(ilam)
              e2mie = e2blend(ilam)
              if (err.eq.1 .or. if.eq.1) then
                 if (rmie/lmie.lt.5000d0) then
                    call MeerhoffMie(rmie,lmie,e1mie,e2mie,csmie,cemie, &
                         Mief11,Mief12,Mief33,Mief34,n_ang)
                 else
                    call MeerhoffMie(rmie,rmie/5000d0,e1mie,e2mie,csmie,cemie, &
                         Mief11,Mief12,Mief33,Mief34,n_ang)
                 endif
              endif

              Mief22 = Mief11
              Mief44 = Mief33

           else
              cemie = qext * pi * rad**2
              csmie = qsca * pi * rad**2
              do j=1,n_ang/2
                 Mief11(j) = (M2(j,1) + M1(j,1))         / csmie/wvno**2*2d0*pi
                 Mief12(j) = (M2(j,1) - M1(j,1))         / csmie/wvno**2*2d0*pi
                 Mief22(j) = (M2(j,1) + M1(j,1))         / csmie/wvno**2*2d0*pi
                 Mief33(j) = (S21(j,1))                  / csmie/wvno**2*2d0*pi
                 Mief34(j) = (-D21(j,1))                 / csmie/wvno**2*2d0*pi
                 Mief44(j) = (S21(j,1))                  / csmie/wvno**2*2d0*pi
                 ! Here we use the assumption that the grid is regular.  An adapted
                 ! grid is not possible if it is not symmetric around pi/2.
                 Mief11(n_ang-j+1) = (M2(j,2) + M1(j,2)) / csmie/wvno**2*2d0*pi
                 Mief12(n_ang-j+1) = (M2(j,2) - M1(j,2)) / csmie/wvno**2*2d0*pi
                 Mief22(n_ang-j+1) = (M2(j,2) + M1(j,2)) / csmie/wvno**2*2d0*pi
                 Mief33(n_ang-j+1) = (S21(j,2))          / csmie/wvno**2*2d0*pi
                 Mief34(n_ang-j+1) = (-D21(j,2))         / csmie/wvno**2*2d0*pi
                 Mief44(n_ang-j+1) = (S21(j,2))          / csmie/wvno**2*2d0*pi
              enddo
           endif

           ! make sure the scattering matrix is properly normalized by adjusting the forward peak.
           ! ASKMICHIEL: You are only adjusting the one point a 0 degrees.
           ! So that means that the basic normalization is already in place at this point,
           ! and you just want to be sure that things are right on your specific grid?
           tot  = 0d0
           tot2 = 0d0
           do j=1,n_ang
              ! This integration assumes that the grid is regular (linear)
              tot  = tot +  Mief11(j)*sin(pi*(real(j)-0.5)/real(n_ang))
              tot2 = tot2 + sin(pi*(real(j)-0.5)/real(n_ang))
           enddo
           Mief11(1) = Mief11(1) + (tot2-tot)/sin(pi*(0.5)/real(n_ang))
           if (Mief11(1).lt.0d0) Mief11(1) = 0d0

           ! Add this contribution with the proper weights
           do j=1,n_ang
              f11(ilam,j) = f11(ilam,j) + wf(if)*nr(is)*Mief11(j)*csmie
              f12(ilam,j) = f12(ilam,j) + wf(if)*nr(is)*Mief12(j)*csmie
              f22(ilam,j) = f22(ilam,j) + wf(if)*nr(is)*Mief22(j)*csmie
              f33(ilam,j) = f33(ilam,j) + wf(if)*nr(is)*Mief33(j)*csmie
              f34(ilam,j) = f34(ilam,j) + wf(if)*nr(is)*Mief34(j)*csmie
              f44(ilam,j) = f44(ilam,j) + wf(if)*nr(is)*Mief44(j)*csmie
           enddo
           cext_ff = cext_ff + wf(if)*nr(is)*cemie
           csca_ff = csca_ff + wf(if)*nr(is)*csmie
           cabs_ff = cabs_ff + wf(if)*nr(is)*(cemie-csmie)
           mass = mass + wf(if)*nr(is)*rho_av*4d0*pi*r1**3/3d0
           vol  = vol  + wf(if)*nr(is)*4d0*pi*r1**3/3d0
        enddo    ! end loop "nf" over form factors

        ! Add the contribution of the current grains size to the overall sum
        cext = cext + cext_ff
        cabs = cabs + cabs_ff
        csca = csca + csca_ff

     enddo   ! end loop "is" over grain sizes

     ! ----------------------------------------------------------------------
     ! Set the cross sections
     ! ----------------------------------------------------------------------
     if (ilam .eq.1) p%rho  = mass/vol
     p%Kext(ilam) = 1d4 * cext / mass
     p%Kabs(ilam) = 1d4 * cabs / mass
     p%Ksca(ilam) = 1d4 * csca / mass

     ! ----------------------------------------------------------------------
     ! Set the elements of the scattering matrix
     ! ----------------------------------------------------------------------
     p%F(ilam)%F11(1:180) = f11(ilam,1:180)/csca
     p%F(ilam)%F12(1:180) = f12(ilam,1:180)/csca
     p%F(ilam)%F22(1:180) = f22(ilam,1:180)/csca
     p%F(ilam)%F33(1:180) = f33(ilam,1:180)/csca
     p%F(ilam)%F34(1:180) = f34(ilam,1:180)/csca
     p%F(ilam)%F44(1:180) = f44(ilam,1:180)/csca

     ! ----------------------------------------------------------------------
     ! Average ofver angles to compute asymmetry factor g
     ! ----------------------------------------------------------------------
     ! FIXME: a regular angular grid is assumed for this computation
     tot = 0.0_dp
     p%g(ilam) = 0.d0  !FIXME we never did this, why did this not cause problems???
     do i=1,180
        p%g(ilam) = p%g(ilam) + p%F(ilam)%F11(i)*cos(pi*(real(i)-0.5)/180d0) &
             *sin(pi*(real(i)-0.5)/180d0)
        tot = tot + p%F(ilam)%F11(i)*sin(pi*(real(i)-0.5)/180d0)
     enddo
     p%g(ilam) = p%g(ilam)/tot

     ! Check the normalization of F11
     tot  = 0d0; tot2 = 0d0
     do j=1,n_ang
        !                  F11                      (sin phi)                      d phi        2pi
        tot  = tot +  p%F(ilam)%F11(j) * (sin(pi*(real(j)-0.5)/real(n_ang))) * (pi/180.d0) * (2.d0*pi)
        tot2 = tot2 + sin(pi*(real(j)-0.5)/real(n_ang))*pi/180.d0*2.d0*pi
     enddo
 !    print *,tot,tot2,tot/4/pi,p%F(ilam)%F11(1),p%F(ilam)%F11(2),p%F(ilam)%F11(3)

     ! do it again, using mu integration
     tot  = 0d0; tot2 = 0d0
     do j=1,n_ang-1
        !                  F11                      (sin phi)                      d phi        2pi
        tot  = tot +  2.d0*pi*p%F(ilam)%F11(j) * (cos(pi*(real(j)/real(n_ang))) - cos(pi*(real(j-1)/real(n_ang))))
     enddo
!     print *,'mu integration gives: ',tot,tot/4.d0/pi
     
  enddo   ! end loop ilam over wavelength
  !$OMP end parallel DO
  if (progress) write(*,*)
  
  deallocate(e1,e2)
  deallocate(e1mantle,e2mantle)

  deallocate(Mief11,Mief12,Mief22,Mief33,Mief34,Mief44)
  deallocate(mu,M1,M2,S21,D21)

  deallocate(vfrac)
  deallocate(epsj)
  deallocate(f11,f12,f22,f33,f34,f44)

  deallocate(r,nr)
  deallocate(f,wf)

  return
end subroutine ComputePart

! ----------------------------------------------------------------------
! ----------------------------------------------------------------------
! ----------------------------------------------------------------------
!
!                          Effective medium routines
!
! ----------------------------------------------------------------------
! ----------------------------------------------------------------------
! ----------------------------------------------------------------------

subroutine brugg(f, nm, e, epsavg)

  !***********************************************************************
  !  This subroutine calculates the average dielectric function.
  !  It uses a generalized version of the Bruggeman dielectric mixing function.
  !  This routine was first coded by:
  !
  !  July 2000 : Benjamin T. Johnson, Atmospheric and Oceanic Sciences Dept.
  !              University of Wisconsin - Madison
  !              jbenjam@aos.wisc.edu
  !
  !  Feb. 2002 : Modifications by Michael A. Walters, SSEC, UW-Madison
  !              walters@rain.aos.wisc.edu
  !              Rewritten to use complex variables in input/output.
  !              Converted to Fortran 90
  !
  !  Jul 2018  : Modified by C. Lefèvre, IRAM, France
  !              lefevre@iram.fr
  !              Generalization to nm grain material and rewritten to accept complex arrays
  !              Integrated to SIGMA code to compute dust properties
  !
  !  Jul 2019   : Generalization to n components by Michiel Min and Charlène Lefèvre
  !  The subroutine will:
  !  1. Accept the parameters f, eps, nm
  !
  !  2. Calculate the average or mixed dielectric constant (epsavg), and
  !     return it back to the calling program.
  !
  !  Variable and parameter descriptions:
  !     nm      = number of grain materials
  !     f       = volume fraction of each component
  !     e       = dielectric constant of each component
  !     epsavg  = averaged dielectric constant, the return values
  !**********************************************************************

  use Defs
  implicit none
  integer           :: nm, i, k, l, m
  integer           :: j(nm+1)
  real(kind =dp)    :: f(nm)
  COMPLEX (kind=dp) :: e(nm)
  COMPLEX (kind=dp) :: epsavg
  COMPLEX (kind=dp) :: c(nm+1)
  COMPLEX (kind=dp) :: prod
  COMPLEX (kind=dp) :: x(nm)
  COMPLEX (kind=dp) :: roots(nm)
  COMPLEX (kind=dp) :: total !to check the result of Bruggeman rule
  logical           :: polish
  polish = .false.

  c = 0d0
  do i=1,nm
     x    = -e/2d0
     x(i) = e(i)

     c(nm+1) = c(nm+1)+f(i)
     do k=1,nm
        do l=1,k
           j(l) = l
        enddo
        j(k+1) = nm+1
1       continue
        prod = 1.0_dp
        do l=1,k
           prod = prod*x(j(l))
        enddo
        c(nm-k+1) = c(nm-k+1) + f(i)*prod*(-1.0_dp)**k
        do l=1,k
           if((j(l)+1) .lt. j(l+1)) then
              j(l) = j(l)+1
              do m=1,l-1
                 j(m) = m
              enddo
              goto 1
           endif
        enddo
        continue
     enddo
  enddo

  call zroots(c,nm,roots,polish)
  do i=1,nm
     if(real(roots(i)).gt.0d0 .and. dimag(roots(i)).gt.0d0) THEN
        epsavg=roots(i)
     else if (roots(i).eq.roots(1) .and. dimag(roots(i)).lt.0d0) then
        write(*,*) "ERROR: Bruggeman rule did not converge: "
        write(*,*) "       no positive solution for effective refractive index: "
        write(*,*) roots(i)
        write(*,*) "Please try with a restricted wavelength range"
        write(*,FMT='(a,f6.1,a,f6.1)') "currently lmin = ", lam(1), ", lmax = ", lam(nlam)
        stop
     endif
  enddo

  total = 0.0_dp
  do i = 1,nm
     total = total + (e(i)-epsavg)/(e(i)+2.0_dp*epsavg)*f(i)
  enddo

  if (abs(dreal(total)).gt.1e-15_dp.or.abs(dimag(total)).gt.1e-15_dp) THEN
     write(*,*) "ERROR Bruggeman rule did not converge"
     write(*,*) total
     stop
  endif
  return
end subroutine brugg

subroutine blender(abun,nm,e_in,e_out)
  ! ----------------------------------------------------------------------
  ! This is the original blender routine used in OpacityTool
  ! ----------------------------------------------------------------------
  implicit none
  integer, parameter :: dp = selected_real_kind(P=15)
  integer            :: nm,j,iter
  real (kind=dp)     ::abun(nm)
  complex (kind=dp)  :: e_in(nm),e_out
  complex (kind=dp)  :: mm,m(nm),me,sum

  mm = dcmplx(1d0,0d0)
  m  = e_in
  do iter=1,100
     sum = 0d0
     do j=1,nm
        sum = sum + ((m(j)**2-mm**2)/(m(j)**2+2d0*mm**2))*abun(j)
     enddo
     me = (2d0*sum+1d0)/(1d0-sum)
     me = mm*cdsqrt(me)
     mm = me
  enddo
  e_out = me
end subroutine blender

subroutine zroots(a,m,roots,polish)
  ! Root finder for complex polynomials.  From Numerical Recipes.
  use Defs
  integer          :: m,MAXM
  real(kind=dp)    :: EPS
  complex(kind=dp) :: a(m+1),roots(m)
  logical          :: polish
  parameter (EPS=1.e-15,MAXM=101) ! A small number and maximum anticipated value of m+1.
  ! USES laguer

  ! Given the degree m and the complex coefficients a(1:m+1) of the
  ! polynomial m+1 a(i)xi−1, i=1 this routine successively calls
  ! laguer and finds all m complex roots. The logical variable polish
  ! should be input as .true. if polishing (also by Laguerre’s method)
  ! is desired, .false. if the roots will be subsequently polished by
  ! other means.
  integer j,jj,its
  COMPLEX(kind=dp) ad(MAXM),x,b,c
  do j=1,m+1
     ad(j) = a(j)
  enddo
  ! Copy of coefficients for successive deflation.
  do j=m,1,-1
     ! Loop over each root to be found.
     ! Start at zero to favor convergence to smallest remaining root.
     x = cmplx(0.,0.)
     call laguer(ad,j,x,its) !Find the root.
     if(abs(aimag(x)).le.2.*EPS**2*abs(real(x))) x=cmplx(real(x),0.)
     roots(j) = x
     b        = ad(j+1)
     do jj=j,1,-1
        c      = ad(jj)
        ad(jj) = b
        b      = x*b+c
     enddo
  enddo

  if (polish) then
     do j=1,m
        ! Polish the roots using the undeflated coefficients.
        call laguer(a,m,roots(j),its)
     enddo
  endif
  ! Adapted for SIGMA by Charlene Lefevre - only root with positive real part is kept
  ! FIXME: this assumes that there is only one, and we keep the last one...
  do j=1,m
     if (dreal(roots(j)).gt.0.0_dp) roots(1)=roots(j)
  enddo
  return
end subroutine zroots

subroutine laguer(a,m,x,its)
  ! Root finder for polynomial with complex coefficients
  ! From Numerical Recipes
  use Defs
  integer m,its,MAXIT,MR,MT
  real (kind=dp) EPSS
  COMPLEX (kind=dp) a(m+1),x
  parameter (EPSS=1.e-15_dp,MR=8,MT=10,MAXIT=MT*MR)
  ! Given the degree m and the complex coefficients a(1:m+1) of the
  ! polynomial and given a complex value x, this routine improves x by
  ! Laguerre’s method until it con- verges, within the achievable
  ! roundoff limit, to a root of the given polynomial. The number of
  ! iterations taken is returned as its.  Parameters: EPSS is the
  ! estimated fractional roundoff error. We try to break (rare) limit
  ! cycles with MR different fractional values, once every MT steps,
  ! for MAXIT total allowed iterations.
  integer iter,j
  real (kind=dp) abx,abp,abm,err,frac(MR)
  COMPLEX (kind=dp) dx,x1,b,d,f,g,h,sq,gp,gm,g2
  SAVE frac
  DATA frac /.5,.25,.75,.13,.38,.62,.88,1./
  do iter=1,MAXIT
     its = iter
     b   = a(m+1)
     err = abs(b)
     d   = cmplx(0.,0.)
     f   = cmplx(0.,0.)
     abx = abs(x)
     do j=m,1,-1
        f = x*f+d
        d = x*d+b
        b = x*b+a(j)
        err = abs(b) + abx*err
     enddo
     err = EPSS*err
     if(abs(b).le.err) then
        return
     else
        ! Fractions used to break a limit cycle. Loop over iterations up to allowed maximum.
        ! Efficient computation of the polynomial and its first two derivatives.
        ! Estimate of roundoff error in evaluating polynomial. We are on the root.
        ! The generic case: use Laguerre’s formula.
        g   = d/b
        g2  = g*g
        h   = g2-2.*f/b
        sq  = sqrt((m-1)*(m*h-g2))
        gp  = g+sq
        gm  = g-sq
        abp = abs(gp)
        abm = abs(gm)
        if(abp.lt.abm) gp = gm
        if (max(abp,abm).gt.0.) then
           dx = m/gp
        else
           dx = exp(cmplx(log(1.+abx),float(iter)))
        endif
     endif
     x1 = x-dx
     if(x.eq.x1)return
     if (mod(iter,MT).ne.0) then
        x = x1
     else
        x = x - dx*frac(iter/MT)
     endif
  enddo
  ! pause "too many iterations in laguer"
end subroutine laguer

subroutine maxgarn_2compo(e1in,e2in,e1mantle_in,e2mantle_in,abun,e1out,e2out)
  ! 2 component Maxwell-Garnet mixing
  use Defs
  implicit none
  real (kind=dp) e1in,e2in,e1out,e2out,e1mantle_in,e2mantle_in,abun !up, down
  complex (kind=dp) m1,m2,me

  ! Maxwell Garnet Rules is used
  m2 = dcmplx(e1mantle_in,e2mantle_in) ! mantle = coating material = matrix
  m1 = dcmplx(e1in,e2in)               ! inner core is the "inclusion"
  ! For the abundance of the core, we have to use (1-abun_mantle)
  me = m2**2*( (2d0*m2**2+m1**2-2d0*(1d0-abun)*(m2**2-m1**2)) &
       /(2d0*m2**2+m1**2+(1d0-abun)*(m2**2-m1**2) ) )
  ! from Mukai & Kraetschmer 1986: Optical Constants of the Mixture of Ices
  ! Earth, Moon and Planets, Volume 36, Issue 2, pp.145-155

  me    = cdsqrt(me)
  e1out = dreal(me)
  e2out = dimag(me)

  return
end subroutine maxgarn_2compo

subroutine gauleg2(x1,x2,x,w,n)
  ! Gauss Legendre integration.  From Numerical Recipes
  use Defs
  integer n
  real (kind=dp) x1,x2,x(n),w(n)
  real (kind=dp) EPS
  parameter (EPS=3.d-14)
  integer i,j,m
  real (kind=dp) p1,p2,p3,pp,xl,xm,z,z1
  m  = (n+1)/2
  xm = 0.5d0*(x2+x1)
  xl = 0.5d0*(x2-x1)
  do i=1,m
     z = cos(pi*(i-.25d0)/(n+.5d0))
1    continue
     p1 = 1.d0
     p2 = 0.d0
     do  j=1,n
        p3 = p2
        p2 = p1
        p1 = ((2.d0*j-1.d0)*z*p2-(j-1.d0)*p3)/j
     enddo
     pp = n*(z*p1-p2)/(z*z-1.d0)
     z1 = z
     z = z1-p1/pp
     if (abs(z-z1).gt.EPS) goto 1
     x(i)     = xm-xl*z
     x(n+1-i) = xm+xl*z
     w(i)     = 2.d0*xl/((1.d0-z*z)*pp*pp)
     w(n+1-i) = w(i)
  enddo
  return
end subroutine gauleg2

! ----------------------------------------------------------------------
! ----------------------------------------------------------------------
! ----------------------------------------------------------------------
!
!                      Helper functions and subroutines
!
! ----------------------------------------------------------------------
! ----------------------------------------------------------------------
! ----------------------------------------------------------------------

function cdlog10(x)
  ! Complex logarith base 10
  implicit none
  complex*16 x, cdlog10
  cdlog10=log(x)/log(10d0)
  return
end function cdlog10

subroutine tellertje(i,n)
  ! ----------------------------------------------------------------------
  ! Show a progress bar on STDOUT
  ! ----------------------------------------------------------------------
  implicit none
  integer :: i,n,f,l,ndots,maxdots=20,mindots=20
  ndots = max(mindots,min(n,maxdots))
  if(i.eq.0) then
     do l=1,ndots
        write(*,'(".",$)')
     enddo
     write(*,*)
  else
     f = int(dble(ndots)*dble(i)/dble(n))
     if(dble(ndots)*real(i-1)/real(n) .lt. real(f) &
          & .and. dble(ndots)*real(i+1)/real(n).GT.real(f)) then
        write(*,'(".",$)')
        call flush(6)
     endif
     if(i.eq.n) write(*,*)
  endif
  return
end subroutine tellertje

subroutine check_for_file (file)
  ! Throw an error if FILE does not exist
  character*(*) file
  logical file_exists
  inquire (file=trim(file),exist=file_exists)
  if (.not. file_exists) then
     write(*,'("ERROR: File ",A, " does not exist")') trim(file)
     stop
  endif
end subroutine check_for_file

subroutine remove_file_if_exists (file)
  ! Remove FILE if it already exists
  character*(*) file
  logical file_exists
  inquire (file=trim(file),exist=file_exists)
  if (file_exists) then
     open(unit=90,file=trim(file))
     close(unit=90,status='delete')
  endif
end subroutine remove_file_if_exists

function make_file_path(directory,file)
  !
  ! Concatenate directory and file, with sanity checks and fixes
  !
  character*(*) :: directory,file
  character*500 :: dir,make_file_path
  integer l
  dir = trim(directory)
  l = len_trim(dir)
  if (l.gt.1) then
     do while ((l .gt. 1) .and. dir(l:l) .eq. '/')
        dir(l:l) = ' '
        l = len_trim(dir)
     enddo
  endif
  if (dir .eq. '') then
     make_file_path = file
  else
     make_file_path = trim(trim(dir) // '/' // file)
  endif
  return
end function make_file_path

subroutine make_directory(dir)
  !
  ! Check if directory exists.  If not, create it.
  !
  character*(*) dir
  logical dir_e
  inquire(file=dir, exist=dir_e)
  if (.not. dir_e) then
     print *,'Creating directory: ',trim(dir)
     call system('mkdir -p '//trim(dir))
  endif
end subroutine make_directory

! ----------------------------------------------------------------------
! ----------------------------------------------------------------------
! ----------------------------------------------------------------------
!
!         Functions to help with the extended command line syntax, where
!                we allow a single switch to have several values.
!
! ----------------------------------------------------------------------
! ----------------------------------------------------------------------
! ----------------------------------------------------------------------

function arg_is_value(i)
  ! Check if command line argument is a value. That means the arg is
  ! there, and it is not a switch (starting with -)
  ! Problems: This still gets it wrong for -.5
  integer i
  logical arg_is_value
  character*2 :: value
  call getarg(i,value)
  if ((value.eq.'') .or. (len_trim(value).ge.2) .and. &
       ((value(1:1).eq.'-') .and. ((value(2:2).lt.'0') .or. value(2:2).gt.'9'))) then
     arg_is_value = .false.
  else
     arg_is_value = .true.
  endif
end function arg_is_value

function arg_is_number(i)
  ! Check if command line arg i is a number, i.e.
  ! starts with [0-9] or .[0-9] or -[0-9]
  ! Problems:
  ! - We do not check for -.5
  integer i
  logical arg_is_number
  character*2 :: value
  call getarg(i,value)
  arg_is_number = .false.
  if (len_trim(value).gt.0) then
     if (len_trim(value).eq.1) then
        if (value(1:1).ge.'0' .and. value(1:1).le.'9') then
           arg_is_number = .true.
        endif
     else
        if ((value(1:1).ge.'0' .and. value(1:1).le.'9') .or. &
             (((value(1:1).eq.'.') .or. value(1:1).eq."-") .and. &
             ((value(2:2).ge.'0' .and. value(2:2).le.'9')))) then
           arg_is_number = .true.
        endif
     endif
  endif
end function arg_is_number

! ----------------------------------------------------------------------
! ----------------------------------------------------------------------
! ----------------------------------------------------------------------
!
!           Reading files and checking for some consistency
!
! ----------------------------------------------------------------------
! ----------------------------------------------------------------------
! ----------------------------------------------------------------------

function count_words (string)
  !
  ! Count the whitespace-separated words in STRING. The separator can
  ! be an arbitrary number of space, tab, or newline characters
  !
  character*(*) :: string
  character*(30) :: s
  integer i, count_words
  logical :: in_space=.false.
  s = trim(string) // ' '
  count_words = 0
  if (.not. (s.eq.'')) then
     count_words = 1
     do i=1,len_trim(s)
        if ((s(i:i).eq.' ') .or. (s(i:i).eq.'	') .or. (s(i:i).eq.NEW_LINE('A'))) then
           if (.not.in_space) then
              in_space = .true.
              count_words = count_words+1
           endif
        else
           if (in_space) then
              in_space = .false.
           end if
        end if
     enddo
  endif
end function count_words

function count_numbers (string)
  !
  ! Count the whitespace-separated words in STRING that look like numbers.
  ! Stop whenever something does not look like a number.
  ! This can be fooled by something like "-e" of "72e34.5", but
  ! otherwise it does a decent job
  !
  implicit none
  character*(*) :: string
  character*(100) :: s
  integer i, count_numbers
  logical :: in_space=.false.
  s = trim(string) // '  '
  count_numbers = 0
  if ((.not. (s.eq.'')) .and. (verify(s(1:1),'0123456789-+eE.').eq.0)) then
     do i=2,len_trim(s)+2
        if ((s(i:i).eq.' ') .or. (s(i:i).eq.'	') .or. (s(i:i).eq.NEW_LINE('A'))) then
           if (.not.in_space) then
              in_space = .true.
              count_numbers = count_numbers+1
           endif
        else
           if (in_space) then
              in_space = .false.
              if (verify(s(i:i),'0123456789-+.').gt.0) goto 1
           endif
           if (verify(s(i:i),'0123456789-+eE.').gt.0) goto 1
        end if
     enddo
  endif
1 continue
end function count_numbers

function count_data_lines (file)
  ! Count the data lines in this file.
  ! We skip over comment lines at the beginning of that file.
  ! The first empty line marks the end of the data.
  character*(*)   :: file
  character*(500) :: line
  integer n,count_data_lines
  open(99,file=trim(file))
1 continue
  read(99,'(A)') line
  line = trim(adjustl(line))
  ! Skip commented lines at the beginning of the file
  if (line(1:1).eq.'#' .or. line(1:1).eq.'!' .or. line(1:1).eq.'*') goto 1
  ! OK, here we have the first real line
  n = 1
2 continue
  read(99,'(A)',end=3) line
  if (line .eq. '') goto 3
  n=n+1
  goto 2
3 continue
  close(99)
  count_data_lines = n
end function count_data_lines

subroutine read_lambda_grid(file)
  ! Read the lambda grid from a file.
  ! The file can start with comment lines (* or ! or #).
  ! First non-comment line needs to have the number of lines as first number
  ! The rest should be lines where lambda is always the first value.
  ! For example, an lnk file would work.
  use Defs
  implicit none
  integer i
  character*(*) file
  character*(500) line
  call check_for_file(file)
  open(99,file=file)
1 continue
  read(99,'(A)') line
  line = trim(adjustl(line))
  ! Skip commented lines at the beginning of the file
  if (line(1:1).eq.'#' .or. line(1:1).eq.'!' .or. line(1:1).eq.'*') goto 1
  read(line,*) nlam
  allocate(lam(nlam))
  do i=1, nlam
     read (99, fmt=* ) lam(i)
     ! Convert to microns
     lam(i) = lam(i)
  end do
  close(99)
end subroutine read_lambda_grid


! ----------------------------------------------------------------------
! ----------------------------------------------------------------------
! ----------------------------------------------------------------------
!
!                        Routines to write output files
!
! ----------------------------------------------------------------------
! ----------------------------------------------------------------------
! ----------------------------------------------------------------------

subroutine write_header (unit,cc,amin,amax,apow,na,lmin,lmax, &
     p_core,p_mantle,fmax,mfrac,nm)
  ! ----------------------------------------------------------------------
  ! Write a header describing the full setup of the calculation
  ! CC is the comment character that should be added in front of each line
  ! ----------------------------------------------------------------------
  use Defs
  implicit none
  integer        :: unit
  integer        :: na,i,nm
  real (kind=dp) :: amin,amax,apow,lmin,lmax,p_core,p_mantle,fmax,mfrac(nm)
  real (kind=dp) :: amean(3)
  character*(*)  :: cc
  cc = trim(cc)
  call plmeans(amin,amax,apow,amean)
  write(unit,'(A,"============================================================================")') cc
  write(unit,'(A," Opacities computed by OpTool          <a^n>=",1p,3e9.2e1)') cc,amean
  write(unit,'(A," Parameters:")') cc
  write(unit,'(A,"   amin [um]=",f11.3," amax [um]=",f10.2,"  na  =",I5,"    apow=",g10.2)') cc,amin, amax, na, apow
  write(unit,'(A,"   lmin [um]=",f11.3," lmax [um]=",f10.2,"  nlam=",I5)') cc,lmin, lmax, nlam
  write(unit,'(A,"   porosity =",f11.3," p_mantle = ",f9.3,"            DHS fmax=",g10.2)') cc,p_core,p_mantle,fmax
  write(unit,'(A," Composition:")') cc

  if (mat_rho(1).eq.0) then
     ! We don't have rho yet, so we don't put it into the header.
     ! This is the situation where we write to the screen.
     write(unit,'(A,"  Where   mfrac  rho   Material")') cc
     write(unit,'(A,"  -----   -----  ----  ----------------------------------------------------")') cc
     do i=1,nm
        write(unit,'(A,"  ",A6,f7.3,f6.2,"  ",A)') cc,mat_loc(i), mfrac(i)/sum(mfrac(1:nm)),mat_rho(i),trim(mat_lnk(i))
     enddo
  else
     write(unit,'(A,"  Where   mfrac  Material")') cc
     write(unit,'(A,"  -----   -----  ----------------------------------------------------")') cc
     do i=1,nm
        write(unit,'(A,"  ",A6,f7.3,"  ",A)') cc,mat_loc(i), mfrac(i)/sum(mfrac(1:nm)),trim(mat_lnk(i))
     enddo
  endif
  write(unit,'(A,"============================================================================")') cc
end subroutine write_header

subroutine write_ascii_file(p,amin,amax,apow,na,lmin,lmax,fmax,p_core,p_mantle,&
     mfrac,nm,label,scatter,for_radmc,progress)
  ! ----------------------------------------------------------------------
  ! Write an ASCII file with opacaties.
  ! The routine will include LABEL in the file name.  With the flag
  ! If the flag SCATTER is set, add the full scattering matrix.
  ! With the flag FOR_RADMC, the extension of the file becomes ".inp"
  ! and the normaliation of the scattering matrix is changed according
  ! to RADMC-3D's convention.
  ! ----------------------------------------------------------------------
  use Defs
  implicit none
  real (kind=dp) :: amin,amax,apow,fmax,p_core,p_mantle,mfrac(nm)
  real (kind=dp) :: lmin,lmax,f
  type(particle) :: p
  integer        :: na,i,ilam,iang,nm
  character*(*)  :: label
  character*(3)  :: ext
  character*(23) :: ml
  logical        :: scatter,for_radmc,progress
  character*500  :: file1,file2
  character*500  :: make_file_path ! Function

  if (for_radmc) then
     ext = 'inp'
     ml = 'Z11 Z12 Z22 Z33 Z34 Z44'
  else
     ext = 'dat'
     ml = 'F11 F12 F22 F33 F34 F44'
  endif

  if (label .eq. '') then
     file1 = "dustkappa"      // '.' // ext
     file2 = "dustkapscatmat" // '.' // ext
  else
     file1 = "dustkappa_"      // trim(label) // '.' // ext
     file2 = "dustkapscatmat_" // trim(label) // '.' // ext
  endif
  file1 = make_file_path(outdir,file1)
  file2 = make_file_path(outdir,file2)

  call remove_file_if_exists(file1)
  call remove_file_if_exists(file2)

  if (.not. scatter) then

     if (progress) write(*,'("Writing dust opacity output to file:  ",A)') trim(file1)
     open(20,file=file1,RECL=100000)
     call write_header(20,'#',amin,amax,apow,na,lmin,lmax, &
          p_core,p_mantle,fmax,mfrac,nm)
     write(20,'("# Output file formatted for RADMC-3D, dustkappa, no scattering matrix")')
     write(20,'("#    iformat")')
     write(20,'("#    nlambda")')
     write(20,'("#    lambda[um]         kabs [cm^2/g]      ksca [cm^2/g]      g_asymmetry")')
     write(20,'("#============================================================================")')
     write(20,*) 3  ! iformat
     write(20,*) nlam  ! number of lambda points
     do ilam=1,nlam
        write(20,'(1p,e19.8,1p,e19.8,1p,e19.8,1p,e19.8)') lam(ilam),p%Kabs(ilam),p%Ksca(ilam),p%g(ilam)
     enddo
     close(20)

  else

     if (progress) write(*,'("Writing full scattering data to file: ",A)') trim(file2)
     open(20,file=file2,RECL=100000)
     call write_header(20,'#',amin,amax,apow,na,lmin,lmax, &
          p_core,p_mantle,fmax,mfrac,nm)
     if (for_radmc) then
        write(20,'("# Output file formatted for RADMC-3D, dustkapscatmat, RADMC normalization")')
     else
        write(20,'("# Standard output file, scattering matrix mormalizattion: <F11>=1/sr")')
     endif
     write(20,'("#    iformat                                     ! format number")')
     write(20,'("#    nlam                                        ! number of wavelengths")')
     write(20,'("#    nang                                        ! number of angles 0-180")')
     write(20,'("#    lam(1)    kabs(1)    ksca(1)    g(1)        ! um, cm^2/g, cm^2/g, none")')
     write(20,'("#    ...")')
     write(20,'("#    lam(nlam) kabs(nlam) ksca(nlam) g(nlam)")')
     write(20,'("#    ang(1)                                      ! ang(1)    must be 0")')
     write(20,'("#    ...")')
     write(20,'("#    ang(nang)                                   ! ang(nang) must be 180")')
     write(20,'("#    ",A,                  "                     ! (ilam=   1,iang=   1)")') ml
     write(20,'("#    ...")')
     write(20,'("#    ",A,                  "                     ! (ilam=   1,iang=nang)")') ml
     write(20,'("#    ",A,                  "                     ! (ilam=   2,iang=   1)")') ml
     write(20,'("#    ... ...")')
     write(20,'("#    ",A,                  "                     ! (ilam=nlam,iang=nang)")') ml
     write(20,'("#============================================================================")')

     if (for_radmc) then
        write(20,*) 1    ! iformat
     else
        write(20,*) 0    ! This is supposed to cause an error when RADMC-3D is reading it
     endif

     write(20,*) nlam ! Number of wavelength points
     write(20,*) 181  ! Number of angular points
     write(20,*)      ! an empty line
     ! The opacities as function of lambda
     do ilam=1,nlam
        write(20,'(1p,e19.8,1p,e19.8,1p,e19.8,1p,e19.8)') lam(ilam),p%Kabs(ilam),p%Ksca(ilam),p%g(ilam)
     enddo
     write(20,*)      ! an empty line
     ! The angular grid
     do iang=0,180
        write(20,'(f8.2)') real(iang)
     enddo
     write(20,*)      ! an empty line
     ! Write the scattering matrix
     do ilam=1,nlam
        ! Set the scaling for the normalization
        if (for_radmc) then
           f = p%ksca(ilam)/(4.d0*pi)
        else
           f = 1.d0
        endif
        do iang=0,180
           ! We have only computed 0-179, but RADMC needs 180 as well
           ! We simply repeat the 179 value
           i = 1 + min(iang,179)
           write(20,'(1p,e19.8,1p,e19.8,1p,e19.8,1p,e19.8,1p,e19.8,1p,e19.8)') &
                p%F(ilam)%F11(i)*f,p%F(ilam)%F12(i)*f,p%F(ilam)%F22(i)*f, &
                p%F(ilam)%F33(i)*f,p%F(ilam)%F34(i)*f,p%F(ilam)%F44(i)*f
        enddo
     enddo
     close(20)
  endif
end subroutine write_ascii_file

#ifdef USE_FITSIO
subroutine write_fits_file(p,amin,amax,apow,na, &
     fmax,p_core,p_mantle, &
     nm,mfrac,fitsfile)
  ! ----------------------------------------------------------------------
  ! Routine to write a FITS file with all the information about
  ! the opacities and scattering properties.
  ! FIXME: Something goes wrong with the keywords
  ! ----------------------------------------------------------------------
  use Defs
  implicit none
  real (kind=dp) :: amin,amax,apow,fmax,p_core,p_mantle
  real (kind=dp) :: mfrac(nm),rho(nm)
  logical blend
  character*6 word
  character*500 fitsfile
  type(particle) p
  integer nm,na,i,j,nm2
  real (kind=dp),allocatable :: array(:,:,:)
  real (kind=dp) :: amean(3)

  integer status,unit,blocksize,bitpix,naxis,naxes(3)
  integer group,fpixel,nelements
  logical simple,extend
  real a0,a1,a2,a3

  call remove_file_if_exists(fitsfile)
  write(*,'("Writing full scattering data to file: ",A)') trim(fitsfile)

  status = 0
  ! Get an unused Logical Unit Number to use to create the FITS file
  call ftgiou(unit,status)
  ! Create the new empty FITS file
  blocksize = 1
  call ftinit(unit,fitsfile,blocksize,status)

  simple=.true.
  extend=.true.
  group=1
  fpixel=1

  bitpix=-64
  naxis=2
  naxes(1)=nlam
  naxes(2)=5
  nelements=naxes(1)*naxes(2)
  allocate(array(nlam,5,1))

  ! Write the required header keywords.
  call ftphpr(unit,simple,bitpix,naxis,naxes,0,1,extend,status)

  ! Write optional keywords to the header

  call ftpkye(unit,'r_min',real(amin),8,'[micron]',status)
  call ftpkye(unit,'r_max',real(amax),8,'[micron]',status)
  call ftpkye(unit,'r_pow',real(apow),8,'',status)
  call ftpkye(unit,'f_max',real(fmax),8,'',status)

  call plmeans(amin,amax,apow,amean)
  a1 = amean(1)
  call ftpkye(unit,'a1',real(a1),8,'[micron]',status)
!  call ftpkye(unit,'density',real(rho_av),8,'[g/cm^3]',status)

  call ftpkye(unit,'porosity',real(p_core),8,'[g/cm^3]',status)
  call ftpkye(unit,'p_mantle',real(p_mantle),8,'[g/cm^3]',status)

  do i=1,nm
     write(word,'("file",i0.2)') i
     call ftpkys(unit,word,trim(mat_lnk(i)),'',status)
  enddo
  do i=1,nm
     write(word,'("frac",i0.2)') i
     call ftpkye(unit,word,real(mfrac(i)),8,'[mass fraction]',status)
  enddo
!  do i=1,nm
!     write(word,'("rho",i0.2)') i
!     call ftpkye(unit,word,real(rho(i)),8,'[g/cm^3]',status)
!  enddo

  call ftpkyj(unit,'n_radii',na,' ',status)
  call ftpkyj(unit,'n_mat',nm,' ',status)

  !  Write the array to the FITS file.

  ! ----------------------------------------------------------------------
  ! HDU 0: opacities
  ! ----------------------------------------------------------------------

  do i=1,nlam
     array(i,1,1)=lam(i)
     array(i,2,1)=p%Kext(i)
     array(i,3,1)=p%Kabs(i)
     array(i,4,1)=p%Ksca(i)
     array(i,5,1)=p%g(i)
  enddo

  call ftpprd(unit,group,fpixel,nelements,array(1:nlam,1:4,1),status)

  deallocate(array)

  ! ----------------------------------------------------------------------
  ! HDU 1: scattering matrix
  ! ----------------------------------------------------------------------
  bitpix=-64
  naxis=3
  naxes(1)=nlam
  naxes(2)=6
  naxes(3)=180
  nelements=naxes(1)*naxes(2)*naxes(3)

  allocate(array(nlam,6,180))

  ! create new hdu
  call ftcrhd(unit, status)

  !  Write the required header keywords.
  call ftphpr(unit,simple,bitpix,naxis,naxes,0,1,extend,status)

  do i=1,nlam
     do j=1,180
        array(i,1,j)=p%F(i)%F11(j)
        array(i,2,j)=p%F(i)%F12(j)
        array(i,3,j)=p%F(i)%F22(j)
        array(i,4,j)=p%F(i)%F33(j)
        array(i,5,j)=p%F(i)%F34(j)
        array(i,6,j)=p%F(i)%F44(j)
     enddo
  enddo

  !  Write the array to the FITS file.
  call ftpprd(unit,group,fpixel,nelements,array,status)

  deallocate(array)

  !  Close the file and free the unit number.
  call ftclos(unit, status)
  call ftfiou(unit, status)

  !  Check for any error, and if so print out error messages
  if (status.gt.0) then
     print*,'error in export to fits file',status
  end if
  return
end subroutine write_fits_file
#endif

subroutine plmeans(a1,a2,p,amean)
  ! Compute the moments of the size disribution f(a) ~ a^(-p)
  ! The results are returned in AMEAN, an array of three: <a>, <a^2>, <a^3>
  implicit none
  integer, parameter     :: dp = selected_real_kind(P=15)
  real (kind=dp) :: a1,a2,p,amean(3),e1,e2,e3
  integer :: n
  do n=1,3
     e1 = 1.d0-p+n
     e2 = 1.d0-p
     e3 = 1.d0/n
     if (abs(e1).lt.1e-3) then
        amean(n) = (  e2     * (log(a2)-log(a1)) / (a2**e2-a1**e2) )**e3
     else
        amean(n) = ( (e2/e1) * (a2**e1-a1**e1)   / (a2**e2-a1**e2) )**e3
     endif
  enddo
end subroutine plmeans

! ----------------------------------------------------------------------
! ----------------------------------------------------------------------
! ----------------------------------------------------------------------
!     Code to compute mean opacities kappa_Rosseland, kappa_Planck
!
!         Adapted from code provided by Kees Dullemond
! ----------------------------------------------------------------------
! ----------------------------------------------------------------------
! ----------------------------------------------------------------------

subroutine mean_opacities(lambda,nlam,kabs,ksca,g,tmin,tmax,ntemp,file)
  ! ----------------------------------------------------------------------
  ! Compute the mean opacities and write them to FILE
  ! ----------------------------------------------------------------------
  implicit none
  integer :: nlam,ntemp
  doubleprecision, external :: bplanck,bplanckdt
  doubleprecision :: lambda(nlam),kabs(nlam),ksca(nlam),g(nlam)
  doubleprecision :: tmin,tmax
  character*(*)   :: file
  doubleprecision, allocatable :: temp(:),kross(:),kplanck(:)
  doubleprecision, allocatable :: nu(:)
  integer :: ilam,itemp
  doubleprecision dumbnu,dumdb
  doubleprecision bnu,dbnudt,dnu
  doubleprecision kaptot_abs,kaptot_scat,kap_p,kap_r

  allocate(temp(ntemp),kross(ntemp),kplanck(ntemp),nu(nlam))
  nu = 2.99792458d10/(lambda*1d-4)  ! includes conversion lambda from micron to cm.
  !
  ! Make the temperature grid
  !
  do itemp=1,ntemp
     temp(itemp) = tmin * (tmax/tmin)**((itemp-1.0)/(ntemp-1.0))
  enddo
  !
  ! Loop over temperatures
  !
  do itemp=1,ntemp
     !
     ! Reset
     !
     kap_p  = 0.d0
     kap_r  = 0.d0
     dumbnu = 0.d0
     dumdb  = 0.d0
     !
     ! Loop over wavelength
     !
     do ilam=1,nlam
        !
        ! Compute the Planck function B_nu(T) and the dB_nu(T)/dT function
        !
        bnu    = bplanck(temp(itemp),nu(ilam))
        dbnudt = bplanckdt(temp(itemp),nu(ilam))
        !
        ! Compute dnu
        !
        if(ilam.eq.1) then
           dnu = 0.5d0*abs(nu(2)-nu(1))
        elseif(ilam.eq.nlam) then
           dnu = 0.5d0*abs(nu(nlam)-nu(nlam-1))
        else
           dnu = 0.5d0*abs(nu(ilam+1)-nu(ilam-1))
        endif
        !
        ! Set the dust opacities. Note that for the scattering I
        ! reduce the scattering opacity by a factor of (1-g), where
        ! g is the scattering anisotropy factor, because according to
        ! Ishimaru (1978) in the diffusion limit non-isotropic scattering
        ! with g is equivalent to reduced isotropic scattering by a factor
        ! of (1-g).
        !
        kaptot_abs  = kabs(ilam)
        kaptot_scat = ksca(ilam) * (1.d0-g(ilam))
        !
        ! Integrate and add to the averages
        !
        kap_p  = kap_p  + kaptot_abs * bnu * dnu
        kap_r  = kap_r  + dbnudt * dnu / ( kaptot_abs + kaptot_scat )
        dumbnu = dumbnu + bnu * dnu
        dumdb  = dumdb  + dbnudt * dnu
     enddo
     !
     ! Compute the mean opacities
     !
     kap_p  = kap_p  / dumbnu     ! Planck mean with local temperature
     kap_r  = dumdb  / kap_r      ! Rosseland mean
     !
     ! Store the result
     !
     kplanck(itemp)      = kap_p
     kross(itemp)        = kap_r
  enddo
  !
  ! Write the output file
  !
  write(*,'("Writing mean opacities to file:       ",A)') trim(file)
  open(unit=1,file=trim(file))
  write(1,*) ntemp
  do itemp=1,ntemp
     write(1,*) temp(itemp),kplanck(itemp),kross(itemp)
  enddo
  close(1)
  !
  ! Deallocate stuff
  !
  deallocate(temp,kross,kplanck,nu)

end subroutine mean_opacities


! ----------------------------------------------------------------------
!                THE BLACKBODY PLANCK FUNCTION B_nu(T)
!
!     This function computes the Blackbody function
!
!                    2 h nu^3 / c^2
!        B_nu(T)  = ------------------    [ erg / cm^2 s ster Hz ]
!                   exp(h nu / kT) - 1
!
!     ARGUMENTS:
!        nu    [Hz]            = Frequency
!        temp  [K]             = Temperature
! ----------------------------------------------------------------------
function bplanck(temp,nu)
  implicit none
  doubleprecision :: temp
  doubleprecision :: nu
  doubleprecision :: bplanck

  if(temp.eq.0.d0) then
     bplanck = 0.d0
     return
  endif

  bplanck = 1.47455d-47 * nu * nu * nu /                   &
            (exp(4.7989d-11 * nu / temp)-1.d0) + 1.d-290

  return
end function bplanck

! ----------------------------------------------------------------------
!           THE TEMPERATURE DERIVATIVE OF PLANCK FUNCTION
!
!      This function computes the temperature derivative of the
!      Blackbody function
!
!         dB_nu(T)     2 h^2 nu^4      exp(h nu / kT)        1
!         --------   = ---------- ------------------------  ---
!            dT          k c^2    [ exp(h nu / kT) - 1 ]^2  T^2
!
!      ARGUMENTS:
!         nu    [Hz]            = Frequency
!         temp  [K]             = Temperature
! ----------------------------------------------------------------------

function bplanckdt(temp,nu)
  implicit none
  doubleprecision :: temp,nu
  doubleprecision :: theexp,bplanckdt

  theexp = exp(4.7989d-11*nu/temp)
  if(theexp.lt.1.d33) then
     bplanckdt = 7.07661334104d-58 * nu**4 * theexp /    &
          ( (theexp-1.d0)**2 * temp**2 ) + 1.d-290
  else
     bplanckdt = 7.07661334104d-58 * nu**4 /             &
          ( theexp * temp**2 ) + 1.d-290
  endif
  return
end function bplanckdt


