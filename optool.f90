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
  write(*,'("-h                        Show this message")')
  write(*,'("===============================================================================")')
end subroutine usage

module Defs
  ! Module with constants and the basic data structures used in the program
  implicit none
  integer, public, parameter     :: dp = selected_real_kind(P=15)
  ! -----------------------------------------------------------------------
  ! Physics and math constants
  ! -----------------------------------------------------------------------
  real (kind=dp), public, parameter :: pi = 3.1415926535897932384_dp
  ! -----------------------------------------------------------------------
  ! Some global switches
  ! -----------------------------------------------------------------------
  logical, public                :: verbose   = .false. ! additional output
  logical, public                :: blendonly = .false. ! only blend materials and write result out
  logical, public                :: CLBlend   = .true.  ! use Charléne Lefévre's new blender
  ! ------------------------------------------------------------------------
  ! Lambda is shared, because multiple routines need it
  ! ------------------------------------------------------------------------
  real (kind=dp), allocatable    :: lam(:)     ! wavelength
  integer                        :: nlam       ! nr of wavelength points
  ! ------------------------------------------------------------------------
  ! Mueller matrix structure, records only non-zero elements of the matrix
  ! ------------------------------------------------------------------------
  type mueller
     real (kind=dp) :: F11(180),F12(180),F22(180),F33(180),F44(180),F34(180)
  end type mueller
  ! ------------------------------------------------------------------------
  ! The particle structure, contains particle and scatteting properties
  ! ------------------------------------------------------------------------
  type particle
     real (kind=dp)              :: rv,rvmin,rvmax ! grain radius with min and max
     real (kind=dp)              :: rho            ! mass density in g.cm^-3
     real (kind=dp), allocatable :: Kabs(:),Ksca(:),Kext(:) ! Opacities
     real (kind=dp), allocatable :: g(:)           ! asymmetry, Henyey Greenstein
     TYPE(MUELLER),  allocatable :: F(:)           ! Mueller matrix elements
  end type particle
end module Defs

program optool
  use Defs
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
  real (kind=dp)  :: tmin,tmax       ! min and max temperature for mean opacities
  integer         :: nt              ! number of temperature steps
  
  integer         :: nm,nm_input     ! nr of grain materials
  
  type(particle)  :: p
  integer         :: i, j
  integer         :: ilam,iang
  character*100   :: tmp
  character*100   :: value
  
  integer,        allocatable  :: number(:)       ! number associated to a component
  character*500,  allocatable  :: location(:)     ! Either 'core' or 'mantle'
  character*500,  allocatable  :: ref_index(:)
  real (kind=dp), allocatable  :: rho(:)          ! specific mass density of material
  real (kind=dp), allocatable  :: mfrac(:)        ! mass fraction of each component
  real (kind=dp), allocatable  :: mfrac_user(:)   ! original values

  logical         :: have_mantle
  logical         :: arg_is_value, arg_is_number  ! functions to test arguments
  character*500   :: fitsfile,meanfile            ! file names for output
  character*50    :: radmclbl   = ""              ! file label for RADMC-3D compatible 
  character*4     :: asciiext   = ".dat"          ! files extension for ASCII output 
  character*50    :: label                        ! for use in file names

  logical         :: split=.false.
  real (kind=dp)  :: asplit,afact,afsub,amaxsplit,aminsplit
  integer         :: nsubgrains = 5,nsub
  integer         :: ia,is

  ! ------------------------------------------------------------------------
  ! Defaults values for parameters and switches
  ! ------------------------------------------------------------------------
  amin          = 0.05       ! micrometer
  amax          = 3000.      ! micrometer
  apow          = 3.50_dp
  na            = 100
  
  lmin          = 0.05_dp    ! micrometers
  lmax          = 10000.0_dp ! micrometers
  nlam          = 300

  write_fits    = .false.
  write_scatter = .false.
  write_mean_kap= .false.
  tmin          = 10d0       ! K
  tmax          = 1d4        ! K
  nt            = 200
  
  fmax          = 0.8_dp     ! volume fraction DHS

  p_core        = 0.0_dp     ! porosity core
  p_mantle      = 0.0_dp     ! porosity mantle

  nm            = 0          ! number of materials - zeor to start with

  ! ------------------------------------------------------------------------
  ! Allocate space for up to 12 different materials
  ! ------------------------------------------------------------------------
  allocate(number(12))
  allocate(location(12))
  allocate(ref_index(12))
  allocate(mfrac(12),mfrac_user(12))
  allocate(rho(12))
  have_mantle = .false.
  ! Initialize rho, because we need the fact that it has not been set
  do i=1,12
     rho(i) = 0.d0
  enddo

  ! ------------------------------------------------------------------------
  ! Process the command line arguments
  ! ------------------------------------------------------------------------

  ! Loop over all command line arguments
  call getarg(1,tmp)
  i = 1
  do while(tmp.ne.' ')
     select case(tmp)

        ! ------------------------------------------------------------------
        ! Definition of a material for the mix
        ! ------------------------------------------------------------------
     case('-c','-m')
        i  = i+1
        nm = nm+1; number(nm) = nm

        ! First value is the material key or refindex file path
        call getarg(i,value); read(value,'(A)') ref_index(nm)

        if (value .eq. '?') then
           call ListBuiltinMaterials()
           stop
        endif
        ! Second value is the volume fraction
        if (.not. arg_is_value(i+1)) then
           print *, "WARNING: 1.0 used for missing mass fraction of material: ",trim(ref_index(nm))
           mfrac(nm) = 1.0
        else
           i = i+1; call getarg(i,value); read(value,*) mfrac(nm)
        endif
        
        if (mfrac(nm).eq.0) then
           print *, "WARNING: Ignoring material with zero mass fraction: ",trim(ref_index(nm))
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
                 location(nm) = 'mantle'
                 have_mantle = .true.
              endif
           else
              ! This is a core material.
              if (have_mantle) then
                 print *,"ERROR: Mantle material must be specified last"; stop
              endif
              location(nm) = 'core'
           endif
        endif

        ! There might be a density, as a third argument
        if (arg_is_value(i+1)) then
           i = i+1; call getarg(i,value); read(value,*) rho(nm)
        endif
        
        ! ------------------------------------------------------------------
        ! Grain size setup
        ! ------------------------------------------------------------------
     case('-a')
        ! ------------------------------------------------------------------
        ! -a expects 2, 3, ro 4 values:  amin amax [na [apow]]
        ! ------------------------------------------------------------------
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
        
        ! ------------------------------------------------------------------
        ! Wavelength setup
        ! ------------------------------------------------------------------
     case('-l')
        ! ------------------------------------------------------------------
        ! -l expects a file name, or 2-3 numbers: lmin lmax [nlam]
        ! ------------------------------------------------------------------
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
        
        ! ------------------------------------------------------------------
        ! Grain geometry
        ! ------------------------------------------------------------------
     case('-p','-porosity','--porosity')
        i = i+1;  call getarg(i,value); read(value,*) p_core
        if (arg_is_value(i+1)) then
           i=i+1; call getarg(i,value); read(value,*) p_mantle
        else
           p_mantle = p_core
        endif
     case('-fmax')
        i = i+1; call getarg(i,value); read(value,*) fmax

        ! ------------------------------------------------------------------
        ! Temperature setup for mean kappas
        ! ------------------------------------------------------------------
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

        ! ------------------------------------------------------------------
        ! Various other switches
        ! ------------------------------------------------------------------
     case('-s','-scatter','-scat')
        write_scatter=.true.
     case('-radmc','-radmc3d')
        asciiext = ".inp"
        if (arg_is_value(i+1)) then
           i=i+1
           call getarg(i,radmclbl)
        endif
     case('-d')
        split = .true.
        if (arg_is_value(i+1)) then
           i=i+1; call getarg(i,value); read(value,*) nsubgrains;
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
     call getarg(i,tmp)
  enddo
  
  ! ------------------------------------------------------------------
  ! Default grain composition if nothing is specified
  ! ------------------------------------------------------------------
  if (nm.eq.0) then
     ! Set default composition will set a default composition here, the DIANA opacities
     write(*,'("No materials specified, using DIANA standard")')
     nm = 2
     ref_index(1) = 'pyr-mg70' ; location(1)  = 'core' ; mfrac(1)     = 0.87d0
     ref_index(2) = 'c-z'      ; location(2)  = 'core' ; mfrac(2)     = 0.1301d0
     p_core       = 0.25d0
  endif
  nm_input = nm
  mfrac_user = mfrac

  if (nm .ge. 10) then
     write(*,*) 'ERROR: Too many materials'
     stop
  endif
    do i = 1, nm
     location(i)   = trim(location(i))
     ref_index(i)  = trim(ref_index(i))
  enddo
  
  ! ------------------------------------------------------------------
  ! Sanity checks
  ! ------------------------------------------------------------------
  
  ! ------------------------------------------------------------------
  ! Write a setup summary to the screen
  ! ------------------------------------------------------------------
  call write_header(6,'',amin,amax,apow,na,lmin,lmax,nlam, &
       p_core,p_mantle,fmax,nm_input,location,mfrac/sum(mfrac(1:nm)),mfrac,ref_index,rho)
 
  meanfile    = "dustkapmean.dat"
  fitsfile    = "dustkappa.fits"

  ! ------------------------------------------------------------------
  ! Make a logarithmic lambda grid, unless read in from file
  ! ------------------------------------------------------------------
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

  ! ----------------------------------------------------------------------
  ! Loop for splitting the output into files by grain size
  ! ----------------------------------------------------------------------
  if (split) then
     nsub = nsubgrains
     if (mod(nsub,2).eq.0) nsub = nsub+1
     afact = (amax/amin)**(1.d0/real(na)) ! factor to next grain size
     afsub = afact**(1.d0/real(nsub-1))     ! Factor to next subgrain size
     do ia=1,na
        asplit    = amin*afact**real(ia-1+0.5)
        aminsplit = asplit*afsub**real(-nsub/2)
        amaxsplit = asplit*afsub**real(+nsub/2)

        write(label,'(I3.3)') ia
        write(fitsfile,'(A,"_",A,".fits")') "dustkappa",trim(label)
        write(*,'("Computing grain size: ",f10.3," micron")') asplit

        call ComputePart(p,aminsplit,amaxsplit,apow,nsub,fmax,p_core,p_mantle, &
             mfrac,location,ref_index,rho,nm)

        if (write_fits) then
           call write_fits_file(p,aminsplit,amaxsplit,apow,nsub, &
                fmax,p_core,p_mantle, &
                nm_input,mfrac,ref_index,fitsfile)
        else
           if (radmclbl .ne. ' ') label = trim(radmclbl) // "_" // label
           call write_ascii_file(p,aminsplit,amaxsplit,apow,nsub,lmin,lmax, &
                fmax,p_core,p_mantle,nm_input, &
                location,mfrac,mfrac_user,ref_index,rho, &
                label,write_scatter,asciiext)
        endif
     enddo
     stop
  endif

  ! ------------------------------------------------------------------
  ! Call the main routine to compute the opacities and scattering matrix
  ! ------------------------------------------------------------------
  call ComputePart(p,amin,amax,apow,na,fmax,p_core,p_mantle,&
       mfrac,location,ref_index,rho,nm)
  ! ------------------------------------------------------------------
  ! Write the output
  ! ------------------------------------------------------------------
  if (write_fits) then
     call write_fits_file(p,amin,amax,apow,nsub, &
          fmax,p_core,p_mantle,&
          nm_input,mfrac,ref_index,fitsfile)
  else
     call write_ascii_file(p,amin,amax,apow,na,lmin,lmax, &
          fmax,p_core,p_mantle,nm_input,&
          location,mfrac,mfrac_user,ref_index,rho, &
          radmclbl,write_scatter,asciiext)
  endif
  ! ------------------------------------------------------------------
  ! Produce and write mean opacities
  ! ------------------------------------------------------------------
  if (write_mean_kap) then 
     call mean_opacities(lam,nlam,p%kabs,p%ksca,p%g,tmin,tmax,nt,meanfile)
  endif

  deallocate(number)
  deallocate(location)
  deallocate(ref_index)
  deallocate(mfrac,mfrac_user)
  deallocate(rho)
  deallocate(p%Kabs)
  deallocate(p%Kext)
  deallocate(p%Ksca)
  deallocate(p%g)
  deallocate(p%F)
  
end program optool

subroutine ComputePart(p,amin,amax,apow,na,fmax,p_c,p_m,mfrac0,loc,ref_index,rho,nm0)
  ! ----------------------------------------------------------------------------
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
  !   loc        Location of material (core or mantle)
  !   ref_index  key of file with refractive index data
  !   nm0        The number of materials, also the size of the arrays mfrac,
  !              loc, ref_index, rho
  ! 
  ! OUTPUT
  !   p          A "particle" structure to return all the information in.
  !              See the module Defs for the definition.
  ! ----------------------------------------------------------------------------
  use Defs
  !use omp_lib
  implicit none
  
  real (kind=dp)                 :: amin,amax,apow   ! min and max grain size, and power law exp
  real (kind=dp)                 :: fmax             ! maximum fraction of vaccum for DHS
  real (kind=dp)                 :: p_c,p_m          ! porosity, core and mantle
  
  integer                        :: nm0,nm,im        ! nr of grain materials in a composite grain
  real (kind=dp)                 :: rho(nm0)         ! specific material dnesities
  real (kind=dp)                 :: mfrac0(nm0)      ! mass fractions, input
  real (kind=dp)                 :: mfrac(nm0)       ! mass fractions, modified
  character*500                  :: loc(nm0)         ! 'core' or 'mantle'
  character*500                  :: ref_index(nm0)   ! Key or file with refractive index data
  real (kind=dp)                 :: vfrac_mantle     ! volume fraction of mantle
  real (kind=dp)                 :: mfrac_mantle     ! mass   fraction of mantle
  
  TYPE (PARTICLE)                :: p
  
  integer                        :: i,j,k
  integer                        :: MAXMAT           ! maximum number of grain material
  integer, parameter             :: n_ang = 180      ! number of angles FIXME let the user set this?
  integer                        :: na               ! Number of grains sizes between amin and amax
  integer                        :: nf,if            ! Number of DHS volume fractions
  integer                        :: ns,is            ! Number of grains sizes
  integer                        :: ilam,il          ! Counter for wavelengths
  integer                        :: i_mantle         ! Index of mantle material, if any
  integer                        :: err,spheres,toolarge 
  
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
  
  real (kind=dp)                 :: maxf
  real (kind=dp)                 :: aminlog, amaxlog
  real (kind=dp)                 :: rad
  real (kind=dp)                 :: r1
  real (kind=dp)                 :: tot, tot2,mtot,vtot
  real (kind=dp)                 :: pow
  real (kind=dp)                 :: mass
  real (kind=dp)                 :: vol
  real (kind=dp)                 :: rho_av,rho_core,rho_mantle
  real (kind=dp)                 :: rcore
  real (kind=dp)                 :: wvno
  
  ! Effective refractive index
  real (kind=dp)                 :: e1blend,        e2blend
  real (kind=dp), allocatable    :: e1(:,:),        e2(:,:)
  real (kind=dp), allocatable    :: e1d(:),         e2d(:)
  real (kind=dp), allocatable    :: e1mantle(:),    e2mantle(:)
  real (kind=dp), allocatable    :: vfrac(:)
  
  COMPLEX (kind=dp), allocatable :: epsj(:)
  COMPLEX (kind=dp)              :: m       ! FIXME rename?
  COMPLEX (kind=dp)              :: eps_eff,eps_eff2
  COMPLEX (kind=dp)              :: min
  
  character (len=3)              :: meth
  character (len=500)            :: mantle

  real (kind=dp) :: cabs_mono,cabs_rgd,cemie_mono,csmie_mono,G,nmono

  nm    = nm0    ! Make copy, to that we can change it
  mfrac = mfrac0 ! Make copy, to that we can change it
  meth = 'DHS'   ! 
  maxf = fmax
  ns   = na      ! number of subgrains to compute
  MAXMAT = nm+1  ! Allocate one more, because vacuum will also be a material
  
  allocate(Mief11(n_ang))
  allocate(Mief12(n_ang))
  allocate(Mief22(n_ang))
  allocate(Mief33(n_ang))
  allocate(Mief34(n_ang))
  allocate(Mief44(n_ang))
  allocate(mu(n_ang))
  allocate(M1(n_ang,2))
  allocate(M2(n_ang,2))
  allocate(S21(n_ang,2))
  allocate(D21(n_ang,2))
  
  allocate(vfrac(MAXMAT))
  allocate(epsj(MAXMAT))
  allocate(f11(nlam,n_ang))
  allocate(f12(nlam,n_ang))
  allocate(f22(nlam,n_ang))
  allocate(f33(nlam,n_ang))
  allocate(f34(nlam,n_ang))
  allocate(f44(nlam,n_ang))
  
  allocate(e1(MAXMAT,nlam))
  allocate(e2(MAXMAT,nlam))

  ! Set the number of f values between 0 and fmax, for DHS
  if (maxf .eq. 0e0) then
     ! DHS is turned off by setting fmax to 0.
     nf = 1
  else
     nf = 20
  endif
  allocate(r(ns))
  allocate(nr(ns))
  allocate(f(nf))
  allocate(wf(nf))
  allocate(e1d(nlam))
  allocate(e2d(nlam))
  allocate(e1mantle(nlam))
  allocate(e2mantle(nlam))

  ! ------------------------------------------------------------------------
  ! Definition of dust components volume fractions
  ! ------------------------------------------------------------------------

  ! Normalize the mass fractions
  tot = 0.0_dp
  do im=1,nm
     tot=tot+mfrac(im)
  enddo
  mfrac = mfrac/tot
  
  ! ------------------------------------------------------------------
  ! Identify the mantle material and take it out of the main list
  ! ------------------------------------------------------------------
  i_mantle = 0
  if (trim(loc(nm)).EQ."mantle") then
     ! Remember where the mantle information is
     mfrac_mantle = mfrac(nm)   ! So this is now the correct mass fraction
     mantle   = ref_index(nm)
     i_mantle = nm
     ! Reduce the number of materials, so that the Bruggeman
     ! rule is only applied to the other materials
     nm = nm-1
  endif

  ! ------------------------------------------------------------------
  ! Create the size distribution
  ! ------------------------------------------------------------------
  aminlog = log10(amin)
  amaxlog = log10(amax)
  pow  = -apow
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

  ! ------------------------------------------------------------------
  ! Get the refractory index data for all materials
  ! ------------------------------------------------------------------
  do im=1,nm
     call GetAndRegridLNK(ref_index(im),lam(1:nlam),e1d(1:nlam),e2d(1:nlam),nlam,.true.,rho(im))
     e1(im,1:nlam)    = e1d(1:nlam)
     e2(im,1:nlam)    = e2d(1:nlam)
  enddo
  if (i_mantle.gt.0) then
     call GetAndRegridLNK(ref_index(nm+1),lam(1:nlam),e1d(1:nlam),e2d(1:nlam),nlam,.true.,rho_mantle)
     e1mantle(1:nlam) = e1d(1:nlam)
     e2mantle(1:nlam) = e2d(1:nlam)
  endif
  
  deallocate(e1d)
  deallocate(e2d)

  ! Turn the mass fractions into volume fractions and compute rho_core
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
  
  min = dcmplx(1d0,0d0)
  ! ------------------------------------------------------------------
  ! Add vacuum as an extra material, to model porosity
  ! ------------------------------------------------------------------
  nm = nm+1
  e1(nm,1:nlam)     = 1.0_dp
  e2(nm,1:nlam)     = 0.0_dp
  rho(nm)           = 0.0_dp
  vfrac(nm)         = p_c
  vfrac(1:nm-1)     = vfrac(1:nm-1)*(1.0_dp-p_c) ! renormalize to 1
  rho_core          = rho_core * (1.d0-p_c)

  ! ------------------------------------------------------------------
  ! A this point, the core is properly normalized, including the
  ! porosity. Note that the mantle is no longer in the normalization
  ! mix at this point.
  ! We now compute the average density.
  ! ------------------------------------------------------------------

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
     ! FIXME: Still needs testing
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

  ! ------------------------------------------------------------------
  ! Now do the mixing
  ! ------------------------------------------------------------------

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
           ! The old Blender from OpacityTool.  Never fails, but not as good
           call Blender(vfrac,nm,epsj,eps_eff)
        endif
        e1(1,il)   = dreal(cdsqrt(eps_eff))
        e2(1,il)   = dimag(cdsqrt(eps_eff))
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
        call maxgarn_2compo(e1(1,il),e2(1,il),e1mantle(il),e2mantle(il),vfrac_mantle, &
             e1blend,e2blend)
        e1(1,il) = e1blend
        e2(1,il) = e2blend
     endif
  enddo ! end of wavelength loop over il
  
  ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! We are done with all the blending, so from now on we have ony one material
  !
  nm = 1
  ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
  ! Write the derived n and k to a file, if that is all we are supposed to do.
  if (blendonly) then
     write(*,'("Writing the blended n and k to blended.lnk, and exiting")')
     call remove_file_if_exists('blended.lnk')
     open(unit=20,file='blended.lnk')
     write(20,'(i5 f5.2)') nlam,rho_av
     do ilam=1,nlam
        write(20,'(3e15.4)') lam(ilam),e1(1,ilam),e2(1,ilam)
     enddo
     close(unit=20)
     stop
  endif

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
  
  if (nf.gt.1 .and. maxf.gt.0.01e0) then
     ! Get the weights for Gauss-Legendre integration
     call gauleg2(0.01e0,maxf,f(1:nf),wf(1:nf),nf)
  else if (maxf.eq.0e0) then
     ! Just a compact sphere, weight is 1
     f(1:nf)  = 0d0
     wf(1:nf) = 1d0/real(nf)
  else
     ! Just one fixed volume fraction of hollow
     f(1)  = maxf
     wf(1) = 1d0
  endif

  ! ------------------------------------------------------------------
  ! Start the main loop over all wavelengths
  ! ------------------------------------------------------------------
  do ilam = 1,nlam
     call tellertje(ilam,nlam)
     csca     = 0d0
     cabs     = 0d0
     cext     = 0d0
     mass     = 0d0
     vol      = 0d0
     
     do j=1,n_ang/2
        theta=(real(j)-0.5)/real(n_ang/2)*pi/2d0
        mu(j)=cos(theta)
     enddo
     
     ! ------------------------------------------------------------------
     ! Start the main loop over all particle sizes
     ! ------------------------------------------------------------------
     do is=1,ns
        r1       = r(is)
        err      = 0
        spheres  = 0
        toolarge = 0
        ! ------------------------------------------------------------------
        ! Start the loop over the DHS f factors
        ! ------------------------------------------------------------------
        cext_ff = 0.d0; cabs_ff = 0.d0; csca_ff = 0.d0
        do if=1,nf
           rad  = r1 / (1d0-f(if))**(1d0/3d0)
           m    = dcmplx(e1(1,ilam),-e2(1,ilam))
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
                   &                   n_ang/2, qext, qsca, qabs, gqsc, &
                   &                   m1, m2, s21, d21, n_ang ,err)
           else
              rcore = rad*0.999
              
              call DMiLay(rcore, rad, wvno, min, m, mu, &
                   &                   n_ang/2, qext, qsca, qabs, gqsc, &
                   &                   m1, m2, s21, d21, n_ang ,err)
           endif
           if (err.eq.1 .or. spheres.eq.1 .or. toolarge.eq.1) then
              rad   = r1
              rcore = rad
              rmie  = rad
              lmie  = lam(ilam)
              e1mie = e1(1,ilam)
              e2mie = e2(1,ilam)
              if (err.eq.1 .or. if.eq.1) then
                 if (rmie/lmie.lt.5000d0) then
                    call MeerhoffMie(rmie,lmie,e1mie,e2mie,csmie,cemie &
                         &      ,Mief11,Mief12,Mief33,Mief34,n_ang)
                 else
                    call MeerhoffMie(rmie,rmie/5000d0,e1mie,e2mie,csmie,cemie &
                         &      ,Mief11,Mief12,Mief33,Mief34,n_ang)
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
                 Mief11(n_ang-j+1) = (M2(j,2) + M1(j,2)) / csmie/wvno**2*2d0*pi
                 Mief12(n_ang-j+1) = (M2(j,2) - M1(j,2)) / csmie/wvno**2*2d0*pi
                 Mief22(n_ang-j+1) = (M2(j,2) + M1(j,2)) / csmie/wvno**2*2d0*pi
                 Mief33(n_ang-j+1) = (S21(j,2))          / csmie/wvno**2*2d0*pi
                 Mief34(n_ang-j+1) = (-D21(j,2))         / csmie/wvno**2*2d0*pi
                 Mief44(n_ang-j+1) = (S21(j,2))          / csmie/wvno**2*2d0*pi
              enddo
           endif
           
           ! make sure the scattering matrix is properly normalized by adjusting the forward peak.
           tot  = 0d0
           tot2 = 0d0
           do j=1,n_ang
              tot  = tot +  Mief11(j)*sin(pi*(real(j)-0.5)/real(n_ang))
              tot2 = tot2 + sin(pi*(real(j)-0.5)/real(n_ang))
           enddo
           Mief11(1) = Mief11(1) + (tot2-tot)/sin(pi*(0.5)/real(n_ang))
           if (Mief11(1).lt.0d0) Mief11(1) = 0d0
           
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
        enddo ! end loop nf over form factors
        
        cext = cext + cext_ff
        cabs = cabs + cabs_ff
        csca = csca + csca_ff

     enddo   ! end loop "is" over grain sizes
     
     ! ------------------------------------------------------------------
     ! Set the cross sections
     ! ------------------------------------------------------------------
     p%rho  = mass/vol
     p%Kext(ilam) = 1d4 * cext / mass
     p%Kabs(ilam) = 1d4 * cabs / mass
     p%Ksca(ilam) = 1d4 * csca / mass

     ! ------------------------------------------------------------------
     ! Set the elements of the scattering matrix
     ! ------------------------------------------------------------------
     p%F(ilam)%F11(1:180) = f11(ilam,1:180)/csca
     p%F(ilam)%F12(1:180) = f12(ilam,1:180)/csca
     p%F(ilam)%F22(1:180) = f22(ilam,1:180)/csca
     p%F(ilam)%F33(1:180) = f33(ilam,1:180)/csca
     p%F(ilam)%F34(1:180) = f34(ilam,1:180)/csca
     p%F(ilam)%F44(1:180) = f44(ilam,1:180)/csca
     tot = 0.0_dp
     do i=1,180
        p%g(ilam) = p%g(ilam) + p%F(ilam)%F11(i)*cos(pi*(real(i)-0.5)/180d0) &
             &                  *sin(pi*(real(i)-0.5)/180d0)
        tot = tot + p%F(ilam)%F11(i)*sin(pi*(real(i)-0.5)/180d0)
     enddo
     p%g(ilam) = p%g(ilam)/tot
  enddo
  deallocate(e1)
  deallocate(e2)
  deallocate(e1mantle)
  deallocate(e2mantle)
  
  deallocate(Mief11)
  deallocate(Mief12)
  deallocate(Mief22)
  deallocate(Mief33)
  deallocate(Mief34)
  deallocate(Mief44)
  deallocate(mu)
  deallocate(M1)
  deallocate(M2)
  deallocate(S21)
  deallocate(D21)
  
  deallocate(vfrac)
  deallocate(epsj)
  deallocate(f11)
  deallocate(f12)
  deallocate(f22)
  deallocate(f33)
  deallocate(f34)
  deallocate(f44)
  
  deallocate(r)
  deallocate(nr)
  deallocate(f)
  deallocate(wf)
  
  return
end subroutine ComputePart

! ------------------------------------------------------------------
! Routines needed to apply mixing rules
! ------------------------------------------------------------------

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
  !     nm      = number of grain material
  !     f       = volume fraction of each component
  !     e       = dielectric constant of each component
  !     epsavg  = averaged dielectric constant
  !**********************************************************************
  
  use Defs
  implicit none
  integer               :: nm, i, k, l, m
  integer               :: j(nm+1)
  real(kind =dp)        :: f(nm)
  COMPLEX (kind=dp)     :: e(nm)
  COMPLEX (kind=dp)     :: epsavg
  COMPLEX (kind=dp)     :: c(nm+1)
  COMPLEX (kind=dp)     :: prod
  COMPLEX (kind=dp)     :: x(nm)
  COMPLEX (kind=dp)     :: roots(nm)
  COMPLEX (kind=dp)     :: total !to check the result of Bruggeman rule
  logical polish
  polish=.false.
  
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
        write(*,*) "ERROR Bruggeman rule did not converge: no positive solution for effective refractive index: "
        write(*,*) roots(i)
        write(*,*) "Please try with a restricted range by using -lmin and -lmax"
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
  IMPLICIT NONE
  integer nm,j,iter
  real abun(nm)
  complex*16 e_in(nm),e_out
  complex*16 mm,m(nm),me,sum
  
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
  integer m,MAXM
  real(kind=dp) EPS
  COMPLEX(kind=dp) a(m+1),roots(m)
  logical polish
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
        !Polish the roots using the undeflated coefficients.
        call laguer(a,m,roots(j),its)
     enddo
  endif
  ! Adapted for SIGMA by Charlene Lefevre - only root with positive real part is kept
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

! ------------------------------------------------------------------
! Helper functions and subroutines
! ------------------------------------------------------------------

function cdlog10(x)
  ! Complex logarith base 10
  implicit none
  complex*16 x, cdlog10
  cdlog10=log(x)/log(10d0)
  return
end function cdlog10

subroutine tellertje(i,n)
  ! ------------------------------------------------------------------
  ! Show a progress bar on STDOUT
  ! ------------------------------------------------------------------
  use Defs
  implicit none
  integer :: i,n,f
  if(i.eq.1) write(*,'("....................")')
  f=int(20.0_dp*dble(i)/dble(n))
  if(20.0_dp*real(i-1)/real(n) .lt. real(f) &
       & .and. 20.0_dp*real(i+1)/real(n).GT.real(f)) then
     write(*,'(".",$)')
     call flush(6)
  endif
  if(i.eq.n) write(*,*)
  return
end subroutine tellertje

function count_words (string)
  !
  ! Count the whitespace-separated words in STRING. The separator can
  ! be an arbitrary number of space, tab, or newline characters
  !
  character*(*) :: string
  character*(30) :: s
  integer i, count_words
  logical in_space
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


subroutine remove_file_if_exists (file)
  ! Remove FILE if it already exists.  So we can cleanly overwrite it.
  character*(*) file
  logical file_exists
  inquire (file=trim(file),exist=file_exists)
  if (file_exists) then
     open(unit=90,file=trim(file))
     close(unit=90,status='delete')
  endif
end subroutine remove_file_if_exists

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

!
! Functions to help with the extended command line syntax, where
! we allow a single switch to have several values.
!
function arg_is_value(i)
  ! Check if command line argument is a value. That means the arg is
  ! there, and it is not a switch (starting with -)
  ! Problems: This gets it wrong for -.5
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

subroutine read_lambda_grid(file)
  ! Read the lambda grid from a file
  ! the file can start with comment lines (* or ! or #).
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


! ------------------------------------------------------------------
! Routines to write output files
! ------------------------------------------------------------------

subroutine write_header (unit,cc,amin,amax,apow,na,lmin,lmax,nlam, &
     p_core,p_mantle,fmax,nm,location,mfrac,mfrac_user,ref_index,rho)
  implicit none
  integer, parameter :: dp = selected_real_kind(P=15)
  integer :: unit
  integer :: na,nlam,nm,i
  real (kind=dp) :: amin,amax,apow,lmin,lmax,p_core,p_mantle,fmax
  real (kind=dp) :: mfrac(nm),mfrac_user(nm),rho(nm),amean(3)
  character*(*) location(nm),ref_index(nm),cc
  cc = trim(cc)
  call plmeans(amin,amax,apow,amean)
  write(unit,'(A,"============================================================================")') cc
  write(unit,'(A," Opacities computed by OpTool          <a^n>=",1p,3e9.2e1)') cc,amean
  write(unit,'(A," Parameters:")') cc
  write(unit,'(A,"   amin [um]=",f11.3," amax [um]=",f10.2,"  na  =",I5,"    apow=",g10.2)') cc,amin, amax, na, apow
  write(unit,'(A,"   lmin [um]=",f11.3," lmax [um]=",f10.2,"  nlam=",I5)') cc,lmin, lmax, nlam
  write(unit,'(A,"   porosity =",f11.3," p_mantle = ",f9.3,"            DHS fmax=",g10.2)') cc,p_core,p_mantle,fmax  
  write(unit,'(A," Composition:")') cc

  if (rho(1).gt.0) then
     ! We don't have rho yet, so we don't put it into the header.
     ! This is the situation where we write to the screen.
     write(unit,'(A,"  Where   mfrac  rho   Material")') cc
     write(unit,'(A,"  -----   -----  ----  ----------------------------------------------------")') cc
     do i=1,nm
        write(unit,'(A,"  ",A6,f7.3,f6.2,"  ",A)') cc,location(i), mfrac(i),rho(i),trim(ref_index(i))
     enddo
  else
     write(unit,'(A,"  Where   mfrac  Material")') cc
     write(unit,'(A,"  -----   -----  ----------------------------------------------------")') cc
     do i=1,nm
        write(unit,'(A,"  ",A6,f7.3,"  ",A)') cc,location(i), mfrac(i),trim(ref_index(i))
     enddo
  endif
  write(unit,'(A,"============================================================================")') cc
end subroutine write_header

subroutine write_ascii_file(p,amin,amax,apow,na,lmin,lmax,fmax,p_core,p_mantle,&
     nm,location,mfrac,mfrac_user,ref_index,rho,label,scatter,ext)
  use Defs
  implicit none
  real (kind=dp) :: amin,amax,apow,fmax,p_core,p_mantle
  real (kind=dp) :: mfrac(nm),mfrac_user(nm),rho(nm),lmin,lmax
  type(particle) p
  integer nm,na,i,j,nm2,ilam,iang
  character*(*) :: location(nm),ref_index(nm)
  character*(*) :: label,ext
  logical scatter
  character*500                :: file1,file2

  if (label .eq. ' ') then
     file1 = "dustkappa"      // ext
     file2 = "dustkapscatmat" // ext
  else
     file1 = "dustkappa_"   // trim(label) // ext
     file2 = "dustkapscatmat_" // trim(label) // ext
  endif
  
  call remove_file_if_exists(file1)
  call remove_file_if_exists(file2)

  if (.not. scatter) then
     write(*,'("Writing dust opacity output to file:  ",A)') trim(file1)
     open(20,file=file1,RECL=100000)
     call write_header(20,'#',amin,amax,apow,na,lmin,lmax,nlam, &
          p_core,p_mantle,fmax,nm,location,mfrac,mfrac_user,ref_index,rho)
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
     write(*,'("Writing full scattering data to file: ",A)') trim(file2)
     open(20,file=file2,RECL=100000)
     call write_header(20,'#',amin,amax,apow,na,lmin,lmax,nlam, &
          p_core,p_mantle,fmax,nm,location,mfrac,mfrac_user,ref_index,rho)
     write(20,'("# Output file formatted for RADMC-3D, dustkapscatmat, with scattering matrix")')
     write(20,'("#    iformat                                 ! iformat   must be 1")')
     write(20,'("#    nlam                                    ! number of wavelengths")')
     write(20,'("#    nang                                    ! number of angles 0-180")')
     write(20,'("#    lam(1)")')
     write(20,'("#    ...")')
     write(20,'("#    lam(nlam)")')
     write(20,'("#    ang(1)                                  ! ang(1)    must be 0")')
     write(20,'("#    ...")')
     write(20,'("#    ang(nang)                               ! ang(nang) must be 180")')
     write(20,'("#    F11 F12 F22 F33 F34 F44                 ! (ilam=   1,iang=   1)")')
     write(20,'("#    ...")')
     write(20,'("#    F11 F12 F22 F33 F34 F44                 ! (ilam=   1,iang=nang))")')
     write(20,'("#    F11 F12 F22 F33 F34 F44                 ! (ilam=   2,iang=   1))")')
     write(20,'("#    ...")')
     write(20,'("#    F11 F12 F22 F33 F34 F44                 ! (ilam=nlam,iang=nang))")')
     write(20,'("#============================================================================")')
     write(20,*) 1    ! iformat
     write(20,*) nlam ! Number of wavelength points
     write(20,*) 181  ! Number of angular points
     write(20,*)
     do ilam=1,nlam
        write(20,'(1p,e19.8,1p,e19.8,1p,e19.8,1p,e19.8)') lam(ilam),p%Kabs(ilam),p%Ksca(ilam),p%g(ilam)
     enddo
     write(20,*)
     do iang=0,180
        write(20,'(f8.2)') real(iang)
     enddo
     write(20,*)
     do ilam=1,nlam
        do iang=0,180
           ! We have only computed 0-179, but RADMC needs 180 as well
           ! We simply repeat the 179 value
           i = 1 + min(iang,179) 
           write(20,'(1p,e19.8,1p,e19.8,1p,e19.8,1p,e19.8,1p,e19.8,1p,e19.8)') &
                p%F(ilam)%F11(i),p%F(ilam)%F12(i),p%F(ilam)%F22(i), &
                p%F(ilam)%F33(i),p%F(ilam)%F34(i),p%F(ilam)%F44(i)
        enddo
     enddo
     close(20)
  endif
end subroutine write_ascii_file

#ifdef USE_FITSIO
subroutine write_fits_file(p,amin,amax,apow,na, &
     fmax,p_core,p_mantle, &
     nm,mfrac,ref_index,fitsfile)
  ! ------------------------------------------------------------------
  ! Routine to write a FITS file with all the information abou the particle.
  ! FIXME: Something goes wrong with the keywords
  ! ------------------------------------------------------------------
  use Defs
  implicit none
  character*500 ref_index(nm)
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
  print *,amean
  a1 = amean(1)
  call ftpkye(unit,'a1',real(a1),8,'[micron]',status)
!  call ftpkye(unit,'density',real(rho_av),8,'[g/cm^3]',status)
  
  call ftpkye(unit,'porosity',real(p_core),8,'[g/cm^3]',status)
  call ftpkye(unit,'p_mantle',real(p_mantle),8,'[g/cm^3]',status)
  
  do i=1,nm
     write(word,'("file",i0.2)') i
     call ftpkys(unit,word,trim(ref_index(i)),'',status)
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
  
  !------------------------------------------------------------------------------
  ! HDU 0: opacities 
  !------------------------------------------------------------------------------
  
  do i=1,nlam
     array(i,1,1)=lam(i)
     array(i,2,1)=p%Kext(i)
     array(i,3,1)=p%Kabs(i)
     array(i,4,1)=p%Ksca(i)
     array(i,5,1)=p%g(i)
  enddo
  
  call ftpprd(unit,group,fpixel,nelements,array(1:nlam,1:4,1),status)
  
  deallocate(array)
  
  !------------------------------------------------------------------------------
  ! HDU 1: scattering matrix
  !------------------------------------------------------------------------------
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
#else
subroutine write_fits_file(p,amin,amax,apow,na, &
     fmax,p_core,p_mantle, &
     nm,mfrac,ref_index,fitsfile)
  ! ------------------------------------------------------------------
  ! Just a dummy routine, in case there is not fits support
  ! We need this to keep the compiler happy
  ! ------------------------------------------------------------------
  use Defs
  implicit none
  character*500 ref_index(nm),fitsfile
  real (kind=dp) :: amin,amax,apow,fmax,p_core,p_mantle,mfrac(nm)
  integer :: nm,na
  type(particle) :: p
  call remove_file_if_exists(fitsfile)
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

  

