
!!! **** The Usage routine and the module for shared variables

subroutine usage()
  write(*,'("")')
  write(*,'("===============================================================================")')
  write(*,'("optool uses DHS (Min+2005) and effective medium theory to compute both the")')
  write(*,'("opacities and the scattering matrix for a size distribution of grains in a")')
  write(*,'("wavelength interval.  Without arguments, it produces the DIANA standard")')
  write(*,'("opacities (Woitke,Min+ 2016). With arguments, set your composition and other")')
  write(*,'("parameters.")')
  write(*,'("")')
  write(*,'("-c                        List available materials")')
  write(*,'("-c KEY-or-FILE [Mfrac]    Add material with mass fraction. -c may be omitted")')
  write(*,'("-m KEY-or-FILE [Mfrac]    Set material and mass fraction of mantle")')
  write(*,'("-p POROSITY [PMANTLE]     Set porosity, possibly different for core and mantle")')
  write(*,'("-fmax VHMAX               Maximum volume fraction of vacuum in DHS computation")')
  write(*,'("-a AMIN AMAX [APOW [NA]]  Set up grain radius radius grid (unit: micron)")')
  write(*,'("-l LMIN LMAX [NLAM]       Set up wavelength grid          (unit: micron)")')
  write(*,'("-l FILE                   Read wavelength grid from file, e.g. some/file.lnk")')
  write(*,'("-d [NSUB]                 Write NA files for specific grain sizes")')
  write(*,'("-s [NANG]                 Add scattering matrix for NANG angles to output")')
  write(*,'("-chop [NDEG]              Remove NDEG degrees from the forward-scattering peak")')
  write(*,'("-radmc [LABEL]; -fits     Special output options")')
  write(*,'("-o [DIRECTORY]            Output to DIRECTORY instead of current working dir")')
  write(*,'("-h; -man                  Show this message; Show the full manual")')
  write(*,'("===============================================================================")')
end subroutine usage

module Defs
  ! Module with constants and the basic data structures used in the program
  implicit none
  integer, public, parameter     :: dp = selected_real_kind(P=15)
  ! ----------------------------------------------------------------------
  ! Physics and math constants
  ! ----------------------------------------------------------------------
  real (kind=dp),public,parameter :: pi = 3.1415926535897932384_dp
  ! ----------------------------------------------------------------------
  ! Some global switches
  ! ----------------------------------------------------------------------
  logical, public                :: blendonly = .false. ! only blend materials and write result out
  logical, public                :: split     = .false. ! split to many files
  logical, public                :: quiet     = .false. ! reduce output to STDOUT
  logical, public                :: justnum   = .false. ! Only give opacities on STDOUT 
  logical, public                :: debug     = .false. ! Additional infor to STDOUT
  ! ----------------------------------------------------------------------
  ! Lambda is shared, because multiple routines need it
  ! ----------------------------------------------------------------------
  real (kind=dp), allocatable    :: lam(:)     ! wavelength array
  integer                        :: nlam       ! nr of wavelength points
  ! ----------------------------------------------------------------------
  ! Number of angles to be considered
  ! ----------------------------------------------------------------------
  integer                        :: nang
  real (kind=dp)                 :: chopangle
  ! ----------------------------------------------------------------------
  ! Material properties
  ! ----------------------------------------------------------------------
  integer         :: mat_nm       ! number of materials specified
  integer         :: mat_nmc      ! number of core materials specified
  integer         :: mat_nmm      ! number of mantle materials specified
  character*500   :: mat_loc(21)  ! either 'core' or 'mantle'
  character*500   :: mat_lnk(21)  ! the lnk key of file path
  real (kind=dp)  :: mat_rho(21)  ! specific mass density of material
  real (kind=dp)  :: mat_mfr(21)  ! mass fraction of each component
  real (kind=dp), allocatable  :: mat_e1(:,:)  ! Real      part of refractive index
  real (kind=dp), allocatable  :: mat_e2(:,:)  ! Imaginary part of refractive index
  ! ----------------------------------------------------------------------
  ! Mueller matrix structure, records only non-zero elements of the matrix
  ! ----------------------------------------------------------------------
  type mueller
     real (kind=dp),allocatable :: F11(:),F12(:),F22(:),F33(:),F44(:),F34(:)
  end type mueller
  ! ----------------------------------------------------------------------
  ! The particle structure, contains particle and scattering properties
  ! ----------------------------------------------------------------------
  type particle
     real (kind=dp)              :: rv,rvmin,rvmax ! grain radius with min and max
     real (kind=dp)              :: rho            ! mass density in g.cm^-3
     real (kind=dp), allocatable :: Kabs(:),Ksca(:),Kext(:) ! Opacities
     real (kind=dp), allocatable :: g(:)           ! asymmetry parameter
     TYPE(MUELLER),  allocatable :: F(:)           ! Mueller matrix elements
  end type particle
  ! ----------------------------------------------------------------------
  ! The output directory
  ! ----------------------------------------------------------------------
  character*500                  :: outdir = ''    ! Output directory
  character*3                    :: method         ! DHS or MMF

end module Defs

!!! **** Main program and ComputePart
program optool
  use Defs
  use omp_lib
  implicit none
  integer         :: na              ! nr of sizes for size distribution
  real (kind=dp)  :: amin,amax       ! min and max size of grains
  real (kind=dp)  :: apow            ! power law index f(a) ~ a^(-apow)
  real (kind=dp)  :: fmax            ! maximum fraction of vaccum for DHS
  real (kind=dp)  :: pcore, pmantle  ! porosity for core and mantle

  real (kind=dp)  :: lmin,lmax       ! min and max wavelength

  logical         :: write_mean_kap  ! Should mean kappas be computed?
  logical         :: write_scatter   ! Should a fits file be written?
  logical         :: write_fits      ! Should a fits file be written?
  logical         :: for_radmc       ! Should the scattering matrix use RADME-3D convention?

  integer         :: nm              ! nr of grain materials
  integer         :: nmant           ! nr of mantle materials
  integer         :: it,il

  type(particle)  :: p
  integer         :: i,ndone         ! counter
  integer         :: im,ia           ! for material, radius
  character*100   :: tmp,value       ! for processing args

  logical         :: arg_is_present  ! functions to test arguments
  logical         :: arg_is_switch   ! functions to test arguments
  logical         :: arg_is_value    ! functions to test arguments
  logical         :: arg_is_number   ! functions to test arguments
  character*500   :: fitsfile        ! file name for FITS output
  character*500   :: meanfile        ! file name for mean opacity output
  character*500   :: make_file_path  ! function
  character*50    :: radmclbl = ""   ! file label for RADMC-3D compatible
  character*50    :: label           ! for use in file names

  character*500   :: dumc            ! temporary storage
  real (kind=dp)  :: dum             ! temporary storage

  
  real (kind=dp)  :: asplit,afact,afsub,amaxsplit,aminsplit
  integer         :: nsubgrains = 5,nsub

  real (kind=dp), allocatable :: e1d(:),e2d(:)

  ! MMF implementation
  real(kind=dp)   :: mmf_a0
  

  ! ----------------------------------------------------------------------
  ! Defaults values for parameters and switches
  ! ----------------------------------------------------------------------
  amin           = 0.05_dp    ! micrometer
  amax           = 3000._dp   ! micrometer
  apow           = 3.50_dp    ! a minus sign will be added internally
  na             = 0          ! will be computed to 10 per decade

  lmin           = 0.05_dp    ! micrometer
  lmax           = 10000.0_dp ! micrometer
  nlam           = 300

  nang           = 180        ! Number of angular point, has to be even
  chopangle      = 0.d0       ! Angle in degree to chop forward peak
  
  fmax           = 0.8_dp     ! maximum volume fraction DHS
  pcore          = 0.0_dp     ! porosity core
  pmantle        = 0.0_dp     ! porosity mantle

  nm             = 0          ! number of materials - zero to start with
  nmant          = 0          ! number of mantle materials - zero to start with

  method         = 'DHS'      ! method to compute the opacities
  
  write_fits     = .false.    ! Default is to write ASCII output
  write_scatter  = .false.    ! Default is to not write scattering matrix
  write_mean_kap = .false.    ! Default is to not compute mean kappas
  for_radmc      = .false.    ! Default is to use optool conventions.

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

     case('?','-h','--help','-help','help')
        call usage(); stop

     case('??','-man','--man','-manual','--manual')
        call manual(); stop

        ! ----------------------------------------------------------------------
        ! Definition of a material for the mix
        ! ----------------------------------------------------------------------
     case('-c','-m')
        i  = i+1
        nm = nm+1;
        if (nm .gt. 20) then
           print *,'ERROR: too many materials'; stop
        endif

        ! First value is the material key or refindex file path
        call getarg(i,value); read(value,'(A)') mat_lnk(nm)

        if ((value .eq. '?') .or. (value .eq. '')) then
           call ListBuiltinMaterials(); stop
        endif
        ! Second value is the mass fraction
        if (.not. arg_is_number(i+1)) then
           if (.not. quiet) print *, "WARNING: 1.0 used for missing mass fraction of material: ",trim(mat_lnk(nm))
           mat_mfr(nm) = 1.0d0
        else
           i = i+1; call getarg(i,value); read(value,*) mat_mfr(nm)
        endif
        if (mat_mfr(nm).eq.0d0) then
           if (.not. quiet) print *, "WARNING: Ignoring material with zero mass fraction: ",trim(mat_lnk(nm))
           nm = nm-1
        else
           ! Set the type, and make sure we have at most one mantle material
           if (tmp.eq.'-m') then
              ! This is the mantle material
              mat_loc(nm) = 'mantle'
              nmant = nmant+1
           else
              ! This is a core material.
              mat_loc(nm) = 'core'
           endif
        endif

        ! There might be a density, as a third argument
        if (arg_is_number(i+1)) then
           i = i+1; call getarg(i,value); read(value,*) mat_rho(nm)
        endif

        ! ----------------------------------------------------------------------
        ! Grain size setup
        ! ----------------------------------------------------------------------
     case('-a')
        ! ----------------------------------------------------------------------
        ! -a expects 1, 2, 3, ro 4 values:  amin [amax [apow [na]]]
        ! ----------------------------------------------------------------------
        if (.not. arg_is_number(i+1)) then
           print *,"ERROR: -a needs 1-4 values: amin [amax [na [apow]]]"; stop
        endif
        i=i+1; call getarg(i,value); read(value,*) amin
        ! Let.s see if there is more, we expect amax
        if (arg_is_number(i+1)) then
           i=i+1; call getarg(i,value); read(value,*) amax
           ! Let.s see if there is more, we expect apow
           if (arg_is_number(i+1)) then
              i=i+1; call getarg(i,value); read(value,*) apow
              ! Let.s see if there is more, we expect na
              if (arg_is_number(i+1)) then
                 i=i+1; call getarg(i,value); read(value,*) na
              endif
           endif
        else
           ! there was only 1 number.  Set up single grain size computation
           ! If na has already been set, do not change it.
           amax = amin;
           if (na.eq.0) na = 1;
        endif
     case('-amin','--amin')
        i = i+1; call getarg(i,value); read(value,*) amin
     case('-amax','--amax')
        i = i+1; call getarg(i,value); read(value,*) amax
     case('-apow','--apow')
        i = i+1; call getarg(i,value); read(value,*) apow
     case('-na')
        i = i+1; call getarg(i,value); read(value,*) na

        ! ----------------------------------------------------------------------
        ! Wavelength setup
        ! ----------------------------------------------------------------------
     case('-l')
        ! ----------------------------------------------------------------------
        ! -l expects a file name, or 1-3 numbers: lmin [lmax [nlam]]
        ! ----------------------------------------------------------------------
        if (.not. arg_is_value(i+1)) then
           print *,"ERROR: -l needs a file or numbers as values"; stop
        else if (.not. arg_is_number(i+1)) then
           ! Could be a file name.  If yes, read the lambda grid from it
           call getarg(i+1,value)
           call check_for_file(trim(value))
           i=i+1
           call read_lambda_grid(trim(value))
        else
           ! We have a number, should be lmin
           i = i+1; call getarg(i,value); read(value,*) lmin
           ! Let's see if there is more, we expect lmax
           if (arg_is_number(i+1)) then
              i = i+1; call getarg(i,value); read(value,*) lmax
              ! Let.s see if there is more, we expect nlam
              if (arg_is_number(i+1)) then
                 i=i+1; call getarg(i,value); if (value.ne.'') read(value,*) nlam
              endif
           else
              ! only lmin was given, set up single lambda computation
              lmax = lmin; nlam = 1;
           endif
        endif
     case('-lmin','--lmin')
        i = i+1; call getarg(i,value); read(value,*) lmin
     case('-lmax','--lmax')
        i = i+1; call getarg(i,value); read(value,*) lmax
     case('-nlam','--nlam','-nl','--nl')
        i = i+1; call getarg(i,value); read(value,*) nlam

        ! ----------------------------------------------------------------------
        ! Grain geometry, including method DHS versus MMF
        ! ----------------------------------------------------------------------
     case('-p','-porosity','--porosity')
        i = i+1;  call getarg(i,value); read(value,*) pcore
        if (arg_is_number(i+1)) then
           i=i+1; call getarg(i,value); read(value,*) pmantle
        else
           pmantle = pcore
        endif
     case('-dhs','-fmax','-mie')
        method = 'DHS'
        if (tmp.eq.'-mie') then
           fmax = 0.
        else if (arg_is_number(i+1)) then
           i = i+1; call getarg(i,value); read(value,*) fmax
        else
           fmax = 0.8
        endif
     case('-mmf')
        method = 'MMF'
        if (arg_is_number(i+1)) then
           i=i+1; call getarg(i,value); read(value,*) mmf_a0
        else
           mmf_a0 = 1d0
        endif

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
        if (arg_is_number(i+1)) then
           i=i+1; call getarg(i,value); read(value,*) nang
        endif
     case('-chop')
        if (.not. arg_is_value(i+1)) then
           chopangle = 2.d0
        else
           i=i+1; call getarg(i,value); read(value,*) chopangle
        endif
     case('-radmc','-radmc3d')
        for_radmc = .true.
        if (arg_is_value(i+1)) then
           i=i+1
           call getarg(i,radmclbl)
        endif
     case('-fits')
        write_fits = .true.
     case('-d')
        split = .true.
        if (arg_is_value(i+1)) then
           i=i+1; call getarg(i,value); read(value,*) nsubgrains
        endif
     case('-b','-blendonly','--blendonly')
        ! Write blended refractive index to file and exit
        blendonly = .true.
     case('-q')
        ! Be less noisy
        quiet = .true.
     case ('-n')
        quiet   = .true.
        justnum = .true.
     case('-debug')
        ! More info to STDOUT
        debug = .true.
     case default
        if (arg_is_switch(i)) then
           write(*,*) "ERROR: Option or Arg: >",trim(tmp),'> not recognized'
           write(*,*) "For help, try: optool -h     ... or find the user guide OpTool.pdf"
           stop
        else
           if (debug) print *,trim(tmp),' will be interpreted as a material'
           tmp = 'FoRcE_-c_FoRcE'
        endif
     end select
     if (tmp .eq. 'FoRcE_-c_FoRcE') then
        tmp = '-c'
        i = i-1
     else
        i = i+1
        call getarg(i,tmp)   ! Get the next argument for the loop.
     endif
  enddo

  ! ----------------------------------------------------------------------
  ! Sanity checks and preparations
  ! ----------------------------------------------------------------------
  ! *** Materials ***
  if (nm .ge. 10) then
     print *,'ERROR: Too many materials'; stop
  endif
  if ( (nm.eq.nmant) .and. (nm.gt.0) ) then
     print *,"ERROR: at least one core material must be specified"; stop
  endif
  ! *** Grain size distribution ***
  if (na .eq. 0) then
     ! set sampling of the grain radius: 10 per decade, min 5
     na = max(5,int((log10(amax)-log10(amin))*10d0+1d0))
  endif
  if ( (amin.le.0d0) .or. (amax.le.0d0) ) then
     print *,'ERROR: Both amin and amax need to be positive numbers',amin,amax; stop
  endif
  if (amin .gt. amax) then
     ! Swap min and max values
     dum = amin; amin = amax; amax = dum
  endif
  if (apow.lt.0d0) then
     print *,'WARNING: Negative value for apow. apow=-3 means f(a)~a^(+3)'
  endif
  ! *** Wavelength grid ***
  if ( (lmin.le.0d0) .or. (lmax.le.0d0) ) then
     print *,'ERROR: Both lmin and lmax need to be positive numbers',lmin,lmax; stop
  endif
  if (lmin .gt. lmax) then
     ! Swap min and max values
     dum = lmin; lmin = lmax; lmax = dum
  endif
  if ( (nlam.le.1) .and. (lmin.ne.lmax)) then
     print *,'ERROR: More than one wavelength point needed to sample a range',lmin,lmax,nlam; stop
  endif
  if ( (lmin.eq.lmax) .and. (nlam.ne.1) ) then
     print *,'WARNING: Setting nlam=1 because lmin=lmax'
     nlam = 1
  endif
  ! *** Other checks
  if (split .and. write_mean_kap) then
     write(*,*) 'ERROR: With both -d and -t options, the dustkapmean.dat file'
     write(*,*) '       would only reflect the final size bin.'
     stop
  endif
  if (split .and. blendonly) then
     if (.not. quiet) write(*,*) 'WARNING: Turning off -s for -blendonly'
     split = .false.
  endif
  ! *** Angular grid ***
  if (mod(nang,2) .eq. 1) then
     write(*,*) 'ERROR: The number of angles in -s NANG must be even'
     stop
  endif
  ! *** Output files ***
#ifndef USE_FITSIO
  if (write_fits) then
     write(*,*) 'ERROR: Support for writing FITS files needs to be compiled in.'
     write(*,*) '       If you want FITS output, make sure cfitsio library is installed."'
     write(*,*) '       Then recompile with: "make clean", and then "make fitsio=true"'
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
     if (.not. quiet) then
        write(*,'("No materials specified, using DIANA standard")')
     endif
     nm = 2
     mat_lnk(1) = 'pyr-mg70' ; mat_loc(1)  = 'core' ; mat_mfr(1)     = 0.87d0
     mat_lnk(2) = 'c-z'      ; mat_loc(2)  = 'core' ; mat_mfr(2)     = 0.1301d0
     pcore      = 0.25d0
  endif
  mat_nm   = nm

  ! ----------------------------------------------------------------------
  ! Sort the materials to make sure that core materials come first
  ! ----------------------------------------------------------------------
  it  = nm  ! target: where to move the next mantle material
  mat_nmc = nm; mat_nmm = 0
  do il = nm,1,-1 ! il loops down, checking all materials
     if (mat_loc(il).eq.'mantle') then
        mat_nmc = mat_nmc-1
        mat_nmm = mat_nmm+1
        if (il.lt.it) then
           dumc = mat_lnk(il); mat_lnk(il)=mat_lnk(it); mat_lnk(it)=dumc;
           dumc = mat_loc(il); mat_loc(il)=mat_loc(it); mat_loc(it)=dumc;
           dum  = mat_mfr(il); mat_mfr(il)=mat_mfr(it); mat_mfr(it)=dum;
           dum  = mat_rho(il); mat_rho(il)=mat_rho(it); mat_rho(it)=dum;
           it = it-1
        endif
     endif
  enddo
  
  ! ----------------------------------------------------------------------
  ! Make a logarithmic lambda grid, unless read in from file
  ! ----------------------------------------------------------------------
  if (allocated(lam)) then
     ! Lam was allocated by reading from a file
     continue
  else
     allocate(lam(nlam))
     if (nlam.eq.1) then
        lam(1) = lmin
     else
        do i = 1,nlam
           lam(i)=10.0_dp**(log10(lmin)+log10(lmax/lmin)*(i-1)/(nlam-1))
        enddo
     endif
  endif
  
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
  ! Write a setup summary to the screen
  ! ----------------------------------------------------------------------
  if (.not. quiet) then
     call write_header(6,'',amin,amax,apow,na,lmin,lmax, &
          pcore,pmantle,0.0d0,fmax,mat_mfr,mat_nm)
  endif

  ! ----------------------------------------------------------------------
  ! Loop for splitting the output into files by grain size
  ! ----------------------------------------------------------------------
  if (split) then
     write(*,'("Computing opacities for ",I3," different grain size bins")') na
     nsub = nsubgrains
     if (mod(nsub,2).eq.0) nsub = nsub+1
     afact = (amax/amin)**(1.d0/real(na))   ! factor to next grain size
     afsub = afact**(1.d0/real(nsub-1))     ! Factor to next subgrain size
     ndone = 0; call tellertje(ndone,na,quiet)
     
     !$OMP parallel do if (split)                                      &
     !$OMP default(none)                                               &
     !$OMP shared(amin,afact,afsub,nsub,apow,fmax,pcore,pmantle)       &
     !$OMP shared(lmin,lmax,ndone,na,mat_mfr,mat_rho,mat_nm,nlam,nang) &
     !$OMP shared(outdir,write_scatter,for_radmc,write_fits,radmclbl)  &
     !$OMP shared(quiet)                                               &
     !$OMP private(ia,asplit,aminsplit,amaxsplit,label,fitsfile,p)
     do ia=1,na

        ! Allocate the particle structure
        ! We do this inside the loop, because iFort in combination with
        ! OpenMP crashes if we expect it to properly handle the composite
        ! and dynamically allocated p as a private variable.
        ! This is a bit wasteful in a non-parallel version of the code,
        ! because the allocation/deallocation is then repeated na times,
        ! unnecessarily.
        allocate(p%Kabs(nlam),p%Kext(nlam),p%Ksca(nlam),p%g(nlam),p%F(nlam))
        do i=1,nlam
           allocate(p%F(i)%F11(nang),p%F(i)%F12(nang),p%F(i)%F22(nang))
           allocate(p%F(i)%F33(nang),p%F(i)%F34(nang),p%F(i)%F44(nang))
        enddo

        asplit    = amin  *afact**real(ia-1d0+0.5d0)
        aminsplit = asplit*afsub**real(-nsub/2)
        amaxsplit = asplit*afsub**real(+nsub/2)
        call ComputePart(p,aminsplit,amaxsplit,apow,nsub,fmax,pcore,pmantle, &
             mat_mfr,mat_nm,.false.)

        ! Outout is done serially, to avoid file handle conflicts
        !$OMP critical
        ndone = ndone + 1
        call tellertje(ndone,na,quiet)
        write(label,'(I3.3)') ia
        if (write_fits) then
#ifdef USE_FITSIO
           write(fitsfile,'(A,"_",A,".fits")') "dustkappa",trim(label)
           fitsfile = make_file_path(outdir,fitsfile)
           call write_fits_file(p,aminsplit,amaxsplit,apow,nsub, &
                fmax,pcore,pmantle,mat_nm,mat_mfr,mat_rho,fitsfile)
#endif
        else
           if (radmclbl .ne. ' ') label = trim(radmclbl) // "_" // label
           call write_ascii_file(p,aminsplit,amaxsplit,apow,nsub,lmin,lmax, &
                fmax,pcore,pmantle,mat_mfr,mat_nm, &
                label,write_scatter,for_radmc,.false.)
        endif
        !$OMP end critical
        do i=1,nlam
           deallocate(p%F(i)%F11,p%F(i)%F12,p%F(i)%F22,p%F(i)%F33,p%F(i)%F34,p%F(i)%F44)
        enddo
        deallocate(p%Kabs,p%Kext,p%Ksca,p%g,p%F)
     enddo
     !$OMP end parallel DO
     stop

  else

     ! Allocate the particle structure
     allocate(p%Kabs(nlam),p%Kext(nlam),p%Ksca(nlam),p%g(nlam),p%F(nlam))
     do i=1,nlam
        allocate(p%F(i)%F11(nang),p%F(i)%F12(nang),p%F(i)%F22(nang))
        allocate(p%F(i)%F33(nang),p%F(i)%F34(nang),p%F(i)%F44(nang))
     enddo
     ! ----------------------------------------------------------------------
     ! Call the main routine to compute the opacities and scattering matrix
     ! ----------------------------------------------------------------------
     call ComputePart(p,amin,amax,apow,na,fmax,pcore,pmantle,mat_mfr,nm,.true.)
     
     ! ----------------------------------------------------------------------
     ! Write the output
     ! ----------------------------------------------------------------------
     if (write_fits) then
#ifdef USE_FITSIO
        call write_fits_file(p,amin,amax,apow,na, &
             fmax,pcore,pmantle,mat_nm,mat_mfr,mat_rho,fitsfile)
#endif
        continue
     else
        if (justnum) then
           call write_to_sdtout(p)
        else
           call write_ascii_file(p,amin,amax,apow,na,lmin,lmax, &
                fmax,pcore,pmantle,mat_mfr,mat_nm, &
                radmclbl,write_scatter,for_radmc,.true.)
        endif
     endif

     ! Deallocate the particle structure
     do i=1,nlam
        deallocate(p%F(i)%F11,p%F(i)%F12,p%F(i)%F22,p%F(i)%F33,p%F(i)%F34,p%F(i)%F44)
     enddo
     deallocate(p%Kabs,p%Kext,p%Ksca,p%g,p%F)
     
  endif
  
end program optool

! **** ComputePart, the central routine avaraging properties over sizes

subroutine ComputePart(p,amin,amax,apow,na,fmax,p_c,p_m,mfrac0,nm,progress)
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
  !   nm        The number of materials, also the size of the array mfrac
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

  integer                        :: nm,im            ! nr of grain materials in a composite grain
  real (kind=dp),allocatable     :: rho(:)           ! specific material dnesities
  real (kind=dp)                 :: mfrac0(nm)       ! mass fractions, input
  real (kind=dp),allocatable     :: mfrac(:)         ! mass fractions, modified
  real (kind=dp)                 :: vfrac_mantle     ! volume fraction of mantle
  real (kind=dp)                 :: mfrac_mantle     ! mass   fraction of mantle

  TYPE (PARTICLE)                :: p

  logical                        :: progress

  integer                        :: i,j              ! counters for various loops
  integer                        :: na               ! Number of grains sizes between amin and amax
  integer                        :: nf,if            ! Number of DHS volume fractions
  integer                        :: ns,is            ! Number of grains sizes
  integer                        :: ilam,il,ndone    ! Counter for wavelengths
  integer                        :: err,spheres,toolarge ! Error control for Mie routines

  real (kind=dp)                 :: cext, csca, cabs
  real (kind=dp)                 :: qext, qsca, qabs, gqsc

  real (kind=dp), allocatable    :: f11(:),    f12(:),    f22(:),    f33(:),    f34(:),    f44(:)
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
  real (kind=dp), allocatable    :: vfrac(:),vfm(:)
  complex (kind=dp), allocatable :: e_in(:)
  complex (kind=dp)              :: m,min,e_out


  integer ichop
  
  ns     = na    ! number of subgrains to compute

  ! ----------------------------------------------------------------------
  ! Allocate the necessary arrays
  ! ----------------------------------------------------------------------
  allocate(Mief11(nang),Mief12(nang),Mief22(nang),Mief33(nang))
  allocate(Mief34(nang),Mief44(nang))
  allocate(mu(nang),M1(nang,2),M2(nang,2),S21(nang,2),D21(nang,2))

  allocate(vfrac(nm),vfm(nm),mfrac(nm),rho(nm))
  allocate(e_in(nm))
  allocate(f11(nang),f12(nang),f22(nang))
  allocate(f33(nang),f34(nang),f44(nang))

  allocate(e1(nm,nlam),e2(nm,nlam))
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
  mfrac(1:nm) = mfrac0(1:nm)
  rho(1:nm)   = mat_rho(1:nm)
  
  ! ----------------------------------------------------------------------
  ! Normalize the mass fractions
  ! ----------------------------------------------------------------------
  tot  = 0.0_dp
  tot2 = 0.0_dp
  do im=1,nm
     tot=tot+mfrac(im)
     if (mat_loc(im).eq.'mantle') tot2 = tot2 + mfrac(im)
  enddo
  mfrac = mfrac/tot
  mfrac_mantle = tot2/tot

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
        nr(is)=1.d0*nr(is)/tot
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

  ! ----------------------------------------------------------------------
  ! Core: Turn mass fractions into volume fractions, compute rho_core
  ! ----------------------------------------------------------------------
  mtot = 0.d0
  vtot = 0.d0
  do im = 1,mat_nmc
     mtot = mtot + mfrac(im)
     vfrac(im) = mfrac(im)/rho(im)
     vtot      = vtot+vfrac(im)
  enddo
  rho_core  = mtot/vtot  ! No porosity included yet, will be done below
  ! Normalize the volume fractions of the core
  vfrac(1:mat_nmc) = vfrac(1:mat_nmc)/vtot
  if (p_c .gt. 0.d0) then
     vfrac(1:mat_nmc) = vfrac(1:mat_nmc)*(1.0_dp-p_c)
     rho_core         = rho_core        *(1.0_dp-p_c)
  endif

  ! ----------------------------------------------------------------------
  ! Mantle: Turn mass fractions into volume fractions, compute rho_mantle
  ! ----------------------------------------------------------------------
  if (mat_nmm .gt. 0) then
     ! we do have mantle stuff
     mtot = 0.d0
     vtot = 0.d0
     do im = mat_nmc+1,nm
        mtot = mtot + mfrac(im)
        vfrac(im) = mfrac(im)/rho(im)
        vtot      = vtot+vfrac(im)
     enddo
     rho_mantle = mtot/vtot  ! No porosity included yet, will be done below
     ! Normalize the volume fractions of the core
     vfrac(mat_nmc+1:nm) = vfrac(mat_nmc+1:nm)/vtot
     if (p_m .gt. 0.d0) then
        vfrac(mat_nmc+1:nm) = vfrac(mat_nmc+1:nm)*(1.0_dp-p_m)
        rho_mantle = rho_mantle *(1.0_dp-p_m)
     endif
  endif

  ! ----------------------------------------------------------------------
  ! Compute the average densities of the whole grain
  ! ----------------------------------------------------------------------

  if (mat_nmm.eq.0) then
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
     rho_av        = rho_core / (1.d0 + mfrac_mantle*(rho_core/rho_mantle-1.d0))
     vfrac_mantle  = mfrac_mantle * rho_av / rho_mantle
  endif

  ! ----------------------------------------------------------------------
  ! Now do the mixing, for all wavelengths
  ! ----------------------------------------------------------------------

  ! Loop over wavelengths
  do il=1,nlam
     ! ----------------
     ! Core
     ! ----------------
     if (nm.eq.1 .and. p_c.eq.0) then
        ! Solid core, single material, nothing to blend for the core
        e1blend(il) = e1(1,il)
        e2blend(il) = e2(1,il)
     else
        ! Blend the core materials
        do im=1,mat_nmc
           e_in(im) = dcmplx(e1(im,il),e2(im,il))
        enddo
        call Blender_vac(vfrac,mat_nmc,e_in,e_out)
        e1blend(il) = dreal(e_out)
        e2blend(il) = dimag(e_out)
     endif
     ! ----------------
     ! Mantle
     ! ----------------
     if (mat_nmm.gt.0) then
        ! We do have a mantle to add
        if ( (mat_nmm.eq.1) .and. (p_m.eq.0) ) then
           ! No Blending needed inside the mantle - just copy e1 and e2
           ! Since it is onyl one material, we know it is index nm
           e1mantle(il) = e1(nm,il)
           e2mantle(il) = e2(nm,il)
        else
           ! Blend the mantle materials
           do im=1,mat_nmm
              e_in(im) = dcmplx(e1(im+mat_nmc,il),e2(im+mat_nmc,il))
              vfm(im)  = vfrac(im+mat_nmc)
           enddo
           call Blender_vac(vfm,mat_nmm,e_in,e_out)
           e1mantle(il) = dreal(e_out)
           e2mantle(il) = dimag(e_out)
        endif

        ! Now we have the mantle material ready - put it on the core
        call Blender_MG(e1blend(il),e2blend(il),e1mantle(il),e2mantle(il), &
             vfrac_mantle,e1mg,e2mg)
        e1blend(il) = e1mg
        e2blend(il) = e2mg
     endif
  enddo ! end of wavelength loop over il

  ! ----------------------------------------------------------------------
  ! Write the derived n and k to a file
  ! ----------------------------------------------------------------------
  if (blendonly) then
     write(*,'("Writing the blended n and k to blended.lnk, and exiting")')
     call remove_file_if_exists('blended.lnk')
     open(unit=20,file='blended.lnk')
     write(20,'(i5,f5.2)') nlam,rho_av
     do ilam=1,nlam
        write(20,'(1p,e15.4,1p,e15.4,1p,e15.4)') lam(ilam),e1blend(ilam),e2blend(ilam)
     enddo
     close(unit=20)
     stop
  endif

  ! ----------------------------------------------------------------------
  ! Check how we are going to average over hollow sphere components
  ! ----------------------------------------------------------------------
  if (nf.gt.1 .and. fmax.gt.0.01e0) then
     ! Get the weights for Gauss-Legendre integration
     call gauleg2(0.01d0,fmax,f(1:nf),wf(1:nf),nf)
  else if (fmax.eq.0e0) then
     ! Just a compact sphere, weight is 1
     f(1:nf)  = 0d0
     wf(1:nf) = 1d0/real(nf)
  else
     ! Just one fixed volume fraction of hollow
     f(1)  = fmax
     wf(1) = 1d0
  endif

  ! ----------------------------------------------------------------------
  ! Initialize mu
  ! ----------------------------------------------------------------------
  do j=1,nang/2
     theta = (real(j)-0.5d0)/real(nang/2)*pi/2d0
     mu(j) = cos(theta)
  enddo
  
  ! ----------------------------------------------------------------------
  ! Start the main loop over all wavelengths
  ! ----------------------------------------------------------------------
  ndone = 0; if (progress) call tellertje(ndone,nlam,quiet)
  !$OMP parallel do if (.not. split)                                      &
  !$OMP default(none)                                                     &
  !$OMP private(f11,f12,f22,f33,f34,f44)                                  &
  !$OMP shared(r,lam,nlam,mu,e1blend,e2blend,p,nr,method,nf,ns,rho_av,wf,f) &
  !$OMP shared(split,progress,ndone,nang,chopangle)                       &
  !$OMP shared(quiet)                                                     &
  !$OMP private(r1,is,if,rcore,rad,ichop)                                 &
  !$OMP private(csca,cabs,cext,mass,vol)                                  &
  !$OMP private(cemie,csmie,e1mie,e2mie,rmie,lmie,qabs,qsca,qext,gqsc)    &
  !$OMP private(err,spheres,toolarge)                                     &
  !$OMP private(m1,m2,d21,s21,m,wvno,min)                                 &
  !$OMP private(Mief11,Mief12,Mief22,Mief33,Mief34,Mief44)                &
  !$OMP private(tot,tot2)
  do ilam = 1,nlam

     ! ----------------------------------------------------------------------
     ! Initialize the scattering matrix elements and summing variables
     ! ----------------------------------------------------------------------
     do j=1,nang
        f11(j) = 0d0; f12(j) = 0d0; f22(j) = 0d0
        f33(j) = 0d0; f34(j) = 0d0; f44(j) = 0d0
     enddo
     csca = 0d0; cabs = 0d0; cext = 0d0
     mass = 0d0; vol  = 0d0

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
        if (method .eq. 'DHS') then
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
              else if (method(1:3) .eq. 'DHS') then 
                 rcore = rad*f(if)**(1d0/3d0)
                 call DMiLay(rcore, rad, wvno, m, min, mu, &
                      nang/2, qext, qsca, qabs, gqsc, &
                      m1, m2, s21, d21, nang ,err)
              else
                 rcore = rad*0.999d0
                 call DMiLay(rcore, rad, wvno, min, m, mu, &
                      nang/2, qext, qsca, qabs, gqsc, &
                      m1, m2, s21, d21, nang ,err)
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
                            Mief11,Mief12,Mief33,Mief34,nang)
                    else
                       call MeerhoffMie(rmie,rmie/5000d0,e1mie,e2mie,csmie,cemie, &
                            Mief11,Mief12,Mief33,Mief34,nang)
                    endif
                 endif
                 
                 Mief22 = Mief11
                 Mief44 = Mief33
                 
              else
                 cemie = qext * pi * rad**2
                 csmie = qsca * pi * rad**2
                 do j=1,nang/2
                    Mief11(j) = (M2(j,1) + M1(j,1))        / csmie/wvno**2*2d0*pi
                    Mief12(j) = (M2(j,1) - M1(j,1))        / csmie/wvno**2*2d0*pi
                    Mief22(j) = (M2(j,1) + M1(j,1))        / csmie/wvno**2*2d0*pi
                    Mief33(j) = (S21(j,1))                 / csmie/wvno**2*2d0*pi
                    Mief34(j) = (-D21(j,1))                / csmie/wvno**2*2d0*pi
                    Mief44(j) = (S21(j,1))                 / csmie/wvno**2*2d0*pi
                    ! Here we use the assumption that the grid is regular.  An adapted
                    ! grid is not possible if it is not symmetric around pi/2.
                    Mief11(nang-j+1) = (M2(j,2) + M1(j,2)) / csmie/wvno**2*2d0*pi
                    Mief12(nang-j+1) = (M2(j,2) - M1(j,2)) / csmie/wvno**2*2d0*pi
                    Mief22(nang-j+1) = (M2(j,2) + M1(j,2)) / csmie/wvno**2*2d0*pi
                    Mief33(nang-j+1) = (S21(j,2))          / csmie/wvno**2*2d0*pi
                    Mief34(nang-j+1) = (-D21(j,2))         / csmie/wvno**2*2d0*pi
                    Mief44(nang-j+1) = (S21(j,2))          / csmie/wvno**2*2d0*pi
                 enddo
              endif    ! (err.eq.1 .or. spheres.eq.1 .or. toolarge.eq.1)
              
              ! make sure the scattering matrix is properly normalized by adjusting the forward peak.
              ! ASKMICHIEL: You are only adjusting the one point a 0 degrees.
              ! So that means that the basic normalization is already in place at this point,
              ! and you just want to be sure that things are right on your specific grid?
              tot  = 0d0
              tot2 = 0d0
              do j=1,nang
                 ! This integration assumes that the grid is regular (linear)
                 tot  = tot +  Mief11(j)*sin(pi*(real(j)-0.5d0)/real(nang))
                 tot2 = tot2 + sin(pi*(real(j)-0.5d0)/real(nang))
              enddo
              Mief11(1) = Mief11(1) + (tot2-tot)/sin(pi*(0.5d0)/real(nang))
              if (Mief11(1) .lt. 0d0) Mief11(1) = 0d0
              
              ! Add this contribution with the proper weights
              do j=1,nang
                 f11(j) = f11(j) + wf(if)*nr(is)*Mief11(j)*csmie
                 f12(j) = f12(j) + wf(if)*nr(is)*Mief12(j)*csmie
                 f22(j) = f22(j) + wf(if)*nr(is)*Mief22(j)*csmie
                 f33(j) = f33(j) + wf(if)*nr(is)*Mief33(j)*csmie
                 f34(j) = f34(j) + wf(if)*nr(is)*Mief34(j)*csmie
                 f44(j) = f44(j) + wf(if)*nr(is)*Mief44(j)*csmie
              enddo
              cext = cext + wf(if)*nr(is)*cemie
              csca = csca + wf(if)*nr(is)*csmie
              cabs = cabs + wf(if)*nr(is)*(cemie-csmie)
              mass = mass + wf(if)*nr(is)*rho_av*4d0*pi*r1**3/3d0
              vol  = vol  + wf(if)*nr(is)*4d0*pi*r1**3/3d0
           enddo    ! end loop "nf" over form factors
        else if (method .eq. 'MMF') then
           !! FIXME: Here we do the MMF magic
           !! We need to fill cext, csca, cabs, mass,vol,f11 ... f44
           !! Make use of nr(is), to be sure
           !! Use equations to compute d_f and k0 from a and p
        else
           print *,"ERROR: invalid method ", method
        endif

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
     p%F(ilam)%F11(1:nang) = f11(1:nang)/csca
     p%F(ilam)%F12(1:nang) = f12(1:nang)/csca
     p%F(ilam)%F22(1:nang) = f22(1:nang)/csca
     p%F(ilam)%F33(1:nang) = f33(1:nang)/csca
     p%F(ilam)%F34(1:nang) = f34(1:nang)/csca
     p%F(ilam)%F44(1:nang) = f44(1:nang)/csca

     ! ----------------------------------------------------------------------
     ! Chop off forward-scattering peak, adjust the scattering cross section
     ! ----------------------------------------------------------------------
     if (chopangle .gt. 0d0) then
        ! Flatten the peak
        ichop =  int( chopangle / (180.d0/nang) )
        if (ichop .gt. 1) then
           do i=1,ichop
              p%F(ilam)%F11(i) = p%F(ilam)%F11(ichop+1)
              p%F(ilam)%F12(i) = p%F(ilam)%F12(ichop+1)
              p%F(ilam)%F22(i) = p%F(ilam)%F22(ichop+1)
              p%F(ilam)%F33(i) = p%F(ilam)%F33(ichop+1)
              p%F(ilam)%F34(i) = p%F(ilam)%F34(ichop+1)
              p%F(ilam)%F44(i) = p%F(ilam)%F44(ichop+1)
           enddo
        endif
        ! Integrate F11
        tot = 0; tot2 = 0
        do j=1,nang
           !                    F11                         sin theta                  d theta
           tot  = tot  +  p%F(ilam)%F11(j) * (sin(pi*(real(j)-0.5d0)/real(nang))) * (pi/dble(nang)) * (2.d0*pi)
           tot2 = tot2 +                      sin(pi*(real(j)-0.5d0)/real(nang))  * (pi/dble(nang)) * (2.d0*pi)
        enddo
        ! Scale the scattering matrix
        do j=1,nang
           p%F(ilam)%F11(j) = p%F(ilam)%F11(j) * tot2/tot
           p%F(ilam)%F12(j) = p%F(ilam)%F12(j) * tot2/tot
           p%F(ilam)%F22(j) = p%F(ilam)%F22(j) * tot2/tot
           p%F(ilam)%F33(j) = p%F(ilam)%F33(j) * tot2/tot
           p%F(ilam)%F34(j) = p%F(ilam)%F34(j) * tot2/tot
           p%F(ilam)%F44(j) = p%F(ilam)%F44(j) * tot2/tot
        enddo
        ! Scale kappa_scat, and change kappa_ext accordingly
        p%Ksca(ilam) = p%Ksca(ilam) * tot/tot2
        p%Kext(ilam) = p%Ksca(ilam) + p%Kabs(ilam)
        
        !tot = 0; tot2 = 0
        !do j=1,nang
        !   tot  = tot  +  p%F(ilam)%F11(j) * (sin(pi*(real(j)-0.5d0)/real(nang))) * (pi/dble(nang)) * (2.d0*pi)
        !   tot2 = tot2 +                      sin(pi*(real(j)-0.5d0)/real(nang))  * (pi/dble(nang)) * (2.d0*pi)
        !enddo
        !print *,'>',tot,tot2
     endif   ! chopangle .gt. 0

     ! ----------------------------------------------------------------------
     ! Average over angles to compute asymmetry factor g
     ! ----------------------------------------------------------------------
     ! A regular angular grid is assumed for this computation
     tot = 0.0_dp
     p%g(ilam) = 0.d0
     do i=1,nang
        p%g(ilam) = p%g(ilam) + p%F(ilam)%F11(i)*cos(pi*(real(i)-0.5d0)/dble(nang)) &
             *sin(pi*(real(i)-0.5d0)/dble(nang))
        tot = tot + p%F(ilam)%F11(i)*sin(pi*(real(i)-0.5d0)/dble(nang))
     enddo
     p%g(ilam) = p%g(ilam)/tot
     
     !$OMP atomic
     ndone = ndone+1
     if (progress) call tellertje(ndone,nlam,quiet)

  enddo   ! end loop ilam over wavelength
  !$OMP end parallel DO

  deallocate(e1,e2)
  deallocate(e1mantle,e2mantle)

  deallocate(Mief11,Mief12,Mief22,Mief33,Mief34,Mief44)
  deallocate(mu,M1,M2,S21,D21)

  deallocate(vfrac,vfm)
  deallocate(e_in)
  deallocate(f11,f12,f22,f33,f34,f44)

  deallocate(r,nr)
  deallocate(f,wf)

  return
end subroutine ComputePart

!!! **** Effective medium routines

subroutine blender(abun,nm,e_in,e_out)
  ! ----------------------------------------------------------------------
  ! This is the original blender routine used in OpacityTool
  ! ----------------------------------------------------------------------
  implicit none
  integer, parameter :: dp = selected_real_kind(P=15)
  integer            :: nm,j,iter
  real (kind=dp)     :: abun(nm)
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
  if ( abs(sum)/abs(mm).gt.1d-6 ) then
     print *,'WARNING: Blender might not be converged (mm,sum)',mm,sum
  endif
  e_out = me
end subroutine blender

subroutine Blender_vac(abun,nm,e_in,e_out)
  ! ----------------------------------------------------------------------
  ! Blend the materials using the Bruggeman rule.
  ! If the abundances do not add up to 1.0, fill the rest with vacuum.
  ! ----------------------------------------------------------------------
  implicit none
  integer, parameter :: dp = selected_real_kind(P=15)
  integer            :: nm,j,iter
  real (kind=dp)     :: abun(nm),abunvac
  complex (kind=dp)  :: e_in(nm),e_out,mvac
  complex (kind=dp)  :: mm,m(nm),me,tot

  ! Compute the volume of vacuum, with sanity check
  abunvac = 1.d0 - sum(abun(1:nm))
  mvac = dcmplx(1d0,0d0)
  if (abunvac .lt. 0.d0) then
     if (abs(abunvac).lt.1d-5) then
        ! Just a rounding error, fix it
        abunvac = 0.d0
     else
        ! Abundances have not been normalized properly, this is bad
        print *,"ERROR: Abundances not normalized in routine blender_vac"
        print *,abun
        print *,abunvac
        stop
     endif
  endif

  ! Do the blending iteratively
  mm   = mvac
  m    = e_in
  do iter=1,100
     tot = ((mvac**2-mm**2)/(mvac**2+2d0*mm**2))*abunvac
     do j=1,nm
        tot = tot + ((m(j)**2-mm**2)/(m(j)**2+2d0*mm**2))*abun(j)
     enddo
     me = (2d0*tot+1d0)/(1d0-tot)
     me = mm*cdsqrt(me)
     mm = me       
  enddo
  if ( abs(tot)/abs(mm).gt.1d-6 ) then
     print *,'WARNING: Blender might not be converged (mm,tot)',mm,tot
  endif
  e_out = me
end subroutine Blender_vac

subroutine Blender_MG(e1in,e2in,e1in_m,e2in_m,vf_m,e1out,e2out)
  ! 2 component Maxwell-Garnet mixing
  ! vf_m is the volume fration of the mantle material
  use Defs
  implicit none
  real (kind=dp) e1in,e2in,e1out,e2out,e1in_m,e2in_m,vf_m,vf_c
  complex (kind=dp) m1,m2,me,sqme

  m2 = dcmplx(e1in_m,e2in_m) ! mantle = coating material = matrix
  m1 = dcmplx(e1in,e2in)     ! inner core is the "inclusion"
  vf_c = 1.d0-vf_m           ! the volume fraction of the core is 1-vf_m
  
  me = m2**2 * (  (2d0*m2**2 + m1**2 - 2d0*vf_c * (m2**2-m1**2) ) &
       &        / (2d0*m2**2 + m1**2 +     vf_c * (m2**2-m1**2) ) )

  sqme  = cdsqrt(me)
  e1out = dreal(sqme); e2out = dimag(sqme)
end subroutine Blender_MG

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

!!! **** Helper functions and subroutines

function cdlog10(x)
  ! Complex logarith base 10
  implicit none
  complex*16 x, cdlog10
  cdlog10=log(x)/log(10d0)
  return
end function cdlog10

subroutine tellertje(i,n,quiet)
  ! ----------------------------------------------------------------------
  ! Show a progress bar on STDOUT
  ! ----------------------------------------------------------------------
  implicit none
  integer :: i,n,f,l,ndots,maxdots=20,mindots=20
  logical :: quiet
  ndots = max(mindots,min(n,maxdots))
  if (quiet) then
     ! do nothong
  else if(i.eq.0) then
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

!!! **** Extended command line argument syntax support

function arg_is_present(i)
  ! Check if command line argument is present.
  implicit none
  integer i
  logical arg_is_present
  character*2 :: value
  call getarg(i,value)
  arg_is_present = (.not. (value .eq. ''))
end function arg_is_present

function arg_is_switch(i)
  ! Check if command line argument is a switch. That means that
  ! the arg is there, starts with dash, followed by a letter
  implicit none
  integer i
  logical arg_is_switch
  character*2 :: value
  call getarg(i,value)
  if ( (len_trim(value).ge.2) &
       .and. (value(1:1).eq.'-') &
       .and. ( ((value(2:2).ge.'a') .and. (value(2:2).le.'z')) &
       .or.    ((value(2:2).ge.'A') .and. (value(2:2).le.'Z')) ) ) then
     arg_is_switch = .true.
  else
     arg_is_switch = .false.
  endif
end function arg_is_switch

function arg_is_value(i)
  ! Check if command line argument is a value as oposed to being
  ! a switch. That means the arg is there, and it is not a switch
  implicit none
  integer i
  logical arg_is_value,arg_is_switch
  character*2 :: value
  call getarg(i,value)
  if ((value.eq.'') .or. arg_is_switch(i)) then
     arg_is_value = .false.
  else
     arg_is_value = .true.
  endif
end function arg_is_value

function arg_is_number(i)
  ! Check if command line arg i is a number, i.e. it is there and
  ! starts with [0-9] or .[0-9] or -[0-9] or -.[0-9]
  implicit none
  integer i,ic
  logical arg_is_number
  character*3 :: value
  call getarg(i,value)
  arg_is_number = .false.
  ic = 1
  if (len_trim(value).gt.0) then
     if (len_trim(value).eq.2) then
        if ((value(ic:ic).eq.'-') .or. (value(ic:ic).eq.'.')) ic = 2
     else ! 3 chars at least
        if (value(ic:ic).eq.'-') ic = ic+1
        if (value(ic:ic).eq.'.') ic = ic+1
     endif
     arg_is_number = ((value(ic:ic).ge.'0') .and. (value(ic:ic).le.'9'))
  endif
end function arg_is_number

!!! **** Reading files and checking for some consistency

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


!!! **** Routines to write output files

subroutine write_to_sdtout(p)
  ! ----------------------------------------------------------------------
  ! Just write the most important numbers to STDOUT
  ! Each line has lambda kabs ksca kext g
  ! ----------------------------------------------------------------------
  use Defs
  implicit none
  type(particle) :: p
  integer i
  do i=1,nlam
     write(6,'(1p,5e12.3)') lam(i),p%kabs(i),p%ksca(i),p%kext(i),p%g(i)
  enddo
end subroutine write_to_sdtout

subroutine write_header (unit,cc,amin,amax,apow,na,lmin,lmax, &
     pcore,pmantle,rho_av,fmax,mfrac,nm)
  ! ----------------------------------------------------------------------
  ! Write a header describing the full setup of the calculation
  ! CC is the comment character that should be added in front of each line
  ! ----------------------------------------------------------------------
  use Defs
  implicit none
  integer        :: unit
  integer        :: na,i,nm
  real (kind=dp) :: amin,amax,apow,lmin,lmax,pcore,pmantle,fmax,mfrac(nm),rho_av
  real (kind=dp) :: amean(3)
  character*(*)  :: cc

  call plmeans(amin,amax,apow,amean)
  write(unit,'(A,"============================================================================")') cc
  write(unit,'(A," Opacities computed by OpTool          <a^n>=",1p,3e9.2e1)') cc,amean
  write(unit,'(A," Parameters:")') cc
  write(unit,'(A,"   amin [um]=",f11.3," amax [um]=",f11.3,"  na  =",I5,"    apow=",g10.2)') cc,amin, amax, na, apow
  write(unit,'(A,"   lmin [um]=",f11.3," lmax [um]=",f11.3,"  nlam=",I5,"    nang=",I6)') cc,lmin, lmax, nlam, nang
  write(unit,'(A,"   porosity =",f11.3," p_mantle =",f11.3,"  fmax=",g9.2,"chop=  ",f4.1)') cc,pcore,pmantle,fmax,chopangle
  write(unit,'(A," Composition:")') cc
  write(unit,'(A,"  Where   mfrac  rho   Material")') cc
  write(unit,'(A,"  -----   -----  ----  ----------------------------------------------------")') cc
  do i=1,nm
     write(unit,'(A,"  ",A6,f7.3,f6.2,"  ",A)') cc,mat_loc(i), mfrac(i)/sum(mfrac(1:nm)),mat_rho(i),trim(mat_lnk(i))
  enddo
  if (rho_av .gt. 0.d0) then
     write(unit,'(A,"  - - -   - - -  -  -  - - - - - - - - - - - - - - - - - - - - - - - - -  -")') cc
     if ( (pcore+pmantle) .gt. 0.d0) then
        write(unit,'(A,"  ",A6,f7.3,f6.2,"  ","mixture of",i3," materials and vacuum")') cc,'grain  ', 1.0,rho_av,nm
     else
        write(unit,'(A,"  ",A6,f7.3,f6.2,"  ","mixture of",i3," materials")') cc,'grain  ', 1.0,rho_av,nm
     endif
  endif
  write(unit,'(A,"============================================================================")') cc
end subroutine write_header

subroutine write_ascii_file(p,amin,amax,apow,na,lmin,lmax,fmax,pcore,pmantle,&
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
  real (kind=dp) :: amin,amax,apow,fmax,pcore,pmantle,mfrac(nm)
  real (kind=dp) :: lmin,lmax,f
  type(particle) :: p
  integer        :: na,i,ilam,iang,nm,i1,i2
  real (kind=dp) :: mu1,mu2,dmu,theta1,theta2,tot
  real (kind=dp),allocatable :: f11(:)
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
          pcore,pmantle,p%rho,fmax,mfrac,nm)
     write(20,'("# Output file formatted for RADMC-3D, dustkappa, no scattering matrix")')
     write(20,'("#    iformat")')
     write(20,'("#    nlambda")')
     write(20,'("#    lambda[um]  kabs [cm^2/g]  ksca [cm^2/g]    g_asymmetry")')
     write(20,'("#============================================================================")')
     write(20,*) 3  ! iformat
     write(20,*) nlam  ! number of lambda points
     do ilam=1,nlam
        write(20,'(1p,e15.6,1p,e15.6,1p,e15.6,1p,e15.6)') lam(ilam),p%Kabs(ilam),p%Ksca(ilam),p%g(ilam)
     enddo
     close(20)

  else

     if (progress) write(*,'("Writing full scattering data to file: ",A)') trim(file2)
     open(20,file=file2,RECL=100000)
     call write_header(20,'#',amin,amax,apow,na,lmin,lmax, &
          pcore,pmantle,p%rho,fmax,mfrac,nm)
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

     write(20,*) nlam   ! Number of wavelength points
     if (for_radmc) then
        write(20,*) nang+1 ! Number of angular points
     else
        write(20,*) nang
     endif
     write(20,*)         ! an empty line
     ! The opacities as function of lambda
     do ilam=1,nlam
        write(20,'(1p,e15.6,1p,e15.6,1p,e15.6,1p,e15.6)') lam(ilam),p%Kabs(ilam),p%Ksca(ilam),p%g(ilam)
     enddo

     write(20,*)      ! an empty line
     ! The angular grid
     if (for_radmc) then
        do iang=0,nang
           write(20,'(f11.5)') dble(iang)/dble(nang)*180.d0
        enddo
     else
        do iang=1,nang
           write(20,'(f11.5)') dble(iang-0.5)/dble(nang)*180.d0
        enddo
     endif
        
     write(20,*)      ! an empty line
     ! Write the scattering matrix
     if (for_radmc) allocate(f11(0:nang))
     do ilam=1,nlam
        if (for_radmc) then
           ! Extra renormalization because of interpolation
           do iang=0,nang
              i1 = max(1,iang); i2=min(iang+1,nang)
              f11(iang) = 0.5d0 * ( p%F(ilam)%F11(i1) + p%F(ilam)%F11(i2) )
           enddo
           tot = 0.d0
           do iang=1,nang
              theta1 = dble(iang-1)/dble(nang) * pi; theta2 = dble(iang)  /dble(nang) * pi
              mu1 = cos(theta1); mu2 = cos(theta2); dmu = mu1-mu2
              tot = tot + 0.5d0 * (f11(iang-1)+f11(iang)) * dmu
           enddo
           tot = 2.d0 * pi * tot
           f = p%ksca(ilam)/tot
           ! print *,tot,p%ksca(ilam),f
           do iang=0,nang
              ! We have only computed 0..nang-1, but RADMC needs a value at
              ! 180 degrees as well. We simply repeat the last value
              i1 = max(1,iang); i2=min(iang+1,nang)
              write(20,'(1p,e15.6,1p,e15.6,1p,e15.6,1p,e15.6,1p,e15.6,1p,e15.6)') &
                   0.5d0 * f * ( p%F(ilam)%F11(i1) + p%F(ilam)%F11(i2) ), &
                   0.5d0 * f * ( p%F(ilam)%F12(i1) + p%F(ilam)%F12(i2) ), &
                   0.5d0 * f * ( p%F(ilam)%F22(i1) + p%F(ilam)%F22(i2) ), &
                   0.5d0 * f * ( p%F(ilam)%F33(i1) + p%F(ilam)%F33(i2) ), &
                   0.5d0 * f * ( p%F(ilam)%F34(i1) + p%F(ilam)%F34(i2) ), &
                   0.5d0 * f * ( p%F(ilam)%F44(i1) + p%F(ilam)%F44(i2) )
           enddo
        else
           do i=1,nang
              write(20,'(1p,e15.6,1p,e15.6,1p,e15.6,1p,e15.6,1p,e15.6,1p,e15.6)') &
                   p%F(ilam)%F11(i),p%F(ilam)%F12(i),p%F(ilam)%F22(i), &
                   p%F(ilam)%F33(i),p%F(ilam)%F34(i),p%F(ilam)%F44(i)
           enddo
        endif
     enddo ! end do ilam=1,nlam
     if (for_radmc) deallocate(f11)
     close(20)
  endif   ! end if (not scatter)
end subroutine write_ascii_file

#ifdef USE_FITSIO
subroutine write_fits_file(p,amin,amax,apow,na, &
     fmax,pcore,pmantle, &
     nm,mfrac,rho,fitsfile)
  ! ----------------------------------------------------------------------
  ! Routine to write a FITS file with all the information about
  ! the opacities and scattering properties.
  ! FIXME: Something goes wrong with the keywords
  ! ----------------------------------------------------------------------
  use Defs
  implicit none
  real (kind=dp) :: amin,amax,apow,fmax,pcore,pmantle
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

  simple = .true.
  extend = .true.
  group  = 1
  fpixel = 1

  bitpix    = -64
  naxis     = 2
  naxes(1)  = nlam
  naxes(2)  = 4
  nelements = naxes(1)*naxes(2)
  allocate(array(nlam,4,1))

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

  call ftpkye(unit,'porosity',real(pcore),8,'[g/cm^3]',status)
  call ftpkye(unit,'p_mantle',real(pmantle),8,'[g/cm^3]',status)

  do i=1,nm
     write(word,'("file",i0.2)') i
     call ftpkys(unit,word,trim(mat_lnk(i)),'',status)
  enddo
  do i=1,nm
     write(word,'("frac",i0.2)') i
     call ftpkye(unit,word,real(mfrac(i)),8,'[mass fraction]',status)
  enddo
  do i=1,nm
     write(word,'("rho",i0.2)') i
     call ftpkye(unit,word,real(rho(i)),8,'[g/cm^3]',status)
  enddo

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
  enddo

  call ftpprd(unit,group,fpixel,nelements,array(1:nlam,1:4,1),status)

  deallocate(array)

  ! ----------------------------------------------------------------------
  ! HDU 1: scattering matrix
  ! ----------------------------------------------------------------------
  bitpix    = -64
  naxis     = 3
  naxes(1)  = nlam
  naxes(2)  = 6
  naxes(3)  = nang
  nelements = naxes(1)*naxes(2)*naxes(3)

  allocate(array(nlam,6,nang))

  ! create new hdu
  call ftcrhd(unit, status)

  !  Write the required header keywords.
  call ftphpr(unit,simple,bitpix,naxis,naxes,0,1,extend,status)

  do i=1,nlam
     do j=1,nang
        array(i,1,j) = p%F(i)%F11(j)
        array(i,2,j) = p%F(i)%F12(j)
        array(i,3,j) = p%F(i)%F22(j)
        array(i,4,j) = p%F(i)%F33(j)
        array(i,5,j) = p%F(i)%F34(j)
        array(i,6,j) = p%F(i)%F44(j)
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
  real (kind=dp) :: a1,a2,p,amean(3),e1,e2,e3,numerator,denominator
  integer :: n
  if (abs((a2-a1)/a1) .lt. 1d-6) then
     amean(1) = a1; amean(2) = a1; amean(3) = a1
  else 
     e2 = 1.d0-p
     if (abs(e2).lt.1e-3) then
        denominator = log(a2)-log(a1)
     else
        denominator = (a2**e2-a1**e2)/e2
     endif
     do n=1,3
        e1 = 1.d0-p+n
        e3 = 1.d0/n
        if (abs(e1).lt.1e-3) then
           numerator = log(a2)-log(a1)
        else
           numerator  = (a2**e1-a1**e1)/e1
        endif
        amean(n) = (numerator/denominator)**e3
     enddo
  endif
end subroutine plmeans


!!! **** File Variables

! Local Variables:
! eval: (outline-minor-mode)
! outline-regexp: "!!!\\|\\(program\\|subroutine\\|function\\|module\\)\\>"
! outline-heading-alist: (("!!!" . 1) ("program" . 2) ("subroutine" . 2) ("module" . 2) ("function" . 2))
! End:
