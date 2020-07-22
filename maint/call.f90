subroutine optool_run (lambda,n_lam,n_ang,a_min,a_max,n_a,a_pow,n_m,m_lnk,m_frac,m_loc,kappa,g_scat,f_scat,fmax,p_core,p_mantle)
  use Defs
  integer        intent(in)          :: n_lam,n_ang,n_m
  integer        intent(in),optional :: n_a
  real (kind=dp),intent(in)          :: lambda(n_lam),m_frac(n_m)
  real (kind=dp),intent(in),optional :: a_min,a_max,a_pow
  real (kind=dp),intent(out)         :: kappa(n_lam,3),g(n_lam),f(n_lam,6,n_ang)
  character*(*) ,intent(in)          :: lnk,loc

  real (kind=dp) :: amin,amax,apow
  integer        :: na

  ! Check for the presence of arguments and supply defaults
  if (present(n_a)) then
     na = n_a
  else
     na = 50
  endif
  if (present(a_min)) then
     amin = a_min
  else
     amin = 0.05
  endif
  if (present(a_max)) then
     amax = a_max
  else
     amax = 0.05
  endif
  if (present(a_pow)) then
     apow = a_pow
  else
     apow = 3.5
  endif
  if (present(fmax)) then
     fmax = f_max
  else
     fmax = 0.8
  endif
  if (present(p_core)) then
     pcore = p_core
  else
     pcore = 0d0
  endif
  if (present(p_mantle)) then
     pmantle = p_mantle
  else
     pmantle = pcore
  endif
  
  ! Fill the stuff that is in the Defs Module
  allocate(lam(n_lam))
  lam(1:n_lam) = lambda(1:n_lam)
  nang = n_ang

  ! allocate the material arrays
  allocate(mat_loc(n_m+1))
  allocate(mat_lnk(n_m+1))
  allocate(mat_mfr(n_m+1))
  allocate(mat_rho(n_m+1))
  have_mantle = .false.
  icore = 0
  do im = 1,n_m
     if (m_loc .eq. "core") then
        icore = icore+1
        iat = icore
     else
        ! mantle material
        if (have_mantle) then
           print *,"ERROR only one mantle material allowed"
           stop
        endif
        iat = nm
     endif
     mat_lnk(iat)  = m_lnk(i)
     mat_loc(iat)  = m_loc(i)
     mat_mfr(iat)  = m_frac(i)
  enddo

  ! Allocate all lambda-dependant arrays
  allocate(p%Kabs(nlam))
  allocate(p%Kext(nlam))
  allocate(p%Ksca(nlam))
  allocate(p%g(nlam))
  allocate(p%F(nlam))
  do i=1,nlam
     allocate(p%F(i)%F11(nang),p%F(i)%F12(nang),p%F(i)%F22(nang))
     allocate(p%F(i)%F33(nang),p%F(i)%F34(nang),p%F(i)%F44(nang))
  enddo
  
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
  
  ! call computepart
  call ComputePart(p,amin,amax,apow,na,fmax,pcore,pmantle,mat_mfr,nm,.false.)

  ! copy results into the correct arrays to return them throught the interface
  kappa(1:nlam,1) = p%kabs(1:nlam)
  kappa(1:nlam,2) = p%ksca(1:nlam)
  kappa(1:nlam,3) = p%ext(1:nlam)
  g_scat = p:g(1:nlam)
  do i=1,nlam
     f_scat(i,1,1:nang) = p%F(i)%F11(1:nang)
     f_scat(i,2,1:nang) = p%F(i)%F12(1:nang)
     f_scat(i,3,1:nang) = p%F(i)%F22(1:nang)
     f_scat(i,4,1:nang) = p%F(i)%F33(1:nang)
     f_scat(i,5,1:nang) = p%F(i)%F34(1:nang)
     f_scat(i,6,1:nang) = p%F(i)%F44(1:nang)
  enddo
  
end subroutine optool_run
