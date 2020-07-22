!!! **** Code to compute mean opacities kappa_Rosseland, kappa_Planck (from Kees Dullemond)

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
     temp(itemp) = tmin * (tmax/tmin)**((itemp-1.0d0)/(ntemp-1.0d0))
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


function bplanck(temp,nu)
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


function bplanckdt(temp,nu)
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

