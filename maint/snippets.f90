  ! ----------------------------------------------------------------------
  ! F11 normalization with fixed 180
  ! ----------------------------------------------------------------------
  ! Check the normalization of F11
  tot  = 0d0; tot2 = 0d0
  do j=1,n_ang
     !                  F11                      (sin phi)                      d phi        2pi
     tot  = tot +  p%F(ilam)%F11(j) * (sin(pi*(real(j)-0.5)/real(n_ang))) * (pi/180.d0) * (2.d0*pi)
     tot2 = tot2 + sin(pi*(real(j)-0.5)/real(n_ang))*pi/180.d0*2.d0*pi
  enddo
  ! print *,'theta integration gives: ',tot,tot/4.d0/pi
  
  ! do it again, using mu integration
  tot  = 0d0; tot2 = 0d0
  do j=1,n_ang-1
     !                  F11                      (sin phi)                      d phi        2pi
     tot  = tot +  2.d0*pi*p%F(ilam)%F11(j) * (cos(pi*(real(j)/real(n_ang))) - cos(pi*(real(j-1)/real(n_ang))))
  enddo
  ! print *,'mu    integration gives: ',tot,tot/4.d0/pi



  
     ! ----------------------------------------------------------------------
     ! F11 normalization after nang variable
     ! ----------------------------------------------------------------------
     ! Check the normalization of F11
     tot  = 0d0; tot2 = 0d0
     do j=1,nang
        !                  F11                      (sin phi)                      d phi        2pi
        tot  = tot  +  p%F(ilam)%F11(j) * (sin(pi*(real(j)-0.5)/real(nang))) * (pi/dble(nang)) * (2.d0*pi)
        tot2 = tot2 +                      sin(pi*(real(j)-0.5)/real(nang))  * (pi/dble(nang)) * (2.d0*pi)
     enddo
     print *,'theta integration gives: ',tot,tot/4.d0/pi
     
     ! do it again, using mu integration
     tot  = 0d0; tot2 = 0d0
     do j=1,nang-1
        !               2pi          F11                                     d mu
        tot  = tot +  2.d0*pi*p%F(ilam)%F11(j) * (cos(pi*(real(j)/real(nang))) - cos(pi*(real(j-1)/real(nang))))
     enddo
     print *,'mu    integration gives: ',tot,tot/4.d0/pi

