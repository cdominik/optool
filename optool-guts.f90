subroutine MeerhoffMie(rad,lam,e1,e2,Csca,Cext,F11,F12,F33,F34,na)
  implicit none
  integer na,nangle
  real*8 rad,lam,e1,e2,csca,cext,F11(na),F12(na),F33(na),F34(na)

  integer nparts,develop,nsubr(1),ngaur(1),idis(1),ndis,i
  real*8 delta,cutoff,thmin,thmax,step,wavel(1),Rem(1),fImm(1)
  real*8 par1(1),par2(1),par3(1),rdis(1,300),nwrdis(1,300)
  real*8 F(4,6000),theta(6000)

  nangle  = na
  nparts  = 1
  develop = 0
  delta   = 1d-8
  cutoff  = 1d-8
  thmin   = 180d0*(real(1)-0.5)/real(na)
  thmax   = 180d0*(real(na)-0.5)/real(na)
  step    = (thmax-thmin)/real(na-1)
  wavel   = lam
  Rem     = e1
  fImm    = e2
  nsubr   = 1
  ngaur   = 1
  idis    = 0
  par1    = rad
  par2    = 0d0
  par3    = 0d0
  ndis    = 1
  rdis(1,1)  =rad
  nwrdis(1,1)=1d0

  call mie(nparts, develop, delta, cutoff, thmin, thmax, step, &
       wavel, Rem, fImm, nsubr , ngaur, idis, &
       par1, par2, par3, ndis, rdis, nwrdis, &
       nangle, theta, F ,Cext, Csca)

  do i=1,nangle
     F11(i) = F(1,nangle-i+1)
     F12(i) = F(2,nangle-i+1)
     F33(i) = F(3,nangle-i+1)
     F34(i) = F(4,nangle-i+1)
  enddo

  return
end subroutine MeerhoffMie


subroutine mie(nparts, develop, delta, cutoff, thmin, thmax, step, &
     wavel, Rem, fImm, nsubr , ngaur, idis, &
     par1, par2, par3, ndis, rdis, nwrdis, &
     nangle, theta, F ,outCext, outCsca)

  implicit double precision (a-h,o-z)
  implicit integer (i-n)
  parameter ( NDang=6000, NDcoef=6000, NDpart = 4, nrunit=88 )
  parameter ( NDdis=300 )
  double precision lambda, nr, ni, miec, nwrdis,outCext,outCsca
  integer develop
  double complex m
  character*60 s1, s2, s3, s4, s5, s6
  character*10 scfile
  dimension F(4,NDang), &
       miec(13),u(NDang),w(NDang),coefs(4,4,0:NDcoef), &
       wavel(*), Rem(*), fImm(*), nsubr(*), &
       ngaur(*), idis(*), &
       par1(*), par2(*), par3(*), &
       rdis(1,*), nwrdis(1,*), &
       scfile(NDpart), theta(NDang)
  pi     = dacos(-1.D0)
  radtod = 180.D0/pi
  scfile(1) = 'mie1_sc  '
  scfile(2) = 'mie2_sc  '
  scfile(3) = 'mie3_sc  '
  scfile(4) = 'mie4_sc  '
  !************************************************************************
  !*  Start a loop over different 'particles' referring to the different  *
  !*  cases specified in the input                                        *
  !************************************************************************
  do iparts=1,nparts
     !************************************************************************
     !*  Determine the integration bounds for the radius integration         *
     !************************************************************************
     call rminmax( idis(iparts), nsubr(iparts), ngaur(iparts), &
          par1(iparts), par2(iparts), par3(iparts), cutoff, &
          rmin, rmax, &
          iparts, ndis, rdis, nwrdis )
     
     lambda = wavel(iparts)
     nr     = Rem(iparts)
     ni     = fImm(iparts)
     m      = dcmplx(nr,-ni)
     !************************************************************************
     !*  Calculate the scattering matrix : this is most of the work !        *
     !************************************************************************
     call scatmat( u, w, m, wavel(iparts), idis(iparts), &
          thmin, thmax, step, develop, &
          nsubr(iparts), ngaur(iparts), rmin, rmax, &
          par1(iparts), par2(iparts), par3(iparts), &
          delta, F, miec, nangle, &
          ndis, rdis, nwrdis, iparts )
     !************************************************************************
     !*  Test if the number of coefficients needed fits in the array coefs   *
     !*  The number of coefficients is equal to the number of integration    *
     !*  points in scattering angle !! (think about this !)                  *
     !************************************************************************
     ncoef = nangle
     if ((ncoef .gt. NDcoef) .and. (develop .eq. 1)) then
        write(*,*) ' main: too many coefficients needed :',ncoef
        write(*,*) '       maximum array size NDcoef = ',NDcoef
        write(*,*) '       Therefore I cannot do the expansion.'
        develop=0
        call clossc(nrunit)
     endif
     !************************************************************************
     !*  If required, calculate the expansion coefficients                   *
     !************************************************************************
     if (develop.eq.1) then
        call devel(ncoef,nangle,u,w,F,coefs)
        npunt = idnint((thmax-thmin)/step)+1
        if (npunt .gt. NDang) then
           write(*,*) ' main: too many angles for table npunt=',npunt
           write(*,*) '       maximum array size NDang = ',NDang
           npunt=61
           write(*,*) '      I have set the number of angles to ',npunt
        endif
     endif
     outCsca=miec(1)
     outCext=miec(2)

     if (nangle .ge. 1) then
        cosbar = coefs(1,1,1)/3.D0
     else
        cosbar = 0.D0
     endif

     do j=1,nangle
        i     = nangle-j+1
        xp    = -100.d0*F(2,i)/F(1,i)
        theta1 = radtod*dacos(u(i))
     enddo
     !************************************************************************
     !*  Print the coefficients on output.                                   *
     !*  Evaluate the expansion in GSF at an equidistant set of scattering   *
     !*  angles and print the resulting scattering matrix.                   *
     !************************************************************************
     if (develop.eq.1) then
        call expand( ncoef, npunt, coefs, u, F, thmin, thmax, step )
        do  i=1,npunt
           xp    = -100.d0*F(2,i)/F(1,i)
           theta1 = radtod*dacos(u(i))
        enddo
     endif
  enddo
  !************************************************************************
  !*  End of loop over 'particles'.                                       *
  !************************************************************************
  return
end subroutine mie

subroutine anbn( m, x, nmax, psi, chi, d, an, bn )
  !************************************************************************
  !*  Calculate the Mie coefficients an and bn.                           *
  !************************************************************************
  implicit double precision (a-h,o-z)
  parameter( NDn=10000000 )
  double complex m, zn, znm1, save, perm
  double complex an(nmax), bn(nmax), D(nmax)
  dimension psi(0:nmax), chi(0:nmax)
  perm = 1.D0/m
  perx = 1.D0/x
  xn   = 0.D0
  do n=1, nmax
     zn   = dcmplx( psi(n),   chi(n))
     znm1 = dcmplx( psi(n-1), chi(n-1))
     xn   = dble(n)*perx
     save = D(n)*perm+xn
     an(n)= (save*psi(n)-psi(n-1)) / (save*zn-znm1)
     save = m*D(n)+xn
     bn(n)= (save*psi(n)-psi(n-1)) / (save*zn-znm1)
  enddo
  return
end subroutine anbn


subroutine clossc(iunit)
  !************************************************************************
  !*  On entry :                                                          *
  !*      iunit     number of the unit to be closed                       *
  !*  On exit :                                                           *
  !*      The file is closed.                                             *
  !*                                                 VLD January 9, 1989  *
  !************************************************************************
  !c      close(unit=iunit,err=999)
  return
  !c  999 write(7,*) ' clossc: error in closing file with unit number',iunit
  !c      stop 'in clossc error in closing file'
end subroutine clossc

subroutine devel(ncoef,nangle,u,w,F,coefs)
  !************************************************************************
  !*  Calculate the expansion coefficients of the scattering matrix in    *
  !*  generalized spherical functions by numerical integration over the   *
  !*  scattering angle.                                                   *
  !************************************************************************
  implicit double precision (a-h,o-z)
  parameter( NDcoef=3000, NDang = 6000 )
  dimension u(NDang),w(NDang),F(4,NDang)
  dimension coefs(4,4,0:NDcoef), P00(NDang,2), P02(NDang,2) &
       , P22(NDang,2), P2m2(NDang,2)
  !************************************************************************
  !*  Initialization                                                      *
  !************************************************************************
  qroot6 = -0.25D0*dsqrt(6.D0)
  do j=0,NDcoef
     do i=1,4
        do ii=1,4
           coefs(ii,i,j)=0.D0
        enddo
     enddo
  enddo
  !************************************************************************
  !*  Multiply the scattering matrix F with the weights w for all angles  *
  !*  We do this here because otherwise it should be done for each l      *
  !************************************************************************
  do k=1,4
     do i=1, nangle
        F(k,i) = w(i)*F(k,i)
     enddo
  enddo

  !************************************************************************
  !*  Start loop over the coefficient index l                             *
  !*  first update generalized spherical functions, then calculate coefs. *
  !*  lold and lnew are pointer-like indices used in recurrence           *
  !************************************************************************
  lnew = 1
  lold = 2

  do l=0, ncoef
     if (l .eq. 0) then
        !************************************************************************
        !*             Adding paper Eq. (77) with m=0                           *
        !************************************************************************
        do i=1, nangle
           P00(i,lold)  = 1.D0
           P00(i,lnew)  = 0.D0
           P02(i,lold)  = 0.D0
           P22(i,lold)  = 0.D0
           P2m2(i,lold) = 0.D0
           P02(i,lnew)  = 0.D0
           P22(i,lnew)  = 0.D0
           P2m2(i,lnew) = 0.D0
        enddo
     else
        fac1 = (2.D0*l-1.d0)/dble(l)
        fac2 = dble(l-1.d0)/dble(l)
        !************************************************************************
        !*             Adding paper Eq. (81) with m=0                           *
        !************************************************************************
        do i=1, nangle
           P00(i,lold) = fac1*u(i)*P00(i,lnew) - fac2*P00(i,lold)
        enddo
     endif
     if (l .eq. 2) then
        !************************************************************************
        !*             Adding paper Eqs. (78) and (80)  with m=2                *
        !*             sql4 contains the factor dsqrt(l*l-4) needed in          *
        !*             the recurrence Eqs. (81) and (82)                        *
        !************************************************************************
        do i=1, nangle
           P02(i,lold)  = qroot6*(1.D0-u(i)*u(i))
           P22(i,lold)  = 0.25D0*(1.D0+u(i))*(1.D0+u(i))
           P2m2(i,lold) = 0.25D0*(1.D0-u(i))*(1.D0-u(i))
           P02(i,lnew)  = 0.D0
           P22(i,lnew)  = 0.D0
           P2m2(i,lnew) = 0.D0
        enddo
        sql41 = 0.D0
     else if (l .gt. 2) then
        !************************************************************************
        !*             Adding paper Eq. (82) with m=0 and m=2                   *
        !************************************************************************
        sql4  = sql41
        sql41 = dsqrt(dble(l*l)-4.d0)
        twol1 = 2.D0*dble(l)-1.d0
        tmp1  = twol1/sql41
        tmp2  = sql4/sql41
        denom = (dble(l)-1.d0)*(dble(l*l)-4.d0)
        fac1  = twol1*(dble(l)-1.d0)*dble(l)/denom
        fac2  = 4.D0*twol1/denom
        fac3  = dble(l)*((dble(l)-1.d0)*(dble(l)-1.d0)-4.d0)/denom
        do i=1, nangle
           P02(i,lold) = tmp1*u(i)*P02(i,lnew) - tmp2*P02(i,lold)
           P22(i,lold) = (fac1*u(i)-fac2)*P22(i,lnew) &
                - fac3*P22(i,lold)
           P2m2(i,lold)= (fac1*u(i)+fac2)*P2m2(i,lnew) &
                - fac3*P2m2(i,lold)
        enddo
     endif
     !************************************************************************
     !*         Switch indices so that lnew indicates the function with      *
     !*         the present index value l, this mechanism prevents swapping  *
     !*         of entire arrays.                                            *
     !************************************************************************
     itmp = lnew
     lnew = lold
     lold = itmp
     !************************************************************************
     !*         Now calculate the coefficients by integration over angle     *
     !*         See de Haan et al. (1987) Eqs. (68)-(73).                    *
     !*         Remember for Mie scattering : F11 = F22 and F33 = F44        *
     !************************************************************************
     alfap = 0.D0
     alfam = 0.D0
     do i=1, nangle
        coefs(1,1,l) = coefs(1,1,l) + P00(i,lnew)*F(1,i)
        alfap = alfap + P22(i,lnew)*(F(1,i)+F(3,i))
        alfam = alfam + P2m2(i,lnew)*(F(1,i)-F(3,i))
        coefs(4,4,l) = coefs(4,4,l) + P00(i,lnew)*F(3,i)
        coefs(1,2,l) = coefs(1,2,l) + P02(i,lnew)*F(2,i)
        coefs(3,4,l) = coefs(3,4,l) + P02(i,lnew)*F(4,i)
     enddo
     !************************************************************************
     !*         Multiply with trivial factors like 0.5D0*(2*l+1)             *
     !************************************************************************
     fl = dble(l)+0.5D0
     coefs(1,1,l) =  fl*coefs(1,1,l)
     coefs(2,2,l) =  fl*0.5D0*(alfap+alfam)
     coefs(3,3,l) =  fl*0.5D0*(alfap-alfam)
     coefs(4,4,l) =  fl*coefs(4,4,l)
     coefs(1,2,l) =  fl*coefs(1,2,l)
     coefs(3,4,l) =  fl*coefs(3,4,l)
     coefs(2,1,l) =     coefs(1,2,l)
     coefs(4,3,l) =    -coefs(3,4,l)
  enddo
  !************************************************************************
  !*     End of loop over index l                                         *
  !************************************************************************

  !************************************************************************
  !*     Remove the weight factor from the scattering matrix              *
  !************************************************************************
  do k=1, 4
     do i=1, nangle
        F(k,i) = F(k,i)/w(i)
     enddo
  enddo
  
  return
end subroutine devel

subroutine expand(ncoef,nangle,coefs,u,F,thmin,thmax,step)

  !     ***********************************************************

  implicit double precision (a-h,o-z)
  parameter( pi=3.141592653589793238462643D0, radfac= pi/180.D0 )
  parameter( NDcoef=3000, NDang=6000 )
  dimension u(NDang),F(4,NDang)
  dimension coefs(4,4,0:NDcoef), P00(NDang,2), P02(NDang,2)
  !
  !     **********************************************************
  !
  !                     initialize

  do i=1, nangle
     u(i) = dcos(radfac*(thmin+dble(i-1)*step))
  enddo
  qroot6 = -0.25D0*dsqrt(6.D0)
  !************************************************************************
  !*  Set scattering matrix F to zero                                     *
  !************************************************************************
  do k=1,4
     do i=1, nangle
        F(k,i) = 0.D0
     enddo
  enddo
  !************************************************************************
  !*  Start loop over the coefficient index l                             *
  !*  first update generalized spherical functions, then calculate coefs. *
  !*  lold and lnew are pointer-like indices used in recurrence           *
  !************************************************************************
  lnew = 1
  lold = 2
  do l=0, ncoef
     if (l .eq. 0) then
        !************************************************************************
        !*             Adding paper Eqs. (76) and (77) with m=0                 *
        !************************************************************************
        do i=1, nangle
           P00(i,lold) = 1.D0
           P00(i,lnew) = 0.D0
           P02(i,lold) = 0.D0
           P02(i,lnew) = 0.D0
        enddo
     else
        fac1 = (2.D0*dble(l)-1.d0)/dble(l)
        fac2 = (dble(l)-1.d0)/dble(l)
        !************************************************************************
        !*             Adding paper Eq. (81) with m=0                           *
        !************************************************************************
        do i=1, nangle
           P00(i,lold) = fac1*u(i)*P00(i,lnew) - fac2*P00(i,lold)
        enddo
     endif
     if (l .eq. 2) then
        !************************************************************************
        !*             Adding paper Eq. (78)                                    *
        !*             sql4 contains the factor dsqrt((l+1)*(l+1)-4) needed in  *
        !*             the recurrence Eqs. (81) and (82)                        *
        !************************************************************************
        do i=1, nangle
           P02(i,lold) = qroot6*(1.D0-u(i)*u(i))
           P02(i,lnew) = 0.D0
        enddo
        sql41 = 0.D0
     else if (l .gt. 2) then
        !************************************************************************
        !*             Adding paper Eq. (82) with m=0                           *
        !************************************************************************
        sql4  = sql41
        sql41 = dsqrt(dble(l*l)-4.d0)
        tmp1  = (2.D0*dble(l)-1.d0)/sql41
        tmp2  = sql4/sql41
        do i=1, nangle
           P02(i,lold) = tmp1*u(i)*P02(i,lnew) - tmp2*P02(i,lold)
        enddo
     endif
     !************************************************************************
     !*         Switch indices so that lnew indicates the function with      *
     !*         the present index value l, this mechanism prevents swapping  *
     !*         of entire arrays.                                            *
     !************************************************************************
     itmp = lnew
     lnew = lold
     lold = itmp
     !************************************************************************
     !*         Now add the l-th term to the scattering matrix.              *
     !*         See de Haan et al. (1987) Eqs. (68)-(73).                    *
     !*         Remember for Mie scattering : F11 = F22 and F33 = F44        *
     !************************************************************************
     do i=1, nangle
        F(1,i) = F(1,i) + coefs(1,1,l)*P00(i,lnew)
        F(2,i) = F(2,i) + coefs(1,2,l)*P02(i,lnew)
        F(3,i) = F(3,i) + coefs(4,4,l)*P00(i,lnew)
        F(4,i) = F(4,i) + coefs(3,4,l)*P02(i,lnew)
     enddo
  enddo
  !************************************************************************
  !*     End of loop over index l                                         *
  !************************************************************************
  return
end subroutine expand
   
subroutine fichid( m, x, nchi, nmax, nD, psi, chi, D )
  !************************************************************************
  !*  Calculate functions psi(x)  chi(x) and D(z) where z = mx.           *
  !*  On entry, the following should be supplied :                        *
  !*      m      : complex index of refraction                            *
  !*      x      : sizeparameter                                          *
  !*      nchi   : starting index for backwards recurrence of chi         *
  !*      nmax   : number of chi, psi and D that must be available        *
  !*      nD     : starting index for backwards recurrence of D           *
  !*  On exit, the desired functions are returned through psi, chi and D  *
  !************************************************************************
  implicit double precision (a-h,o-z)
  double complex D,m,z, perz, zn1
  dimension psi(0:nchi), chi(0:nmax+1), D(nd)

  z = m*x
  perz = 1.D0/z
  perx = 1.D0/x
  sinx = dsin(x)
  cosx = dcos(x)
  psi(nchi)=0d0
  !************************************************************************
  !*  (mis-) use the array psi to calculate the functions rn(x)
  !*  De Rooij and van der Stap Eq. (A6)
  !************************************************************************
  do n=nchi-1, 0, -1
     psi(n) = 1.D0 / (dble(2*n+1)/x - psi(n+1))
  enddo
  !************************************************************************
  !*  Calculate functions D(z) by backward recurrence
  !*  De Rooij and van der Stap Eq. (A11)
  !************************************************************************
  D(nD) = dcmplx(0.D0,0.D0)
  do n=nD - 1, 1, -1
     zn1 = dble(n+1)*perz
     D(n) = zn1 - 1.D0/(D(n+1)+zn1)
  enddo
  !************************************************************************
  !*  De Rooij and van der Stap Eqs. (A3) and (A1)
  !*  De Rooij and van der Stap Eq. (A5) where psi(n) was equal to r(n)
  !*  and Eq. (A2)
  !************************************************************************
  psi(0) = sinx
  psi1   = psi(0)*perx - cosx
  if (dabs(psi1) .gt. 1.d-4) then
     psi(1) = psi1
     do n=2,nmax
        psi(n) = psi(n)*psi(n-1)
     enddo
  else
     do n=1,nmax
        psi(n) = psi(n)*psi(n-1)
     enddo
  endif
  !************************************************************************
  !*  De Rooij and van der Stap Eqs. (A4) and (A2)
  !************************************************************************
  chi(0) = cosx
  chi(1) = chi(0)*perx + sinx
  do n=1, nmax-1
     chi(n+1) = dble(2*n+1)*chi(n)*perx - chi(n-1)
  enddo
  !* DEBUG
  !  write(7,*) ' fichid: x = ',x
  !  write(7,12)
  !  do n=0, nchi
  !     write(7,11) n,psi(n),chi(n),D(n)
  !  enddo
  !11 format(i4,4e24.14)
  !12 format('   n',t20,'psi(n)',t44,'chi(n)',t68 &
  !        ,'Re(D(n))',t92,'Im(D(n))',/ &
  !        ,' ----------------------------------------------------' &
  !        ,'-----------------------------------------------------')
  !* END DEBUG
  return
end subroutine fichid

subroutine gaulegCP(ndim,ngauss,a,b,x,w)
  !************************************************************************
  !*   Given the lower and upper limits of integration a and b, and given *
  !*   the number of Gauss-Legendre points ngauss, this routine returns   *
  !*   through array x the abscissas and through array w the weights of   *
  !*   the Gauss-Legendre quadrature formula. Eps is the desired accuracy *
  !*   of the abscissas. This routine is documented further in :          *
  !*   W.H. Press et al. 'Numerical Recipes' Cambridge Univ. Pr. (1987)   *
  !*   page 125 ISBN 0-521-30811-9                                        *
  !************************************************************************
  implicit double precision (a-h,o-z)
  parameter( eps = 1.d-14 )
  double precision x(ndim), w(ndim)
  double precision a,b,xm,xl,z,p1,p2,p3,pp,z1,pi
  pi = 4.D0*datan(1.D0)
  m  = (ngauss+1)/2
  xm = 0.5D0*(a+b)
  xl = 0.5D0*(b-a)
  do i=1,m
     !* THIS IS A REALLY CLEVER ESTIMATE :
     z = dcos(pi*(dble(i)-0.25D0)/(dble(ngauss)+0.5D0))
1    continue
     p1 = 1.D0
     p2 = 0.D0
     do j=1,ngauss
        p3 = p2
        p2 = p1
        p1 = ((dble(2*j)-1.d0)*z*p2-(dble(j)-1.d0)*p3)/dble(j)
     enddo
     pp = ngauss*(z*p1-p2)/(z*z-1.d0)
     z1 = z
     z  = z1-p1/pp
     if (dabs(z-z1) .gt. eps) goto 1
     x(i) = xm-xl*z
     x(ngauss+1-i)= xm+xl*z
     w(i) = 2.D0*xl/((1.D0-z*z)*pp*pp)
     !******          write(7,*) ' gaulegCP: ',i,'  x=',x(i),' w=',w(i)
     w(ngauss+1-i)= w(i)
  enddo
  return
end subroutine gaulegCP

subroutine gausspt(ndim,ngauss,a,b,x,w)
  !
  !     ***********************************************************
  !
  !     put the gauss-legendre points of order ndim in array x,
  !     the weights in array w. the points are transformed to the
  !     interval [a,b]
  !
  !     ***********************************************************
  
  implicit double precision (a-h,o-z)
  dimension x(ndim),w(ndim)

  !     ***********************************************************
  !
  !                     find starting values
  !
  gn     = 0.5D0/dble(ngauss)
  extra  = 1.0D0/(.4D0*dble(ngauss)*dble(ngauss)+5.0D0)
  xz     = -gn
  nt     = 0
  nteken = 0
5 pnm2   = 1.0D0
  pnm1   = xz
  do i=2,ngauss
     pnm1xz = pnm1*xz
     pn     = 2.0D0*pnm1xz-pnm2-(pnm1xz-pnm2)/dble(i)
     pnm2   = pnm1
     pnm1   = pn
  enddo
  mteken=1
  if (pn .le. 0.0D0) mteken=-1
  if ((mteken+nteken) .eq. 0) then
     nt    = nt+1
     x(nt) = xz
  endif
  nteken   = mteken
  if ((1.0D0-xz) .le. extra) go to 30
  xz = xz+(1.D0-xz*xz)*gn+extra
  go to 5
30 continue

  !     ***********************************************************
  !
  !                determine zero's and weights

  do i=1,nt
     xz     = x(i)
     delta2 = 1.D0
35   pnm2   = 1.0D0
     pnm1   = xz
     pnm1af = 1.0D0
     z      = .5D0 + 1.5D0*xz*xz
     do k=2,ngauss
        pnm1xz = pnm1*xz
        pn     = 2.0D0*pnm1xz-pnm2-(pnm1xz-pnm2)/dble(k)
        pnaf   = xz*pnm1af+k*pnm1
        z      = z+(dble(k)+0.5D0)*pn*pn
        pnm2   = pnm1
        pnm1   = pn
        pnm1af = pnaf
     enddo
     delta1    = pn/pnaf
     xz        = xz-delta1
     if(delta1.lt.0.0D0) delta1=-delta1
     if( (delta1.ge.delta2) .and. (delta2.lt.1.d-14) ) go to 50
     delta2 = delta1
     go to 35
50   x(i) = xz
     w(i) = 1.0D0/z
  enddo

  !     ***********************************************************
  !
  !                 transform to the interval [a,b]

  nghalf = ngauss/2
  ngp1   = ngauss+1
  ntp1   = nt+1
  apb    = a+b
  bmag2  = (b-a)/2.0D0
  do i=1,nghalf
     x(ngp1-i) = b-bmag2*(1.0D0-x(ntp1-i))
     w(ngp1-i) = bmag2*w(ntp1-i)
  enddo
  if (nghalf .ne. nt) then
     x(nt) = apb/2.0D0
     w(nt) = w(1)*bmag2
  endif
  do i=1,nghalf
     x(i) = apb-x(ngp1-i)
     w(i) = w(ngp1-i)
  enddo
  return
end subroutine gausspt

subroutine pitau(u,nmax,pi,tau)
  !     ***********************************************************
  !     calculates pi,n(u) and tau,n(u) with upward recursion
  !
  !     ***********************************************************
  implicit double precision (a-h,o-z)
  dimension pi(nmax),tau(nmax)
  !
  !       starting values:
  pi(1) = 1.D0
  pi(2) = 3.D0*u
  delta = 3.D0*u*u-1.d0
  tau(1)= u
  tau(2)= 2.D0*delta-1.d0

  !        upward recursion:
  do n=2, nmax-1
     pi(n+1) = dble(n+1)/dble(n)*delta + u*pi(n)
     delta   = u*pi(n+1) - pi(n)
     tau(n+1)= dble(n+1)*delta - pi(n)
  enddo
  return
end subroutine pitau

subroutine rminmax( idis, nsubr, ngaur, par1, par2, par3, cutoff &
     , rmin, rmax, iparts, ndis, rdis, nwrdis )
  !************************************************************************
  !*  Find the integration bounds rmin and rmax for the integration over  *
  !*  a size distribution. These bounds are chosen such that the size     *
  !*  distribution falls below the user specified cutoff. It is essential *
  !*  that the size distribution is normalized such that the integral     *
  !*  over all r is equal to one !                                        *
  !*  This is programmed rather clumsy and will in the future be changed  *
  !************************************************************************
  implicit double precision (a-h,o-z)
  parameter( eps = 1.d-10, NDpart = 4, NDdis=300 )
  double precision nwithr, rdis(NDpart,NDdis), nwrdis(NDpart,NDdis)
  dimension r(1), nwithr(1)
  !*
  if (idis.eq.0) then
     rmin= par1
     rmax= par1
  else
     goto (10,20,30,40,50,60,70,80,90) idis
     write(*,*) ' rminmax: illegal size distribution index :',idis
     stop 'in rminmax illegal size distribution index'
     !************************************************************************
10   sef = 1.D0/dsqrt(par2+3.D0)
     ref = 1.D0/(sef*sef*par2)
     rref= ref
     goto 100

20   ref = par1
     sef = dsqrt(par2)
     rref= ref
     goto 100

30   sef = dsqrt(par3)
     ref = dmax1(par1,par2)+sef
     rref= dmin1(par1,par2)
     goto 100

40   sef = dsqrt(dexp(dlog(par2)**2)-1.d0)
     ref = par1*(1.D0+sef*sef)**0.4D0
     rref= ref
     goto 100

50   ref = par1
     sef = dsqrt(ref)
     rref= ref
     goto 100

60   rmin= par2
     rmax= par3
     goto 999

70   ref = par2
     sef = 2.D0*ref
     rref=0.5D0*ref
     goto 100

80   ref = (par1/(par2*par3))**par3
     sef = 2.D0*ref
     rref= 0.5D0*ref
     goto 100

90   rmin = par1
     rmax = par2
     goto 999

     !************************************************************************
     !*  search for a value of r such that the size distribution
     !*  is less than the cutoff. Start the search at ref+sef which          *
     !*  guarantees that such a value will be found on the TAIL of the       *
     !*  distribution.                                                       *
     !************************************************************************
100  r(1) = ref+sef
     r0   = ref
200  call sizedis( idis, par1, par2, par3, r, 1, nwithr,  &
          iparts, ndis, rdis, nwrdis )
     if (nwithr(1) .gt. cutoff) then
        r0   = r(1)
        r(1) = 2.D0*r(1)
        goto 200
     endif
     r1 = r(1)
     !************************************************************************
     !*  Now the size distribution assumes the cutoff value somewhere        *
     !*  between r0 and r1  Use bisection to find the corresponding r        *
     !************************************************************************
300  r(1) = 0.5D0*(r0+r1)
     call sizedis( idis, par1, par2, par3, r, 1, nwithr, &
          iparts, ndis, rdis, nwrdis )
     if (nwithr(1) .gt. cutoff) then
        r0 = r(1)
     else
        r1 = r(1)
     endif
     if ((r1-r0) .gt. eps) goto 300
     rmax = 0.5D0*(r0+r1)
     !************************************************************************
     !*  Search for a value of r on the low end of the size distribution     *
     !*  such that the distribution falls below the cutoff. There is no      *
     !*  guarantee that such a value exists, so use an extra test to see if  *
     !*  the search comes very near to r = 0                                 *
     !************************************************************************
     r1 = rref
     r0 = 0.D0
400  r(1) = 0.5D0*r1
     call sizedis( idis, par1, par2, par3, r, 1, nwithr, &
          iparts, ndis, rdis, nwrdis )
     if (nwithr(1) .gt. cutoff) then
        r1 = r(1)
        if (r1 .gt. eps) goto 400
     else
        r0 = r(1)
     endif
     !************************************************************************
     !*  Possibly the size distribution goes through cutoff between r0       *
     !*  and r1 try to find the exact value of r where this happens by       *
     !*  bisection.                                                          *
     !*  In case there is no solution, the algorithm will terminate soon.    *
     !************************************************************************
500  r(1) = 0.5D0*(r0+r1)
     call sizedis( idis, par1, par2, par3, r, 1, nwithr, &
          iparts, ndis, rdis, nwrdis )
     if (nwithr(1) .gt. cutoff) then
        r1 = r(1)
     else
        r0 = r(1)
     endif
     if ((r1-r0) .gt. eps) goto 500
     if(r1 .le. eps) then
        rmin = 0.D0
     else
        rmin = 0.5D0*(r0+r1)
     endif
  endif

999 continue!write(7,*) ' rminmax: found rmin = ',rmin,' rmax = ',rmax
  return
end subroutine rminmax


subroutine scatmat( u, wth, m, lambda, idis, &
     thmin, thmax, step, develop, &
     nsub, ngauss, rmin, rmax, &
     par1, par2, par3, &
     delta, F, miec, nangle, &
     ndis, rdis, nwrdis, iparts )
  !************************************************************************
  !*  Calculate the scattering matrix of an ensemble of homogenous        *
  !*  spheres. On entry, the following must be supplied :                 *
  !*     m            : complex index of refraction                       *
  !*     lambda       : wavelength                                        *
  !*     idis         : index of the size distribution                    *
  !*     nsub         : number of subintervals for integration over r     *
  !*     ngauss       : number of Gauss points used per subinterval       *
  !*     rmin         : lower bound for integration over r                *
  !*     rmax         : upper bound for integration over r                *
  !*     par1,2,3     : parameters of the size distribution               *
  !*     delta        : cutoff used in truncation of the Mie sum          *
  !*     thmin        : minimum scattering angle in degrees               *
  !*     thmax        : maximum scattering angle in degrees               *
  !*     step         : step in scattering angle in degrees               *
  !*     develop      : expansion in GSF (1) or not (0)                   *
  !*  On exit, the following results are returned :                       *
  !*     u            : cosines of scattering angles                      *
  !*     wth          : Gaussian weights associated with u                *
  !*     F            : scattering matrix for all cosines in u            *
  !*     miec         : array containing cross sections etc.              *
  !*     nangle       : the number of scattering angles                   *
  !*  When develop=1 u contains Gauss points, when develop=0 u contains   *
  !*  cosines of equidistant angles between thmin and thmax with step.    *
  !************************************************************************
  implicit double precision (a-h,o-z)
  parameter( NDn=10000000, NDr=1000, NDang=6000, NDdis=300,NDpart=4)
  double complex   m, ci, Splusf, Sminf, cSplusf
  double complex   cSminf, Splusb, Sminb, cSplusb, cSminb
  
  double complex, allocatable :: an(:), bn(:),D(:)
  double precision,allocatable :: pi(:), tau(:), fi(:), chi(:)
  double precision,allocatable :: facf(:), facb(:)
  
  double precision lambda, nwithr, miec, numpar, thmin, thmax, step, nwrdis
  integer     develop
  dimension   u(NDang), wth(NDang), F(4,NDang),rdis(NDpart, NDdis)
  dimension   miec(13), nwrdis(NDpart, NDdis),r(NDr), w(NDr), nwithr(NDr)
  logical     symth
  !************************************************************************
  !*  Initialization                                                      *
  !************************************************************************
  do j=1,NDang
     do k=1,4
        F(k,j)=0.D0
     enddo
  enddo
  Csca  = 0.D0
  Cext  = 0.D0
  numpar= 0.D0
  G     = 0.D0
  reff  = 0.D0
  nfou  = 0
  fac90 = 1.D0
  ci    = dcmplx(0.D0,1.D0)
  call tstsym( develop, thmin, thmax, step, symth )
  !************************************************************************
  !*  Constants                                                           *
  !************************************************************************
  pie   = dacos(-1.d0)
  radfac= pie/180.D0
  rtox  = 2.D0*pie/lambda
  fakt  = lambda*lambda/(2.D0*pie)
  !* nfac is the number of precomputed factors (2n+1)/(n*(n+1))
  nfac  = 0
  !************************************************************************
  !*  distinguish between distribution or not                             *
  !************************************************************************
  if (idis.eq.0) then
     w(1)     = 1.D0
     r(1)     = rmin
     nwithr(1)= 1.D0
     nsub     = 1
     ngauss   = 1
     dr       = 0.D0
  else
     dr = (rmax-rmin)/dble(nsub)
     !*        call gausspt( ngauss, ngauss, (rmax-dr) , rmax ,r ,w )
     call gaulegCP(  ngauss, ngauss, (rmax-dr) , rmax ,r ,w )
     call sizedis( idis, par1, par2, par3, r, ngauss, nwithr, &
          iparts, ndis, rdis, nwrdis )
  endif
  !************************************************************************
  !*  Start integration over radius r with largest radius !               *
  !************************************************************************
  do l=nsub,1,-1
     do i=ngauss,1,-1
        
        sw   = nwithr(i)*w(i)
        x    = rtox*r(i)
        nmax = x + 4.05D0*x**(1.D0/3.D0) + 20
        nfi  = nmax+60
        zabs = x*cdabs(m)
        nD   = zabs + 4.05D0*zabs**(1.D0/3.D0) + 70
        
	allocate(an(max(nD,nfi,nmax)))
        allocate(bn(max(nD,nfi,nmax)))
        allocate(pi(max(nD,nfi,nmax)))
        allocate(tau(max(nD,nfi,nmax)))
        allocate(fi(0:max(nD,nfi,nmax)))
        allocate(chi(0:max(nD,nfi,nmax)))
        allocate(D(max(nD,nfi,nmax)))
        allocate(facf(max(nD,nfi,nmax)))
        allocate(facb(max(nD,nfi,nmax)))
        
        if ((nD.gt.NDn) .or. (nfi.gt.NDn)) then
           write(*,*) ' scatmat: estimated number of Mie-terms:',nD
           write(*,*) '          for particle sizeparameter   :',x
           write(*,*) '          maximum NDn is only          : ',NDn
           stop 'in scatmat no room for all Mie terms'
        endif
        call fichid( m, x, nfi, nmax, nD, fi, chi, D )
        call anbn( m, x, nmax, fi, chi, D, an, bn )
        !************************************************************************
        !*  Precompute the factor (2n+1)/(n*(n+1)) needed in Mie sum over n     *
        !************************************************************************
        if (nmax .gt. nfac) then
           do n=nfac+1, nmax
              facf(n) = dble(2*n+1)/dble(n*(n+1))
              facb(n) = facf(n)
              if (mod(n,2) .eq. 1) facb(n) = -facb(n)
           enddo
           nfac = nmax
        endif
        !************************************************************************
        !*  Calculate extinction and scattering cross section                   *
        !*  Use the convergence criterion to determine the number of terms that *
        !*  will later be used in the mie sum for the scattering matrix itself  *
        !************************************************************************
        Cextsum = 0.D0
        Cscasum = 0.D0
        nstop = nmax
        do  n=1, nmax
           aux = (2.D0*dble(n)+1.D0) * &
                dabs(dble(an(n)*conjg(an(n)) + bn(n)*conjg(bn(n))))
           Cscasum = Cscasum + aux
           Cextsum = Cextsum + (2.D0*n+1.D0)*dble(an(n)+bn(n))
           if (aux .lt. delta) then
              nstop = n
              goto 53
           endif
        enddo
53      nfou = nstop
        if (nfou .ge. nmax) then
           write(*,*) ' WARNING from scatmat : Mie sum not converged for', &
                ' scattering cross section'
           write(*,*) '   radius r = ',r(i),' sizeparameter x = ',x, &
                ' sizedistribution nr. ',idis
           write(*,*) '   Re(m) = ',dble(m),' Im(m) = ',dimag(m)
           write(*,*) '   a priori estimate of number of Mie terms:',nmax
           write(*,*) '   term ',nmax,' for Csca was ',aux
           write(*,*) '   should have been less than ',delta
           write(*,*) '   the a priori estimate will be used'
        endif
        !************************************************************************
        !*  Only for the first run through the loop set points in u= dcos(th)   *
        !************************************************************************
        if ((l.eq.nsub) .and. (i.eq.ngauss)) then
           !************************************************************************
           !*  In case of expansion in GSF : set Gauss points for dcos(th)         *
           !************************************************************************
           if (develop .eq. 1) then
              !************************************************************************
              !*  Ensure accurate integrations: add two terms: nangle = 2*nfou+2      *
              !*  One should be sufficient, but total should be even !                *
              !************************************************************************
              nangle = 2*nfou+2
              if (nangle .gt. NDang) then
                 write(*,*) ' scatmat: need too many integration angles', &
                      ' nangle=',nangle
                 write(*,*) '       maximum array size NDang = ',NDang
                 stop 'too many integration angles'
              endif
              !*             call gausspt(nangle,nangle,-1.d0,1.D0,u,wth)
              call gaulegCP(nangle,nangle,-1.d0,1.D0,u,wth)
           else
              !************************************************************************
              !*  In case no expansion in GSF is desired : set u= dcos(th) for        *
              !*  for equidistant angles between thmin and thmax.                     *
              !************************************************************************
              nangle = idnint((thmax-thmin)/step) + 1
              if (nangle .gt. NDang) then
                 write(*,*) ' scatmat : need too many scattering angles', &
                      ' nangle=',nangle
                 write(*,*) '       maximum array size NDang = ',NDang
                 stop 'too many scattering angles'
              endif
              wfac = 2.d0/dble(nangle)
              do iang=1, nangle
                 th = thmin + dble(iang-1)*step
                 u(nangle+1-iang) = dcos( radfac*th )
                 wth(iang) = wfac
              enddo
           endif
        endif
        !************************************************************************
        !*  Integration for normalization of size distibution, geometrical      *
        !*  cross section and effective radius                                  *
        !************************************************************************
        numpar = numpar+sw
        G      = G     +sw*r(i)*r(i)
        reff   = reff  +sw*r(i)*r(i)*r(i)
        if (symth) then
           !************************************************************************
           !*  Start loop over scattering angles, do only half and use symmetry    *
           !*  between forward and backward scattering angles                      *
           !*  The factor fac90 will later be used to correct for the fact that    *
           !*  for a symmetrical set of scattering angles with an odd number of    *
           !*  angles the scattering matrix is a factor 2 too big at 90 degrees    *
           !*  due to the way we programmed the symmetry relations                 *
           !************************************************************************
           if (mod(nangle,2) .eq. 1) then
              nhalf = (nangle+1)/2
              fac90 = 0.5D0
           else
              nhalf = nangle/2
           endif
           !*
           do j=1, nhalf
              call pitau( u(j), nmax, pi, tau )
              Splusf = dcmplx(0.D0,0.D0)
              Sminf  = dcmplx(0.D0,0.D0)
              Splusb = dcmplx(0.D0,0.D0)
              Sminb  = dcmplx(0.D0,0.D0)
              !*  THIS IS THE INNERMOST LOOP !! (Mie sum)
              !*  can be programmed more efficiently by taking the facf multiplication
              !*  outside the angle loop over index j
              do n=1,nfou
                 Splusf = Splusf + facf(n)*(an(n)+bn(n)) * (pi(n)+tau(n))
                 Sminf  = Sminf  + facf(n)*(an(n)-bn(n)) * (pi(n)-tau(n))
                 Splusb = Splusb + facb(n)*(an(n)+bn(n)) * (pi(n)-tau(n))
                 Sminb  = Sminb  + facb(n)*(an(n)-bn(n)) * (pi(n)+tau(n))
              enddo
              cSplusf = conjg(Splusf)
              cSminf  = conjg(Sminf )
              cSplusb = conjg(Splusb)
              cSminb  = conjg(Sminb )
              k = nangle-j+1
              !*  the forward scattering elements
              F(1,j) = F(1,j) +    sw*(Splusf*cSplusf + Sminf *cSminf)
              F(2,j) = F(2,j) -    sw*(Sminf *cSplusf + Splusf*cSminf)
              F(3,j) = F(3,j) +    sw*(Splusf*cSplusf - Sminf *cSminf)
              F(4,j) = F(4,j) + ci*sw*(Sminf *cSplusf - Splusf*cSminf)
              !*  the backward scattering elements
              F(1,k) = F(1,k) +    sw*(Splusb*cSplusb + Sminb *cSminb)
              F(2,k) = F(2,k) -    sw*(Sminb *cSplusb + Splusb*cSminb)
              F(3,k) = F(3,k) +    sw*(Splusb*cSplusb - Sminb *cSminb)
              F(4,k) = F(4,k) + ci*sw*(Sminb *cSplusb - Splusb*cSminb)
           enddo
        else
           !************************************************************************
           !*  Start loop over scattering angles, do all angles                    *
           !************************************************************************
           do j=1, nangle
              call pitau( u(j), nmax, pi, tau )
              Splusf = dcmplx(0.D0,0.D0)
              Sminf  = dcmplx(0.D0,0.D0)
              !*  THIS IS THE INNERMOST LOOP !! (Mie sum)
              do n=1,nfou
                 Splusf = Splusf + facf(n)*(an(n)+bn(n)) * (pi(n)+tau(n))
                 Sminf  = Sminf  + facf(n)*(an(n)-bn(n)) * (pi(n)-tau(n))
              enddo
              cSplusf = conjg(Splusf)
              cSminf  = conjg(Sminf )
              k = nangle-j+1
              !*  the forward scattering elements
              F(1,j) = F(1,j) +      sw*(Splusf*cSplusf + Sminf *cSminf)
              F(2,j) = F(2,j) -      sw*(Sminf *cSplusf + Splusf*cSminf)
              F(3,j) = F(3,j) +      sw*(Splusf*cSplusf - Sminf *cSminf)
              F(4,j) = F(4,j) + ci*sw*(Sminf *cSplusf - Splusf*cSminf)
           enddo
        endif
        !************************************************************************
        !*  Integration for cross sections, shift radius to next subinterval    *
        !************************************************************************
        Csca = Csca + sw*Cscasum
        Cext = Cext + sw*Cextsum
        r(i) = r(i) - dr
     enddo
     if (l .ne. 1) &
          call sizedis( idis, par1, par2, par3, r, ngauss, nwithr, &
          iparts, ndis, rdis, nwrdis )
  enddo
  !************************************************************************
  !*  End of integration over size distribution                           *
  !*  Some final corrections :                                            *
  !************************************************************************
  do j=1,nangle
     do k=1,4
        F(k,j)= F(k,j)/(2.D0*Csca)
     enddo
  enddo
  if (symth) then
     do k=1,4
        F(k,nhalf) = fac90*F(k,nhalf)
     enddo
  endif
  G     = pie*G
  Csca  = Csca*fakt
  Cext  = Cext*fakt
  Qsca  = Csca/G
  Qext  = Cext/G
  albedo= Csca/Cext
  volume= (4.d0/3.d0)*pie*reff
  reff  = pie*reff/G
  xeff  = rtox*reff
  !*
  miec(1) = Csca
  miec(2) = Cext
  miec(3) = Qsca
  miec(4) = Qext
  miec(5) = albedo
  miec(6) = G
  miec(7) = reff
  miec(8) = xeff
  miec(9) = numpar
  miec(10)= volume
  !*
  
  deallocate(an)
  deallocate(bn)
  deallocate(pi)
  deallocate(tau)
  deallocate(fi)
  deallocate(chi)
  deallocate(D)
  deallocate(facf)
  deallocate(facb)
  
  return
end subroutine scatmat

subroutine sizedis( idis, par1, par2, par3, r, numr, nwithr, &
     iparts, ndis, rdis, nwrdis )
  !************************************************************************
  !*  Calculate the size distribution n(r) for the numr radius values     *
  !*  contained in array r and return the results through the array nwithr*
  !*  The size distributions are normalized such that the integral over   *
  !*  all r is equal to one.                                              *
  !************************************************************************
  implicit double precision (a-h,o-z)
  parameter ( NDdis=300, NDpart = 4 )
  double precision nwithr, logC, logC1, logC2, nwrdis, nwrdisp
  dimension r(numr),nwithr(numr),nwrdis(NDpart,NDdis), &
       nwrdisp(NDdis), rdis(NDpart,NDdis), rdisp(NDdis), &
       y2(NDdis)
  
  pi     = dacos(-1.d0)
  root2p = dsqrt(pi+pi)
  if (idis .eq. 0) return
  write(*,*) ' sizedis: illegal index : ',idis
  stop 'illegal index in sizedis'
  select case(idis)
  case(1)
     !************************************************************************
     !*  1 TWO PARAMETER GAMMA with alpha and b given                        *
     !************************************************************************
     alpha = par1
     b     = par2
     alpha1= alpha+1.D0
     logC  = alpha1*dlog(b)-gammlnCP(alpha1)
     do  i=1, numr
        nwithr(i) = dexp(logC+alpha*dlog(r(i))-b*r(i))
     enddo
  case(2)
     !************************************************************************
     !*  2 TWO PARAMETER GAMMA with par1= reff and par2= veff given          *
     !************************************************************************
     alpha = 1.D0/par2 - 3.D0
     b     = 1.D0/(par1*par2)
     alpha1= alpha+1.D0
     logC  = alpha1*dlog(b)-gammlnCP(alpha1)
     do i=1, numr
        nwithr(i) = dexp(logC+alpha*dlog(r(i))-b*r(i))
     enddo
  case(3)
     !************************************************************************
     !*  3 BIMODAL GAMMA with equal mode weights                             *
     !************************************************************************
     alpha = 1.D0/par3 - 3.D0
     b1    = 1.D0/(par1*par3)
     b2    = 1.D0/(par2*par3)
     gamlna= gammlnCP(alpha+1.D0)
     logC1 = (alpha+1.D0)*dlog(b1)-gamlna
     logC2 = (alpha+1.D0)*dlog(b2)-gamlna
     do i=1, numr
        nwithr(i) = 0.5D0*( dexp(logC1+alpha*dlog(r(i))-b1*r(i)) &
             + dexp(logC2+alpha*dlog(r(i))-b2*r(i)) )
     enddo
  case(4)
     !************************************************************************
     !*  4 LOG NORMAL with rg and sigma given                                *
     !************************************************************************
     flogrg = dlog(par1)
     flogsi = dabs(dlog(par2))
     C      = 1.D0/(root2p*flogsi)
     fac    = -0.5D0/(flogsi*flogsi)
     do i=1, numr
        nwithr(i) = C * dexp( fac*(dlog(r(i))-flogrg)**2 ) / r(i)
     enddo
  case(5)
     !************************************************************************
     !*  5 LOG NORMAL with reff and veff given                               *
     !************************************************************************
     rg     = par1/(1.D0+par2)**2.5D0
     flogrg = dlog(rg)
     flogsi = dsqrt(dlog(1.D0+par2))
     C      = 1.D0/(root2p*flogsi)
     fac    = -0.5D0/(flogsi*flogsi)
     do i=1, numr
        nwithr(i) = C * dexp( fac*(dlog(r(i))-flogrg)**2 ) / r(i)
     enddo
  case(6)
     !************************************************************************
     !*  6 POWER LAW                                                         *
     !************************************************************************
     alpha = par1
     rmin  = par2
     rmax  = par3
     !* DEBUG
     !*     write(7,*) ' sizedis : power law with alpha = ',alpha
     !*     write(7,*) '                          rmin  = ',rmin
     !*     write(7,*) '                          rmax  = ',rmax
     !* END DEBUG
     if (dabs(alpha+1.D0) .lt. 1.d-10) then
        C = 1.D0/dlog(rmax/rmin)
     else
        alpha1 = alpha-1.d0
        C = alpha1 * rmax**alpha1 / ((rmax/rmin)**alpha1-1.d0)
     endif
     do i=1, numr
        if ((r(i) .lt. rmax) .and. (r(i) .gt. rmin)) then
           nwithr(i) = C*r(i)**(-alpha)
        else
           nwithr(i) = 0.D0
        endif
     enddo
  case(7)
     !************************************************************************
     !*  7 MODIFIED GAMMA with alpha, rc and gamma given                     *
     !************************************************************************
     alpha = par1
     rc    = par2
     gamma = par3
     b     = alpha / (gamma*rc**gamma)
     aperg = (alpha+1.D0)/gamma
     logC  = dlog(gamma) + aperg*dlog(b) - gammlnCP(aperg)
     do i=1, numr
        nwithr(i) = dexp( logC + alpha*dlog(r(i)) - b*r(i)**gamma )
     enddo
  case(8)
     !************************************************************************
     !*  8 MODIFIED GAMMA with alpha, b and gamma given                      *
     !************************************************************************
     alpha = par1
     b     = par2
     gamma = par3
     aperg = (alpha+1.D0)/gamma
     logC  = dlog(gamma) + aperg*dlog(b) - gammlnCP(aperg)
     do i=1, numr
        nwithr(i) = dexp( logC + alpha*dlog(r(i)) - b*r(i)**gamma )
     enddo
  case(9)
     !************************************************************************
     !*  9 DISCRETE VALUED DISTRIBUTION                                      *
     !************************************************************************
     do k=1, ndis
        rdisp(k)   = rdis(iparts,k)
        nwrdisp(k) = nwrdis(iparts,k)
     enddo
     call splineCP(rdisp,nwrdisp,ndis,1.d30,1.d30,y2,NDdis,NDdis)
     do j=1,numr
        call splintCP(rdisp,nwrdisp,y2,ndis,r(j),nwithr(j),NDdis,NDdis)
     enddo
  case default
  end select
  !************************************************************************
  !* DEBUG
  !*     write(7,*) ' sizedis:'
  !*     write(7,*) '   i             r(i)               n(r(i))'
  !*     write(7,*) ' ----------------------------------------------------'
  !*     do 1000 i=1, numr
  !*         write(7,1001) i,r(i),nwithr(i)
  !*1000 continue
  !*1001 format(i4,1pe24.14,1pe22.14)
  !* END DEBUG
  return
end subroutine sizedis


function gammlnCP(xarg)
  !************************************************************************
  !*  Return the value of the natural logarithm of the gamma function.    *
  !*  The argument xarg must be real and positive.                        *
  !*  This function is documented in :                                    *
  !*                                                                      *
  !*  W.H. Press et al. 1986, 'Numerical Recipes' Cambridge Univ. Pr.     *
  !*  page 157 (ISBN 0-521-30811)                                         *
  !*                                                                      *
  !*  When the argument xarg is between zero and one, the relation (6.1.4)*
  !*  on page 156 of the book by Press is used.                           *
  !*                                         V.L. Dolman April 18 1989    *
  !************************************************************************
  implicit double precision (a-h,o-z)
  parameter( eps = 1.d-7, one = 1.D0, two = 2.D0, half = 0.5D0, fpf = 5.5D0 )
  dimension cof(6)
  data cof,stp/76.18009173D0,-86.50532033D0, 24.01409822D0, &
       -1.231739516D0, 0.120858003D-2, -0.536382D-5, 2.50662827465D0/
  pi = 4.D0*datan(1.D0)
  if (xarg .le. 0.D0) then
     write(*,*) ' gammlnCP: called with negative argument xarg = ',xarg
     stop 'function gammlnCP called with negative value'
  endif
  if (dabs(xarg-one) .lt. eps) then
     write(*,*) ' gammlnCP: argument too close to one for algorithm'
     stop ' in function gammlnCP argument too close to one'
  endif
  if (xarg .ge. one) then
     xx = xarg
  else
     xx = xarg+two
  endif
  x = xx-one
  tmp = x+fpf
  tmp = (x+half)*dlog(tmp)-tmp
  ser = one
  do j=1, 6
     x = x+one
     ser = ser+cof(j)/x
  enddo
  gtmp = tmp+dlog(stp*ser)
  if  (xarg .gt. one) then
     gammlnCP = gtmp
  else
     pix = pi*(one-xarg)
     gammlnCP = dlog(pix/dsin(pix))-gtmp
  endif
  return
end function gammlnCP


subroutine splineCP(x, y, n, yp1, ypn, y2,NDmui,nn)
  !***********************************************************************
  !** Spline interpolation routine from Press et al. (1986, p.88).      **
  !**                                                                   **
  !** Given arrays x and y of length n containing a tabulated function, **
  !** i.e. y(i)=f(x(i)), with x(1)<x(2)<...<x(n), and given values yp1  **
  !** and ypn for the first derivative of the interpolating function at **
  !** points 1 and n respectively, this routine returns an array y2 of  **
  !** length n which contains the second derivatives of the interpola-  **
  !** ting function at the tabulated points x(i).                       **
  !** If yp1 and/or yp2 are equal to 1x10^30 or larger, the routine is  **
  !** signalled to set the corresponding boundary condition for a natu- **
  !** ral spline, with zero second derivative on that boundary.         **
  !***********************************************************************
  implicit double precision (a-h,o-z)
  parameter (nmax=1500)
  dimension x(NDmui),y(nn),y2(nn),u(nmax)
  
  if (nn .ne. NDmui) then
     write(*,*) ' splineCP : uncomfortable error occurred, '
     write(*,*) ' nn .ne. NDmui '
     stop ' dimension error in splineCP '
  endif
  if (nmax .lt. nn) then
     write(*,*) ' splineCP : nmax may not be smaller than nn '
     write(*,*) ' nmax = ',nmax,' nn = ',nn
     stop ' dimension error in splineCP '
  endif
  if (yp1.gt..99d30) then
     y2(1)=0.d0
     u(1)=0.d0
  else
     y2(1)=-0.5d0
     u(1)=(3.d0/(x(2)-x(1)))*((y(2)-y(1))/(x(2)-x(1))-yp1)
  endif
  do i=2,n-1
     sig=(x(i)-x(i-1))/(x(i+1)-x(i-1))
     p=sig*y2(i-1)+2.d0
     y2(i)=(sig-1.d0)/p
     u(i)=(6.d0*((y(i+1)-y(i))/(x(i+1)-x(i))-(y(i)-y(i-1)) &
          /(x(i)-x(i-1)))/(x(i+1)-x(i-1))-sig*u(i-1))/p
  enddo
  if (ypn.gt..99d30) then
     qn=0.d0
     un=0.d0
  else
     qn=0.5d0
     un=(3.d0/(x(n)-x(n-1)))*(ypn-(y(n)-y(n-1))/(x(n)-x(n-1)))
  endif
  y2(n)=(un-qn*u(n-1))/(qn*y2(n-1)+1.d0)
  do k=n-1,1,-1
     y2(k)=y2(k)*y2(k+1)+u(k)
  enddo
  return
end subroutine splineCP

subroutine splintCP(xa, ya, y2a, n, x, y,NDmui,nn)
  !***********************************************************************
  !** Spline interpolation routine from Press et al. (1986, p.88).      **
  !**                                                                   **
  !** Given the arrays xa and ya of length n, which tabulate a function **
  !** (with the xa(i)'s in order), and given the array y2a, which is    **
  !** the output from splineCP above, and given a value of x, this        **
  !** routine returns a cubic-spline interpolated value y.              **
  !***********************************************************************
  implicit double precision (a-h,o-z)
  
  dimension xa(NDmui),ya(nn),y2a(nn)
  
  if (nn .ne. NDmui) then
     write(*,*) ' a very stupid an unnecessary error occurred '
     write(*,*) ' in splintCP: nn .ne. NDmui '
     stop ' dimension error in splintCP '
  endif
  klo=1
  khi=n
  
1 if (khi-klo.gt.1) then
     k=(klo+khi)/2.d0
     if (xa(k).gt.x) then
        khi=k
     else
        klo=k
     endif
     goto 1
  endif
  h=xa(khi)-xa(klo)
  if (dabs(h).lt.1.d-10) write (7,10)
  a=(xa(khi)-x)/h
  b=(x-xa(klo))/h
  y=a*ya(klo)+b*ya(khi)+ &
       ((a**3-a)*y2a(klo)+(b**3-b)*y2a(khi))*(h**2)/6.d0
  return
10 format('ERROR in routine splintCP: Bad XA input.')
end subroutine splintCP


subroutine tstsym( develop, thmin, thmax, step, symth )
  !************************************************************************
  !*  Test if the set of theta points is symmetrical around 90 degrees    *
  !*  and return this information through logical symth                   *
  !*  In case of development in GSF we have a symmetrical Gauss set !!    *
  !************************************************************************
  implicit double precision(a-h,o-z)
  parameter( eps=1.d-6, heps=0.5d0*eps )
  integer develop
  logical symth
  symth = .false.
  if (develop .eq. 1) then
     symth = .true.
  else
     if ( (dabs( 180.d0 - thmin - thmax ) .lt. eps)  .and. &
          (dmod( thmax-thmin+heps, step ) .lt. eps) ) &
          symth = .true.
  endif
  !* DEBUG
  !c      if (symth) then
  !c          write(7,*) ' tstsym: theta points symmetrical'
  !c      else
  !c          write(7,*) ' tstsym: theta points NOT symmetrical !!'
  !c      endif
  !* END DEBUG
  return
end subroutine tstsym


subroutine DMiLay( rcore, rshell, wvno, rindsh, rindco, mu, &
     numang, qext, qsca, qbs, gqsc, &
     m1, m2, s21, d21, maxang, err)
  
  ! **********************************************************************
  ! DOUBLE PRECISION version of MieLay, which computes electromagnetic
  ! scattering by a stratified sphere, i.e. a particle consisting of a
  ! spherical core surrounded by a spherical shell.  The surrounding
  ! medium is assumed to have refractive index unity.  The formulas,
  ! manipulated to avoid the ill-conditioning that plagued earlier
  ! formulations, were published in:
  
  ! Toon, O. and T. Ackerman, Applied Optics 20, 3657 (1981)
  
  ! Further documentation, including definitons of input and output
  ! arguments, is inside the single precision version of this program
  ! (SUBROUTINE MieLay, available by anonymous ftp from
  ! climate.gsfc.nasa.gov in directory pub/wiscombe).
  
  ! It is recommended to use this DOUBLE PRECISION version for IEEE
  ! arithmetic (32-bit floating point) computers, just to be safe.  If
  ! computer time is critical, back-to-back tests with the single
  ! precision version should be done for ranges of radii and
  ! refractive index relevant to your particular problem, before
  ! adopting the single precision version.  This version is also
  ! recommended for cases of large size parameter (bigger than 10 or
  ! so) and/or large imaginary refractive index (bigger than 1 or so)
  ! and also whenever overflows or strange behavior are encountered in
  ! running the single precision version.  Sometimes the bigger
  ! exponent range in DOUBLE PRECISION is as important as the added
  ! precision.
  
  ! This version is designed to be interchangeable with the single
  ! precision version: all the floating-point input arguments are
  ! still single precision.  Only the name of the routine has been
  ! changed to prevent confusion (and it is strongly urged not to
  ! change it to MieLay for the same reason).
  
  ! Parameters

  integer   mxang, ll, err
  double precision zero, one, two
  parameter ( mxang = 1440, ll = 300000, zero = 0.0d0, one = 1.0d0, &
       two = 2.0d0 )

  ! Scalar Arguments

  integer   maxang, numang
  integer, parameter        :: dp      = selected_real_kind(p=15)
  real (kind=dp)      gqsc, qbs, qext, qsca, rcore, rshell, wvno
  complex (kind=dp)  rindco, rindsh

  ! Array Arguments

  real (kind=dp) mu( numang ), d21( maxang, 2 ), m1( maxang, 2 ), &
       m2( maxang, 2 ), s21( maxang, 2 )

  ! Local Scalars

  logical   inperr, pass1
  integer   j, k, m, n, nmx1, nmx2, nn

  double precision  aa, aim, am1im, am1re, are, bb, bim, bm1im, &
       bm1re, bre, cc, cosx1, cosx4, dd, denom, &
       dgqsc, dqext, dqsca, e2y1, &
       ey1, ey1my4, ey1py4, ey4, fourpi, pinum, &
       rmm, rx, sinx1, sinx4, toler, x1, x4, &
       xcore, xshell, y1, y4 
  
  double complex  ac, acoe, acoem1, bc, bcoe, bcoem1, ci, czero, &
       dh1, dh2, dh4, dummy, dumsq, k1, k2, k3, &
       p24h21, p24h24, rrfx, sback, wm1

  ! Local Arrays

  double precision  pi( mxang, 3 ), si2tht( mxang ), t( 5 ), &
       ta( 4 ), tau( mxang, 3 )

  double complex  s1( mxang, 2 ), s2( mxang, 2 ), &
       u( 8 ), wfn( 2 ), z( 4 )
  double complex,allocatable :: w(:,:),acap(:)
  
  ! External Functions

  logical   wrtbad, wrtdim
  external  wrtbad, wrtdim

  ! External Subroutines

  external  errmsg

  ! Intrinsic Functions

  intrinsic abs, dimag, asin, dcmplx, cos, exp, mod, dble, sin

  ! Save statement

  save  pinum, pass1

  ! Data statements

  data pass1 / .true. / , toler / 1.d-6 / , &
       czero / ( 0.d0, 0.d0 ) / , ci / ( 0.d0, 1.d0 ) /

  allocate(w(3,ll))
  allocate(acap(ll))

  if( pass1 ) then

     pinum  = two*asin( one )
     pass1  = .false.

  endif

  xshell = rshell*wvno
  xcore  = rcore*wvno
  t( 1 ) = xshell*abs( rindsh )
  nmx1   = 1.1d0*t( 1 )
  nmx2   = t( 1 )
  
  if( nmx1.le.150 ) then

     nmx1   = 150
     nmx2   = 135

  endif

  ! ** Check input arguments for gross errors
  inperr = .False.

  if( wvno.le.0.0 ) inperr = wrtbad( 'wvno' )

  if( rshell.le.0.0 ) inperr = wrtbad( 'rshell' )

  if( rcore.le.0.0 .or. rcore.gt.rshell ) &
       inperr = wrtbad( 'rcore' )

  if( real(rindsh).le.0.0 .or. aimag(rindsh).gt.0.0 ) &
       inperr = wrtbad( 'rindsh' )

  if( real(rindco).le.0.0 .or. aimag(rindco).gt.0.0 ) &
       inperr = wrtbad( 'rindco' )

  if( numang.lt.0 ) inperr = wrtbad( 'numang' )

  if( numang.gt.mxang ) inperr = wrtdim( 'mxang', numang )

  if( numang.gt.maxang ) inperr = wrtdim( 'maxang', numang )

  if( nmx1 + 1 .gt. ll ) inperr = wrtdim( 'll', nmx1 + 1 )

  do j = 1, numang
     if( mu(j).lt.- toler .or. mu(j).gt. 1.0+toler ) &
          inperr = wrtbad( 'mu' )
  enddo

  if( inperr ) then
     err=1
     return
     call errmsg( &
          'MIELAY--Input argument errors.  Aborting...', .True. )
  endif
  
  k1     = rindco*wvno
  k2     = rindsh*wvno
  k3     = dcmplx( wvno )
  z( 1 ) = rindsh*xshell
  z( 2 ) = xshell
  z( 3 ) = rindco*xcore
  z( 4 ) = rindsh*xcore
  x1     = dble( z(1) )
  y1     = dimag( z(1) )
  x4     = dble( z(4) )
  y4     = dimag( z(4) )
  rx     = one / xshell
  
  ! ** Down-recurrence for A function
  acap( nmx1 + 1 ) = czero
  do m = 1, 3
     w( m, nmx1 + 1 ) = czero
  enddo
  
  rrfx  = one / ( rindsh*xshell)
  do nn = nmx1, 1, - 1
     acap( nn ) = ( ( nn + 1)*rrfx ) - &
          one / ( ( (nn + 1)*rrfx) + acap( nn + 1) )
     do m = 1, 3
        w( m, nn ) = ( ( nn + 1) / z( m + 1) ) - &
             one / ( ( (nn + 1)/z(m + 1)) + w( m, nn + 1) )
     enddo
  enddo
  
  do j = 1, numang
     si2tht( j ) = one - mu( j )**2
     pi( j, 1 )  = zero
     pi( j, 2 )  = one
     tau( j, 1 ) = zero
     tau( j, 2 ) = mu( j )
  enddo
  
  ! ** Initialization of homogeneous sphere
  
  t( 1 )   = cos( xshell )
  t( 2 )   = sin( xshell )
  wm1      = dcmplx( t(1), - t(2) )
  wfn( 1 ) = dcmplx( t(2), t(1) )
  ta( 1 )  = t( 2 )
  ta( 2 )  = t( 1 )
  wfn( 2 ) = rx*wfn( 1 ) - wm1
  ta( 3 )  =  dble( wfn(2) )
  ta( 4 )  = dimag( wfn(2) )
  
  ! ** Initialization procedure for stratified sphere
  n      = 1
  sinx1  = sin( x1 )
  sinx4  = sin( x4 )
  cosx1  = cos( x1 )
  cosx4  = cos( x4 )
  ey1    = exp( y1 )
  e2y1   = ey1**2
  ey4    = exp( y4 )
  ey1my4 = exp( y1 - y4 )
  ey1py4 = ey1*ey4
  aa     = sinx4*( ey1py4 + ey1my4 )
  bb     = cosx4*( ey1py4 - ey1my4 )
  cc     = sinx1*( e2y1 + one )
  dd     = cosx1*( e2y1 - one )
  denom  = one + e2y1*( 4.0d0*sinx1**2 - two + e2y1 )
  dummy  = dcmplx( ( aa*cc + bb*dd) / denom, &
       ( bb*cc - aa*dd) / denom )
  dummy  = dummy*( acap(n) + n / z(1) ) / ( w(3, n) + n / z(4) )
  dumsq  = dummy**2
  
  p24h24 = 0.5d0 + dcmplx( sinx4**2 - 0.5d0, cosx4*sinx4 )*ey4**2
  p24h21 = 0.5d0*dcmplx( sinx1*sinx4 - cosx1*cosx4, &
       sinx1*cosx4 + cosx1*sinx4 )*ey1py4 &
       + 0.5d0*dcmplx( sinx1*sinx4 + cosx1*cosx4, &
       - sinx1*cosx4 + cosx1*sinx4 )*ey1my4
  dh1    = z( 1 ) / ( one + ci*z( 1) ) - one / z( 1 )
  dh2    = z( 2 ) / ( one + ci*z( 2) ) - one / z( 2 )
  dh4    = z( 4 ) / ( one + ci*z( 4) ) - one / z( 4 )
  p24h24 = p24h24 / ( ( dh4 + n/z(4))*( w(3, n) + n/z(4)) )
  p24h21 = p24h21 / ( ( dh1 + n/z(1))*( w(3, n) + n/z(4)) )
  
  u( 1 ) = k3*acap( n ) - k2*w( 1, n )
  u( 2 ) = k3*acap( n ) - k2*dh2
  u( 3 ) = k2*acap( n ) - k3*w( 1, n )
  u( 4 ) = k2*acap( n ) - k3*dh2
  u( 5 ) = k1*w( 3, n ) - k2*w( 2, n )
  u( 6 ) = k2*w( 3, n ) - k1*w( 2, n )
  u( 7 ) = - ci*( dummy*p24h21 - p24h24 )
  u( 8 ) = ta( 3 ) / wfn( 2 )
  
  acoe  = u( 8 )*( u(1)*u(5)*u(7) + k1*u(1) - dumsq*k3*u(5) ) / &
       ( u(2)*u(5)*u(7) + k1*u(2) - dumsq*k3*u(5) )
  
  bcoe  = u( 8 )*( u(3)*u(6)*u(7) + k2*u(3) - dumsq*k2*u(6) ) / &
       ( u(4)*u(6)*u(7) + k2*u(4) - dumsq*k2*u(6) )
  
  acoem1 = acoe
  bcoem1 = bcoe
  are    =  dble( acoe )
  aim    = dimag( acoe )
  bre    =  dble( bcoe )
  bim    = dimag( bcoe )

  dqext  = 3.d0*( are + bre )
  dqsca  = 3.d0*( are**2 + aim**2 + bre**2 + bim**2 )
  dgqsc  = zero
  sback  = 3.d0*( acoe - bcoe )
  rmm    = one
  
  ac  = 1.5d0*acoe
  bc  = 1.5d0*bcoe
  do j = 1, numang
     s1( j, 1 ) = ac*pi( j, 2 ) + bc*tau( j, 2 )
     s1( j, 2 ) = ac*pi( j, 2 ) - bc*tau( j, 2 )
     s2( j, 1 ) = bc*pi( j, 2 ) + ac*tau( j, 2 )
     s2( j, 2 ) = bc*pi( j, 2 ) - ac*tau( j, 2 )
  enddo
  
  ! ***************** Start of Mie summing loop ******************

  n  = 2
70 continue
  !                        ** Recurrences for functions little-pi,
  !                           little-tau of Mie theory
  t( 1 ) = 2*n - 1
  t( 2 ) = n - 1
  do  j = 1, numang
     pi( j, 3 ) = ( t( 1)*pi( j, 2)*mu( j) - n*pi( j, 1) ) / t( 2 )
     
     tau( j, 3 ) = mu( j )*( pi( j, 3) - pi( j, 1) ) - &
          t( 1 )*si2tht( j )*pi( j, 2 ) + tau( j, 1 )
  enddo
  
  !                               ** Here set up homogeneous sphere
  wm1      = wfn( 1 )
  wfn( 1 ) = wfn( 2 )
  wfn( 2 ) = t( 1 )*rx*wfn( 1 ) - wm1
  ta( 1 )  =  dble( wfn( 1) )
  ta( 2 )  = dimag( wfn( 1) )
  ta( 3 )  =  dble( wfn( 2) )
  ta( 4 )  = dimag( wfn( 2) )

  !                               ** Here set up stratified sphere

  dh1    = - n / z( 1 ) + one / ( n / z( 1) - dh1 )
  dh2    = - n / z( 2 ) + one / ( n / z( 2) - dh2 )
  dh4    = - n / z( 4 ) + one / ( n / z( 4) - dh4 )
  p24h24 = p24h24 / ( ( dh4 + n/z(4))*( w(3, n) + n/z(4)) )
  p24h21 = p24h21 / ( ( dh1 + n/z(1))*( w(3, n) + n/z(4)) )
  dummy  = dummy*( acap(n) + n / z(1) ) / ( w(3, n) + n / z(4) )
  dumsq  = dummy**2

  u( 1 ) = k3*acap( n ) - k2*w( 1, n )
  u( 2 ) = k3*acap( n ) - k2*dh2
  u( 3 ) = k2*acap( n ) - k3*w( 1, n )
  u( 4 ) = k2*acap( n ) - k3*dh2
  u( 5 ) = k1*w( 3, n ) - k2*w( 2, n )
  u( 6 ) = k2*w( 3, n ) - k1*w( 2, n )
  u( 7 ) = - ci*( dummy*p24h21 - p24h24 )
  u( 8 ) = ta( 3 ) / wfn( 2 )
  
  acoe  = u( 8 )*( u(1)*u(5)*u(7) + K1*U(1) - DUMSQ*K3*U(5) ) / &
       ( U(2)*U(5)*U(7) + K1*U(2) - DUMSQ*K3*U(5) )
  
  bcoe  = u( 8 )*( u(3)*u(6)*u(7) + k2*u(3) - dumsq*k2*u(6) ) / &
       ( u(4)*u(6)*u(7) + k2*u(4) - dumsq*k2*u(6) )
  are  =  dble( acoe )
  aim  = dimag( acoe )
  bre  =  dble( bcoe )
  bim  = dimag( bcoe )

  !                        ** Increment sums for efficiency factors

  am1re  =  dble( acoem1 )
  am1im  = dimag( acoem1 )
  bm1re  =  dble( bcoem1 )
  bm1im  = dimag( bcoem1 )
  t( 4 ) = (2*n - one) / ( n*(n - one) )
  t( 2 ) = (n - one)*(n + one) / n
  dgqsc  = dgqsc + t( 2 )*( am1re*are + am1im*aim + &
       bm1re*bre + bm1im*bim ) + &
       t( 4 )*( am1re*bm1re + am1im*bm1im )
  
  t( 3 )  = 2*n + 1
  dqext   = dqext + t( 3 )*( are + bre )
  t( 4 )  = are**2 + aim**2 + bre**2 + bim**2
  dqsca   = dqsca + t( 3 )*t( 4 )
  rmm     = - rmm
  sback  = sback + t( 3 ) * rmm *( acoe - bcoe )
  
  t( 2 ) = n*( n + 1 )
  t( 1 ) = t( 3 ) / t( 2 )

  ac  = t( 1 )*acoe
  bc  = t( 1 )*bcoe
  do j = 1, numang
     s1( j, 1 ) = s1( j, 1 ) + ac*pi( j, 3 ) + bc*tau( j, 3 )
     s2( j, 1 ) = s2( j, 1 ) + bc*pi( j, 3 ) + ac*tau( j, 3 )
  enddo

  !                            ** Scattering matrix elements for
  !                               supplements of 0-90 degree scattering
  !                               angles submitted by user
  if( mod(n, 2).eq.0 ) then

     do j = 1, numang
        s1( j, 2 ) = s1( j, 2 ) - ac*pi( j, 3 ) + bc*tau( j, 3 )
        s2( j, 2 ) = s2( j, 2 ) - bc*pi( j, 3 ) + ac*tau( j, 3 )
     enddo
     
  else

     do j = 1, numang
        s1( j, 2 ) = s1( j, 2 ) + ac*pi( j, 3 ) - bc*tau( j, 3 )
        s2( j, 2 ) = s2( j, 2 ) + bc*pi( j, 3 ) - ac*tau( j, 3 )
     enddo
     
  endif

  !                                   ** Test for convergence of sums
  if( t(4).ge.1.0d-14 ) then

     n  = n + 1

     if( n.gt.nmx2 ) then
        err=1
        return
        call errmsg( &
             'MIELAY--Dimensions for W,ACAP not enough. Suggest'// &
             ' get detailed output, modify routine', .True. )
     endif

     do j = 1, numang
        
        pi( j, 1 ) = pi( j, 2 )
        pi( j, 2 ) = pi( j, 3 )
        tau( j, 1 ) = tau( j, 2 )
        tau( j, 2 ) = tau( j, 3 )
        
     enddo

     acoem1 = acoe
     bcoem1 = bcoe

     go to 70

  endif

  !  ***************** End of summing loop ******************

  !                            ** Transform complex scattering amplitudes
  !                               into elements of real scattering matrix

  do j = 1, numang
     
     do k = 1, 2
        
        m1( j, k )  = dble( s1(j, k) )**2 + dimag( s1(j, k) )**2
        m2( j, k )  = dble( s2(j, k) )**2 + dimag( s2(j, k) )**2
        s21( j, k ) = dble(  s1(j, k) )*dble(  s2(j, k) ) + &
             dimag( s1(j, k) )*dimag( s2(j, k) )
        d21( j, k ) = dimag( s1(j, k) )*dble( s2(j, k) ) - &
             dimag( s2(j, k) )*dble( s1(j, k) )
        
     enddo
     
  enddo

  t( 1 ) = two*rx**2
  qext   = t( 1 )*dqext
  qsca   = t( 1 )*dqsca
  gqsc   = two*t( 1 )*dgqsc
  sback  = 0.5*sback
  qbs    = ( dble(sback)**2 + dimag(sback)**2 ) / (pinum*xshell**2)
  
  return
end subroutine DMiLay

subroutine errmsg( messag, fatal )

  ! Print out a warning or error message;  abort if error
  ! after making symbolic dump (machine-specific)
  
  ! Provenance:  the 3 error-handlers ErrMsg, WrtBad, WrtDim are
  !               borrowed from MIEV, the Wiscombe Mie program
  
  ! Scalar Arguments

  character messag*( * )
  logical   fatal

  ! Local Scalars

  logical   msglim
  integer   maxmsg, nummsg

  ! External Subroutines

  ! cccc EXTERNAL  SYMDUMP

  ! Save statement

  save      maxmsg, nummsg, msglim

  ! Data statements

  data      nummsg / 0 / , maxmsg / 100 / , msglim / .false. /

  if( fatal ) then

     write( *, '(//,2A,//)' ) ' ****** ERROR *****  ', messag

     !        ** Example symbolic dump call for Cray
     ! cccc    CALL SYMDUMP( '-B -c3' )

     write(*,*) 'I should actually stop, but... whatever'!STOP

  endif

  nummsg = nummsg + 1

  if( msglim ) return

  if( nummsg.le.maxmsg ) then
     
     write( *, '(/,2A,/)' ) ' ****** WARNING *****  ', messag

  else

     write( *, '(//,A,//)' ) &
          ' ****** TOO MANY WARNING MESSAGES --  ' // &
          'They will no longer be printed *******'
     
     msglim = .True.

  endif

end subroutine errmsg

logical function WrtBad( VARNAM )

  ! Write names of erroneous variables and return 'TRUE'
  
  ! INPUT :   VarNam = Name of erroneous variable to be written
  !                    ( CHARACTER, any length )
  
  ! Scalar Arguments

  character varnam*( * )
  
  ! Local Scalars
  
  integer   maxmsg, nummsg
  
  ! External Subroutines
  
  external  errmsg
  
  ! Save statement
  
  save      nummsg, maxmsg

  ! Data statements
  
  data      nummsg / 0 / , maxmsg / 50 /
  
  wrtbad = .true.
  nummsg = nummsg + 1
  write( *, '(3A)' ) ' ****  Input variable  ', varnam, &
       '  in error  ****'
  
  if( nummsg.eq.maxmsg ) call errmsg( &
      'Too many input errors.  Aborting...', .true. )
  
end function WrtBad
   

logical function WrtDim( dimnam, minval )

  ! Write name of too-small symbolic dimension and
  ! the value it should be increased to;  return 'TRUE'

  ! INPUT :  DimNam = Name of symbolic dimension which is too small
  !                   ( CHARACTER, any length )
  !          Minval = Value to which that dimension should be
  !                   increased (at least)
  
  ! Scalar Arguments

  character dimnam*( * )
  integer   minval

  write( *, '(3A,I7)' ) ' ****  Symbolic dimension  ', &
       dimnam, '  should be increased to at least ', minval
  write (*,*) 'Take a look at subroutine dmilay()'
  wrtdim = .true.
  
end function WrtDim

!!! **** File Variables

! Local Variables:
! eval: (outline-minor-mode)
! outline-regexp: "!!!\\|\\(program\\|subroutine\\|function\\|module\\)\\>"
! outline-heading-alist: (("!!!" . 1) ("program" . 2) ("subroutine" . 2) ("module" . 2) ("function" . 2))
! End:
