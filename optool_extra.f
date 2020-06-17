	subroutine MeerhoffMie(rad,lam,e1,e2,Csca,Cext,F11,F12,F33,
     >	F34,na)
	IMPLICIT NONE
	integer na,nangle
	real*8 rad,lam,e1,e2,csca,cext,F11(na),F12(na),F33(na),F34(na)

	integer nparts,develop,nsubr(1),ngaur(1),idis(1),ndis,i
	real*8 delta,cutoff,thmin,thmax,step,wavel(1),Rem(1),fImm(1)
	real*8 par1(1),par2(1),par3(1),rdis(1,300),nwrdis(1,300)
	real*8 F(4,6000),theta(6000)

	nangle=na
	nparts=1
	develop=0
	delta=1d-8
	cutoff=1d-8
	thmin=180d0*(real(1)-0.5)/real(na)
	thmax=180d0*(real(na)-0.5)/real(na)
	step=(thmax-thmin)/real(na-1)
	wavel=lam
	Rem=e1
	fImm=e2
	nsubr=1
	ngaur=1
	idis=0
	par1=rad
	par2=0d0
	par3=0d0
	ndis=1
	rdis(1,1)=rad
	nwrdis(1,1)=1d0

      call mie(nparts, develop, delta, cutoff, thmin, thmax, step
     +               , wavel, Rem, fImm, nsubr , ngaur, idis
     +               , par1, par2, par3, ndis, rdis, nwrdis
     +               , nangle, theta, F ,Cext, Csca)

	do i=1,nangle
		F11(i)=F(1,nangle-i+1)
		F12(i)=F(2,nangle-i+1)
		F33(i)=F(3,nangle-i+1)
		F34(i)=F(4,nangle-i+1)
	enddo

	return
	end




      subroutine mie(nparts, develop, delta, cutoff, thmin, thmax, step
     +               , wavel, Rem, fImm, nsubr , ngaur, idis
     +               , par1, par2, par3, ndis, rdis, nwrdis
     +               , nangle, theta, F ,outCext, outCsca)


      implicit double precision (a-h,o-z)
      implicit integer (i-n)
      parameter ( NDang=6000, NDcoef=6000, NDpart = 4, nrunit=88 )
      parameter ( NDdis=300 )
      double precision lambda, nr, ni, miec, nwrdis,outCext,outCsca
      integer develop
      double complex m
      character*60 s1, s2, s3, s4, s5, s6
      character*10 scfile
      dimension F(4,NDang)
     +        , miec(13),u(NDang),w(NDang),coefs(4,4,0:NDcoef)
     +        , wavel(*), Rem(*), fImm(*), nsubr(*)
     +        , ngaur(*), idis(*)
     +        , par1(*), par2(*), par3(*)
     +        , rdis(1,*), nwrdis(1,*)
     +        , scfile(NDpart), theta(NDang)
      pi     = dacos(-1.D0)
      radtod = 180.D0/pi
      scfile(1) = 'mie1_sc  '
      scfile(2) = 'mie2_sc  '
      scfile(3) = 'mie3_sc  '
      scfile(4) = 'mie4_sc  '
************************************************************************
*  Start a loop over different 'particles' referring to the different  *
*  cases specified in the input                                        *
************************************************************************
      do 50 iparts=1,nparts
c      write(7,920) iparts
************************************************************************
*  In case of development in generalized spherical functions a file    *
*  is opened on which the expansion coefficients will be written       *
************************************************************************
c      if (develop .ne. 0)
c     +    call opensc(scfile(iparts),nrunit,'new')
************************************************************************
*  Determine the integration bounds for the radius integration         *
************************************************************************
      call rminmax( idis(iparts), nsubr(iparts), ngaur(iparts)
     +   , par1(iparts), par2(iparts), par3(iparts), cutoff
     +   , rmin, rmax
     +   , iparts, ndis, rdis, nwrdis )
*
      lambda = wavel(iparts)
      nr     = Rem(iparts)
      ni     = fImm(iparts)
      m      = dcmplx(nr,-ni)
************************************************************************
*  Calculate the scattering matrix : this is most of the work !        *
************************************************************************
      call scatmat( u, w, m, wavel(iparts), idis(iparts)
     +            , thmin, thmax, step, develop
     +            , nsubr(iparts), ngaur(iparts), rmin, rmax
     +            , par1(iparts), par2(iparts), par3(iparts)
     +            , delta, F, miec, nangle
     +            , ndis, rdis, nwrdis, iparts )
************************************************************************
*  Test if the number of coefficients needed fits in the array coefs   *
*  The number of coefficients is equal to the number of integration    *
*  points in scattering angle !! (think about this !)                  *
************************************************************************
      ncoef = nangle
      if ((ncoef .gt. NDcoef) .and. (develop .eq. 1)) then
         write(*,*) ' main: too many coefficients needed :',ncoef
         write(*,*) '       maximum array size NDcoef = ',NDcoef
         write(*,*) '       Therefore I cannot do the expansion.'
         develop=0
         call clossc(nrunit)
      endif
************************************************************************
*  If required, calculate the expansion coefficients                   *
************************************************************************
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
************************************************************************
*  Print to standard output the following data :                       *
*                                                                      *
*  miec(1) = Csca : the average scattering cross section               *
*  miec(2) = Cext : the average extinction cross section               *
*  miec(3) = Qsca : the scattering efficiency factor                   *
*  miec(4) = Qext : the extinction efficiency factor                   *
*  miec(5) = a    : the single scattering albedo                       *
*  miec(6) = G    : the average geometrical cross section              *
*  miec(7) = reff : the average effective radius                       *
*  miec(8) = xeff : the average effective size parameter               *
*  miec(9) = num  : the integrated number of particles                 *
*  miec(10)= V    : the average volume                                 *
************************************************************************
c      write(7,922) miec(1)
c      write(7,923) miec(2)
c      write(7,924) miec(3)
c      write(7,925) miec(4)
c      write(7,926) miec(5)
c      write(7,927) miec(6)
c      write(7,928) miec(7)
c      write(7,929) miec(8)
c      write(7,930) miec(9)
c      write(7,931) miec(10)

	outCsca=miec(1)
	outCext=miec(2)

c      if (develop .eq. 0) goto 1020

************************************************************************
*  In case a development in GSF is required, prepare the heading of    *
*  coefficient file.                                                   *
************************************************************************
* ruler: '.........1.........2.........3.........4.........5.........6'
c      s1='PRODUCED BY THE MEERHOFF MIE PROGRAM VERSION 3.1            '
c      write(s2,1001) wavel(iparts), nr, ni
c 1001 format('lambda= ',f11.7,' Re(m) = ',f11.7,' Im(m) = ',f11.7)
************************************************************************
c      if (idis(iparts) .ne. 0) then
c          write(s5,1007) rmin, rmax
c          write(s6,1009) nsubr(iparts), ngaur(iparts)
c      endif
c 1007 format('rmin  = ',f11.7,' rmax  = ',f11.7,'                     ')
c 1009 format('r integration with ',i3,' interval(s) of '
c     +                                            ,i3,' Gauss points  ')
************************************************************************
c      if (idis(iparts) .eq. 0) then
c      s5='                                                            '
c      s6='                                                            '
c      s3='No size distribution (single particle)                      '
c      write(s4,1003) par1(iparts)
c 1003 format('Radius = ',1f11.7
c     +                       ,'                                       ')
************************************************************************
c      else if (idis(iparts) .eq. 1) then
c      s3='Size distribution : two parameter gamma                     '
c      write(s4,1011) par1(iparts), par2(iparts)
c 1011 format('alpha = ',f11.7,' b    = ',f11.7,'                    ')
************************************************************************
c      else if (idis(iparts) .eq. 2) then
c      s3='Size distribution : two parameter gamma                     '
c      write(s4,1012) par1(iparts), par2(iparts)
c 1012 format('reff  = ',f11.7,' veff  = ',f11.7,'                   ')
************************************************************************
c      else if (idis(iparts) .eq. 3) then
c      s3='Size distribution : bimodal gamma with equal mode weights   '
c      write(s4,1013) par1(iparts), par2(iparts), par3(iparts)
c 1013 format('reff1 = ',f11.7,' reff2 = ',f11.7,' veff = ',f11.7)
************************************************************************
c      else if (idis(iparts) .eq. 4) then
c      s3='Size distribution : log normal                              '
c      write(s4,1014) par1(iparts), par2(iparts)
c 1014 format('rg     = ',f11.7,' sigma = ',f11.7,'                    ')
************************************************************************
c      else if (idis(iparts) .eq. 5) then
c      s3='Size distribution : log normal                              '
c      write(s4,1015) par1(iparts), par2(iparts)
c 1015 format('reff   = ',f11.7,' veff  = ',f11.7,'                    ')
************************************************************************
c      else if (idis(iparts) .eq. 6) then
c      s3='Size distribution : power law                               '
c      write(s4,1016) par1(iparts), par2(iparts), par3(iparts)
c 1016 format('alpha = ',f11.7,' rmin  = ',f11.7,' rmax  = ',f11.7)
************************************************************************
c      else if (idis(iparts) .eq. 7) then
c      s3='Size distribution : modified gamma                          '
c      write(s4,1017) par1(iparts), par2(iparts), par3(iparts)
c 1017 format('alpha = ',f11.7,' rc    = ',f11.7,' gamma = ',f11.7)
************************************************************************
c      else if (idis(iparts) .eq. 8) then
c      s3='Size distribution : modified gamma                          '
c      write(s4,1018) par1(iparts), par2(iparts), par3(iparts)
c 1018 format('alpha = ',f11.7,' b     = ',f11.7,' gamma = ',f11.7)
************************************************************************
c      else if (idis(iparts) .eq. 9) then
c      s3='Size distribution : discrete values                         '
c      write(s4,1019) par1(iparts), par2(iparts), par3(iparts)
c 1019 format('min radius = ',f11.7,' max = ',f11.7,
c     +       ' ngaus = ',f5.1)
************************************************************************
c      else
c          write(7,*) ' MAIN: illegal size distribution nr:',idis(iparts)
c      endif
      if (nangle .ge. 1) then
          cosbar = coefs(1,1,1)/3.D0
      else
          cosbar = 0.D0
      endif
************************************************************************
*  Write the coefficients out onto file                                *
************************************************************************
c      call writsc(nrunit,coefs,NDcoef,ncoef,miec(5),cosbar
c     +        ,miec(4),miec(3),miec(6),miec(10),s1, s2, s3, s4, s5, s6 )
************************************************************************
*  Print a table of the scattering matrix and the linear polarization  *
************************************************************************
c 1020 write (7,911)
      do 30 j=1,nangle
          i     = nangle-j+1
          xp    = -100.d0*F(2,i)/F(1,i)
          theta1 = radtod*dacos(u(i))
c          write(7,913) j,theta1,F(1,i),F(2,i),F(3,i),F(4,i),xp
   30 continue
************************************************************************
*  Print the coefficients on output.                                   *
*  Evaluate the expansion in GSF at an equidistant set of scattering   *
*  angles and print the resulting scattering matrix.                   *
************************************************************************
      if (develop.eq.1) then
c          write (7,914)
          do 40 i=0,ncoef
c              write(7,915) i,coefs(1,1,i),coefs(2,2,i),coefs(3,3,i),
c     +                       coefs(4,4,i),coefs(1,2,i),coefs(3,4,i)
   40    continue
         call expand( ncoef, npunt, coefs, u, F, thmin, thmax, step )
c         write (7,911)
         do 45 i=1,npunt
             xp    = -100.d0*F(2,i)/F(1,i)
             theta1 = radtod*dacos(u(i))
c             write(7,913) i,theta1,F(1,i),F(2,i),F(3,i),F(4,i),xp
   45    continue
      endif
   50 continue
************************************************************************
*  End of loop over 'particles'.                                       *
************************************************************************
c      if (develop .ne. 0) call clossc(nrunit)

c  800 format(1x )
c  801 format(51x,i2)
c  802 format(51x,e8.2)
c  803 format(51x,i1)
c  804 format(51x,f5.2)
c  805 format(8x,3(2x,f7.4),4x,i1,4x,i2,4x,i4)
c  806 format(8x,5(2x,f7.2))

c  911 format(/,'    num   scat.angle     F11              F21   ',
c     +           '           F33             F34         pol',
c     +'arization',/,t90,'[proc]')
c  913 format(4x,i3,4x,f6.2,4(3x,e14.8),3x,f7.3)
c  914 format(1h1,13x,'alpha1',11x,'alpha2',11x,'alpha3',11x,
c     +'alpha4',12x,'beta1',12x,'beta2',/)
c  915 format(1h ,'l= ',i3,6f17.12)
c  920 format('1PARTICLE NUMBER ',i3,/)
c  921 format(1h ,'the number of computed scattering angles was :'
c     +,i7)
c  922 format(1h ,'Csca : the average scattering cross section  :'
c     +,e24.14)
c  923 format(1h ,'Cext : the average extinction cross section  :'
c     +,e24.14)
c  924 format(1h ,'Qsca : the scattering efficiency factor      :'
c     +,e24.14)
c  925 format(1h ,'Qext : the extinction efficiency factor      :'
c     +,e24.14)
c  926 format(1h ,'a    : the single scattering albedo          :'
c     +,e24.14)
c  927 format(1h ,'G    : the average geometrical cross section :'
c     +,e24.14)
c  928 format(1h ,'reff : the average effective radius          :'
c     +,e24.14)
c  929 format(1h ,'xeff : the average effective size parameter  :'
c     +,e24.14)
c  930 format(1h ,'num  : the integrated number of particles    :'
c     +,e24.14)
c  931 format(1h ,'V    : the average volume                    :'
c     +,e24.14)

c      stop 'program mie terminated'
	return
      end

      subroutine anbn( m, x, nmax, psi, chi, d, an, bn )
************************************************************************
*  Calculate the Mie coefficients an and bn.                           *
************************************************************************
      implicit double precision (a-h,o-z)
      parameter( NDn=10000000 )
      double complex m, zn, znm1, save, perm
      double complex an(nmax), bn(nmax), D(nmax)
      dimension psi(0:nmax), chi(0:nmax)
      perm = 1.D0/m
      perx = 1.D0/x
      xn   = 0.D0
* DEBUG
*     write(7,*) ' anbn:'
*     write(7,*) '     Re(an)           Im(an)       Re(bn)      Im(bn)'
* END DEBUG
      do 100 n=1, nmax
          zn   = dcmplx( psi(n),   chi(n))
          znm1 = dcmplx( psi(n-1), chi(n-1))
          xn   = dble(n)*perx
          save = D(n)*perm+xn
          an(n)= (save*psi(n)-psi(n-1)) / (save*zn-znm1)
          save = m*D(n)+xn
          bn(n)= (save*psi(n)-psi(n-1)) / (save*zn-znm1)
* DEBUG
*          write(72,*) n,an(n),bn(n),psi(n),chi(n),D(n)
* END DEBUG
  100 continue
      return
      end
      subroutine clossc(iunit)
************************************************************************
*  On entry :                                                          *
*      iunit     number of the unit to be closed                       *
*  On exit :                                                           *
*      The file is closed.                                             *
*                                                 VLD January 9, 1989  *
************************************************************************
c      close(unit=iunit,err=999)
      return
c  999 write(7,*) ' clossc: error in closing file with unit number',iunit
c      stop 'in clossc error in closing file'
      end
      subroutine devel(ncoef,nangle,u,w,F,coefs)
************************************************************************
*  Calculate the expansion coefficients of the scattering matrix in    *
*  generalized spherical functions by numerical integration over the   *
*  scattering angle.                                                   *
************************************************************************
      implicit double precision (a-h,o-z)
      parameter( NDcoef=3000, NDang = 6000 )
      dimension u(NDang),w(NDang),F(4,NDang)
      dimension coefs(4,4,0:NDcoef), P00(NDang,2), P02(NDang,2)
     +                             , P22(NDang,2), P2m2(NDang,2)
************************************************************************
*  Initialization                                                      *
************************************************************************
      qroot6 = -0.25D0*dsqrt(6.D0)
	do 22 j=0,NDcoef
	   do 21 i=1,4
              do 20 ii=1,4
		 coefs(ii,i,j)=0.D0
 20	      continue
 21	   continue
 22	continue
 30	continue
 40	continue
************************************************************************
*  Multiply the scattering matrix F with the weights w for all angles  *
*  We do this here because otherwise it should be done for each l      *
************************************************************************
      do 60 k=1,4
          do 50 i=1, nangle
              F(k,i) = w(i)*F(k,i)
   50     continue
   60 continue
************************************************************************
*  Start loop over the coefficient index l                             *
*  first update generalized spherical functions, then calculate coefs. *
*  lold and lnew are pointer-like indices used in recurrence           *
************************************************************************
      lnew = 1
      lold = 2

      do 70 l=0, ncoef
               if (l .eq. 0) then
************************************************************************
*             Adding paper Eq. (77) with m=0                           *
************************************************************************
              do 80 i=1, nangle
                  P00(i,lold) = 1.D0
                  P00(i,lnew) = 0.D0
                  P02(i,lold) = 0.D0
                  P22(i,lold) = 0.D0
                  P2m2(i,lold)= 0.D0
                  P02(i,lnew) = 0.D0
                  P22(i,lnew) = 0.D0
                  P2m2(i,lnew)= 0.D0
   80         continue
          else
              fac1 = (2.D0*l-1.d0)/dble(l)
              fac2 = dble(l-1.d0)/dble(l)
************************************************************************
*             Adding paper Eq. (81) with m=0                           *
************************************************************************
              do 90 i=1, nangle
                  P00(i,lold) = fac1*u(i)*P00(i,lnew) - fac2*P00(i,lold)
   90         continue
          endif
          if (l .eq. 2) then
************************************************************************
*             Adding paper Eqs. (78) and (80)  with m=2                *
*             sql4 contains the factor dsqrt(l*l-4) needed in          *
*             the recurrence Eqs. (81) and (82)                        *
************************************************************************
              do 100 i=1, nangle
                  P02(i,lold) = qroot6*(1.D0-u(i)*u(i))
                  P22(i,lold) = 0.25D0*(1.D0+u(i))*(1.D0+u(i))
                  P2m2(i,lold)= 0.25D0*(1.D0-u(i))*(1.D0-u(i))
                  P02(i,lnew) = 0.D0
                  P22(i,lnew) = 0.D0
                  P2m2(i,lnew)= 0.D0
  100         continue
              sql41 = 0.D0
          else if (l .gt. 2) then
************************************************************************
*             Adding paper Eq. (82) with m=0 and m=2                   *
************************************************************************
              sql4  = sql41
              sql41 = dsqrt(dble(l*l)-4.d0)
              twol1 = 2.D0*dble(l)-1.d0
              tmp1  = twol1/sql41
              tmp2  = sql4/sql41
              denom = (dble(l)-1.d0)*(dble(l*l)-4.d0)
              fac1  = twol1*(dble(l)-1.d0)*dble(l)/denom
              fac2  = 4.D0*twol1/denom
              fac3  = dble(l)*((dble(l)-1.d0)*(dble(l)-1.d0)-4.d0)/denom
              do 110 i=1, nangle
                  P02(i,lold) = tmp1*u(i)*P02(i,lnew) - tmp2*P02(i,lold)
                  P22(i,lold) = (fac1*u(i)-fac2)*P22(i,lnew)
     +                                                - fac3*P22(i,lold)
                  P2m2(i,lold)= (fac1*u(i)+fac2)*P2m2(i,lnew)
     +                                               - fac3*P2m2(i,lold)
  110         continue
          endif
************************************************************************
*         Switch indices so that lnew indicates the function with      *
*         the present index value l, this mechanism prevents swapping  *
*         of entire arrays.                                            *
************************************************************************
          itmp = lnew
          lnew = lold
          lold = itmp
************************************************************************
*         Now calculate the coefficients by integration over angle     *
*         See de Haan et al. (1987) Eqs. (68)-(73).                    *
*         Remember for Mie scattering : F11 = F22 and F33 = F44        *
************************************************************************
          alfap = 0.D0
          alfam = 0.D0
          do 120 i=1, nangle
              coefs(1,1,l) = coefs(1,1,l) + P00(i,lnew)*F(1,i)
              alfap = alfap + P22(i,lnew)*(F(1,i)+F(3,i))
              alfam = alfam + P2m2(i,lnew)*(F(1,i)-F(3,i))
              coefs(4,4,l) = coefs(4,4,l) + P00(i,lnew)*F(3,i)
              coefs(1,2,l) = coefs(1,2,l) + P02(i,lnew)*F(2,i)
              coefs(3,4,l) = coefs(3,4,l) + P02(i,lnew)*F(4,i)
  120     continue
************************************************************************
*         Multiply with trivial factors like 0.5D0*(2*l+1)             *
************************************************************************
          fl = dble(l)+0.5D0
          coefs(1,1,l) =  fl*coefs(1,1,l)
          coefs(2,2,l) =  fl*0.5D0*(alfap+alfam)
          coefs(3,3,l) =  fl*0.5D0*(alfap-alfam)
          coefs(4,4,l) =  fl*coefs(4,4,l)
          coefs(1,2,l) =  fl*coefs(1,2,l)
          coefs(3,4,l) =  fl*coefs(3,4,l)
          coefs(2,1,l) =     coefs(1,2,l)
          coefs(4,3,l) =    -coefs(3,4,l)
   70 continue
************************************************************************
*     End of loop over index l                                         *
************************************************************************

************************************************************************
*     Remove the weight factor from the scattering matrix              *
************************************************************************
      do 140 k=1, 4
          do 130 i=1, nangle
              F(k,i) = F(k,i)/w(i)
  130     continue
  140 continue
      return
      end
      subroutine expand(ncoef,nangle,coefs,u,F,thmin,thmax,step)
c
c     ***********************************************************
c
      implicit double precision (a-h,o-z)
      parameter( pi=3.141592653589793238462643D0, radfac= pi/180.D0 )
      parameter( NDcoef=3000, NDang=6000 )
      dimension u(NDang),F(4,NDang)
      dimension coefs(4,4,0:NDcoef), P00(NDang,2), P02(NDang,2)
c
c     **********************************************************
c
c                           initialize

      do 10 i=1, nangle
          u(i) = dcos(radfac*(thmin+dble(i-1)*step))
   10 continue
      qroot6 = -0.25D0*dsqrt(6.D0)
************************************************************************
*  Set scattering matrix F to zero                                     *
************************************************************************
      do 30 k=1,4
          do 20 i=1, nangle
              F(k,i) = 0.D0
   20     continue
   30 continue
************************************************************************
*  Start loop over the coefficient index l                             *
*  first update generalized spherical functions, then calculate coefs. *
*  lold and lnew are pointer-like indices used in recurrence           *
************************************************************************
      lnew = 1
      lold = 2
      do 90 l=0, ncoef
               if (l .eq. 0) then
************************************************************************
*             Adding paper Eqs. (76) and (77) with m=0                 *
************************************************************************
              do 40 i=1, nangle
                  P00(i,lold) = 1.D0
                  P00(i,lnew) = 0.D0
                  P02(i,lold) = 0.D0
                  P02(i,lnew) = 0.D0
   40         continue
          else
              fac1 = (2.D0*dble(l)-1.d0)/dble(l)
              fac2 = (dble(l)-1.d0)/dble(l)
************************************************************************
*             Adding paper Eq. (81) with m=0                           *
************************************************************************
              do 50 i=1, nangle
                  P00(i,lold) = fac1*u(i)*P00(i,lnew) - fac2*P00(i,lold)
   50         continue
          endif
          if (l .eq. 2) then
************************************************************************
*             Adding paper Eq. (78)                                    *
*             sql4 contains the factor dsqrt((l+1)*(l+1)-4) needed in  *
*             the recurrence Eqs. (81) and (82)                        *
************************************************************************
              do 60 i=1, nangle
                  P02(i,lold) = qroot6*(1.D0-u(i)*u(i))
                  P02(i,lnew) = 0.D0
   60         continue
              sql41 = 0.D0
          else if (l .gt. 2) then
************************************************************************
*             Adding paper Eq. (82) with m=0                           *
************************************************************************
              sql4  = sql41
              sql41 = dsqrt(dble(l*l)-4.d0)
              tmp1  = (2.D0*dble(l)-1.d0)/sql41
              tmp2  = sql4/sql41
              do 70 i=1, nangle
                  P02(i,lold) = tmp1*u(i)*P02(i,lnew) - tmp2*P02(i,lold)
   70         continue
          endif
************************************************************************
*         Switch indices so that lnew indicates the function with      *
*         the present index value l, this mechanism prevents swapping  *
*         of entire arrays.                                            *
************************************************************************
          itmp = lnew
          lnew = lold
          lold = itmp
************************************************************************
*         Now add the l-th term to the scattering matrix.              *
*         See de Haan et al. (1987) Eqs. (68)-(73).                    *
*         Remember for Mie scattering : F11 = F22 and F33 = F44        *
************************************************************************
          do 80 i=1, nangle
              F(1,i) = F(1,i) + coefs(1,1,l)*P00(i,lnew)
              F(2,i) = F(2,i) + coefs(1,2,l)*P02(i,lnew)
              F(3,i) = F(3,i) + coefs(4,4,l)*P00(i,lnew)
              F(4,i) = F(4,i) + coefs(3,4,l)*P02(i,lnew)
   80     continue
   90 continue
************************************************************************
*     End of loop over index l                                         *
************************************************************************
      return
      end
      subroutine fichid( m, x, nchi, nmax, nD, psi, chi, D )
************************************************************************
*  Calculate functions psi(x)  chi(x) and D(z) where z = mx.           *
*  On entry, the following should be supplied :                        *
*      m      : complex index of refraction                            *
*      x      : sizeparameter                                          *
*      nchi   : starting index for backwards recurrence of chi         *
*      nmax   : number of chi, psi and D that must be available        *
*      nD     : starting index for backwards recurrence of D           *
*  On exit, the desired functions are returned through psi, chi and D  *
************************************************************************
      implicit double precision (a-h,o-z)
      double complex D,m,z, perz, zn1
      dimension psi(0:nchi), chi(0:nmax+1), D(nd)
*
      z = m*x
      perz = 1.D0/z
      perx = 1.D0/x
      sinx = dsin(x)
      cosx = dcos(x)
      psi(nchi)=0d0
************************************************************************
*  (mis-) use the array psi to calculate the functions rn(x)
*  De Rooij and van der Stap Eq. (A6)
************************************************************************
      do 10 n=nchi-1, 0, -1
          psi(n) = 1.D0 / (dble(2*n+1)/x - psi(n+1))
   10 continue
************************************************************************
*  Calculate functions D(z) by backward recurrence
*  De Rooij and van der Stap Eq. (A11)
************************************************************************
      D(nD) = dcmplx(0.D0,0.D0)
      do 20 n=nD - 1, 1, -1
          zn1 = dble(n+1)*perz
          D(n) = zn1 - 1.D0/(D(n+1)+zn1)
   20 continue
************************************************************************
*  De Rooij and van der Stap Eqs. (A3) and (A1)
*  De Rooij and van der Stap Eq. (A5) where psi(n) was equal to r(n)
*  and Eq. (A2)
************************************************************************
      psi(0) = sinx
      psi1   = psi(0)*perx - cosx
      if (dabs(psi1) .gt. 1.d-4) then
          psi(1) = psi1
          do 30 n=2,nmax
              psi(n) = psi(n)*psi(n-1)
   30     continue
      else
          do 35 n=1,nmax
              psi(n) = psi(n)*psi(n-1)
   35     continue
      endif
************************************************************************
*  De Rooij and van der Stap Eqs. (A4) and (A2)
************************************************************************
      chi(0) = cosx
      chi(1) = chi(0)*perx + sinx
      do 40 n=1, nmax-1
          chi(n+1) = dble(2*n+1)*chi(n)*perx - chi(n-1)
   40 continue
* DEBUG
*     write(7,*) ' fichid: x = ',x
*     write(7,12)
*     do 26 n=0, nchi
*         write(7,11) n,psi(n),chi(n),D(n)
*  26 continue
*  11 format(i4,4e24.14)
*  12 format('   n',t20,'psi(n)',t44,'chi(n)',t68
*    +      ,'Re(D(n))',t92,'Im(D(n))',/
*    +      ,' ----------------------------------------------------'
*    +      ,'-----------------------------------------------------')
* END DEBUG
      return
      end
      subroutine gaulegCP(ndim,ngauss,a,b,x,w)
************************************************************************
*   Given the lower and upper limits of integration a and b, and given *
*   the number of Gauss-Legendre points ngauss, this routine returns   *
*   through array x the abscissas and through array w the weights of   *
*   the Gauss-Legendre quadrature formula. Eps is the desired accuracy *
*   of the abscissas. This routine is documented further in :          *
*   W.H. Press et al. 'Numerical Recipes' Cambridge Univ. Pr. (1987)   *
*   page 125 ISBN 0-521-30811-9                                        *
************************************************************************
      implicit double precision (a-h,o-z)
      parameter( eps = 1.d-14 )
      double precision x(ndim), w(ndim)
      double precision a,b,xm,xl,z,p1,p2,p3,pp,z1,pi
      pi=4.D0*datan(1.D0)
      m=(ngauss+1)/2
      xm=0.5D0*(a+b)
      xl=0.5D0*(b-a)
      do 12 i=1,m
*         THIS IS A REALLY CLEVER ESTIMATE :
          z= dcos(pi*(dble(i)-0.25D0)/(dble(ngauss)+0.5D0))
    1     continue
              p1=1.D0
              p2=0.D0
              do 11 j=1,ngauss
                  p3= p2
                  p2= p1
                  p1=((dble(2*j)-1.d0)*z*p2-(dble(j)-1.d0)*p3)/dble(j)
   11         continue
              pp=ngauss*(z*p1-p2)/(z*z-1.d0)
              z1= z
              z= z1-p1/pp
          if (dabs(z-z1) .gt. eps) goto 1
          x(i)= xm-xl*z
          x(ngauss+1-i)= xm+xl*z
          w(i)=2.D0*xl/((1.D0-z*z)*pp*pp)
******          write(7,*) ' gaulegCP: ',i,'  x=',x(i),' w=',w(i)
          w(ngauss+1-i)= w(i)
   12 continue
      return
      end
      subroutine gausspt(ndim,ngauss,a,b,x,w)
c
c     ***********************************************************
c
c     put the gauss-legendre points of order ndim in array x,
c     the weights in array w. the points are transformed to the
c     interval [a,b]
c
c     ***********************************************************
c
      implicit double precision (a-h,o-z)
      dimension x(ndim),w(ndim)
c
c     ***********************************************************
c
c                     find starting values
c
      gn=0.5D0/dble(ngauss)
      extra=1.0D0/(.4D0*dble(ngauss)*dble(ngauss)+5.0D0)
      xz=-gn
      nt=0
      nteken=0
    5 pnm2=1.0D0
      pnm1= xz
      do 10 i=2,ngauss
          pnm1xz= pnm1*xz
          pn=2.0D0*pnm1xz-pnm2-(pnm1xz-pnm2)/dble(i)
          pnm2= pnm1
          pnm1= pn
   10 continue
      mteken=1
      if (pn .le. 0.0D0) mteken=-1
      if ((mteken+nteken) .eq. 0) then
          nt=nt+1
          x(nt)= xz
      endif
      nteken=mteken
      if ((1.0D0-xz) .le. extra) go to 30
      xz= xz+(1.D0-xz*xz)*gn+extra
      go to 5
   30 continue
c
c     ***********************************************************
c
c                determine zero's and weights
c
      do 60 i=1,nt
          xz= x(i)
          delta2=1.D0
   35         pnm2=1.0D0
              pnm1= xz
              pnm1af=1.0D0
              z=.5D0+1.5D0*xz*xz
              do 40 k=2,ngauss
                  pnm1xz= pnm1*xz
                  pn=2.0D0*pnm1xz-pnm2-(pnm1xz-pnm2)/dble(k)
                  pnaf= xz*pnm1af+k*pnm1
                  z= z+(dble(k)+0.5D0)*pn*pn
                  pnm2= pnm1
                  pnm1= pn
                  pnm1af= pnaf
   40         continue
              delta1= pn/pnaf
              xz= xz-delta1
              if(delta1.lt.0.0D0) delta1=-delta1
              if((delta1.ge.delta2).and.(delta2.lt.1.d-14)) go to 50
              delta2= delta1
          go to 35
   50     x(i)= xz
          w(i)=1.0D0/z
   60 continue
c
c     ***********************************************************
c
c                 transform to the interval [a,b]
c
      nghalf=ngauss/2
      ngp1=ngauss+1
      ntp1=nt+1
      apb= a+b
      bmag2=(b-a)/2.0D0
	do 90 i=1,nghalf
          x(ngp1-i)= b-bmag2*(1.0D0-x(ntp1-i))
	  w(ngp1-i)= bmag2*w(ntp1-i)
 90	continue 
      if (nghalf .ne. nt) then
          x(nt)= apb/2.0D0
          w(nt)= w(1)*bmag2
      endif
      do 120 i=1,nghalf
          x(i)= apb-x(ngp1-i)
          w(i)= w(ngp1-i)
  120 continue
      return
      end


      subroutine pitau(u,nmax,pi,tau)
c     ***********************************************************
c     calculates pi,n(u) and tau,n(u) with upward recursion
c
c     ***********************************************************
      implicit double precision (a-h,o-z)
      dimension pi(nmax),tau(nmax)
c
c       starting values:
      pi(1) = 1.D0
      pi(2) = 3.D0*u
      delta = 3.D0*u*u-1.d0
      tau(1)= u
      tau(2)= 2.D0*delta-1.d0

c       upward recursion:
      do 10 n=2, nmax-1
          pi(n+1) = dble(n+1)/dble(n)*delta + u*pi(n)
          delta   = u*pi(n+1) - pi(n)
          tau(n+1)= dble(n+1)*delta - pi(n)
   10 continue
      return
      end
      subroutine rminmax( idis, nsubr, ngaur, par1, par2, par3, cutoff
     +                  , rmin, rmax
     +                        , iparts, ndis, rdis, nwrdis )
************************************************************************
*  Find the integration bounds rmin and rmax for the integration over  *
*  a size distribution. These bounds are chosen such that the size     *
*  distribution falls below the user specified cutoff. It is essential *
*  that the size distribution is normalized such that the integral     *
*  over all r is equal to one !                                        *
*  This is programmed rather clumsy and will in the future be changed  *
************************************************************************
      implicit double precision (a-h,o-z)
      parameter( eps = 1.d-10, NDpart = 4, NDdis=300 )
      double precision nwithr, rdis(NDpart,NDdis), nwrdis(NDpart,NDdis)
      dimension r(1), nwithr(1)
*
      if (idis.eq.0) then
         rmin= par1
         rmax= par1
      else
          goto (10,20,30,40,50,60,70,80,90) idis
          write(*,*) ' rminmax: illegal size distribution index :',idis
          stop 'in rminmax illegal size distribution index'
************************************************************************
   10     sef = 1.D0/dsqrt(par2+3.D0)
          ref = 1.D0/(sef*sef*par2)
          rref= ref
          goto 100

   20     ref = par1
          sef = dsqrt(par2)
          rref= ref
          goto 100

   30     sef = dsqrt(par3)
          ref = dmax1(par1,par2)+sef
          rref= dmin1(par1,par2)
          goto 100

   40     sef = dsqrt(dexp(dlog(par2)**2)-1.d0)
          ref = par1*(1.D0+sef*sef)**0.4D0
          rref= ref
          goto 100

   50     ref = par1
          sef = dsqrt(ref)
          rref= ref
          goto 100

   60     rmin= par2
          rmax= par3
          goto 999

   70     ref = par2
          sef = 2.D0*ref
          rref=0.5D0*ref
          goto 100

   80     ref = (par1/(par2*par3))**par3
          sef = 2.D0*ref
          rref= 0.5D0*ref
          goto 100

   90     rmin = par1
          rmax = par2
          goto 999

************************************************************************
*  search for a value of r such that the size distribution
*  is less than the cutoff. Start the search at ref+sef which          *
*  guarantees that such a value will be found on the TAIL of the       *
*  distribution.                                                       *
************************************************************************
  100     r(1) = ref+sef
          r0   = ref
  200          call sizedis( idis, par1, par2, par3, r, 1, nwithr
     +                                , iparts, ndis, rdis, nwrdis )
          if (nwithr(1) .gt. cutoff) then
              r0   = r(1)
              r(1) = 2.D0*r(1)
              goto 200
          endif
          r1 = r(1)
************************************************************************
*  Now the size distribution assumes the cutoff value somewhere        *
*  between r0 and r1  Use bisection to find the corresponding r        *
************************************************************************
  300     r(1) = 0.5D0*(r0+r1)
          call sizedis( idis, par1, par2, par3, r, 1, nwithr
     +                                , iparts, ndis, rdis, nwrdis )
          if (nwithr(1) .gt. cutoff) then
              r0 = r(1)
          else
              r1 = r(1)
          endif
          if ((r1-r0) .gt. eps) goto 300
          rmax = 0.5D0*(r0+r1)
************************************************************************
*  Search for a value of r on the low end of the size distribution     *
*  such that the distribution falls below the cutoff. There is no      *
*  guarantee that such a value exists, so use an extra test to see if  *
*  the search comes very near to r = 0                                 *
************************************************************************
          r1 = rref
          r0 = 0.D0
  400     r(1) = 0.5D0*r1
          call sizedis( idis, par1, par2, par3, r, 1, nwithr
     +                                , iparts, ndis, rdis, nwrdis )
          if (nwithr(1) .gt. cutoff) then
              r1 = r(1)
              if (r1 .gt. eps) goto 400
          else
              r0 = r(1)
          endif
************************************************************************
*  Possibly the size distribution goes through cutoff between r0       *
*  and r1 try to find the exact value of r where this happens by       *
*  bisection.                                                          *
*  In case there is no solution, the algorithm will terminate soon.    *
************************************************************************
  500     r(1) = 0.5D0*(r0+r1)
          call sizedis( idis, par1, par2, par3, r, 1, nwithr
     +                                , iparts, ndis, rdis, nwrdis )
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
      end

      subroutine scatmat( u, wth, m, lambda, idis
     +                  , thmin, thmax, step, develop
     +                  , nsub, ngauss, rmin, rmax
     +                  , par1, par2, par3
     +                  , delta, F, miec, nangle
     +                  , ndis, rdis, nwrdis, iparts )
************************************************************************
*  Calculate the scattering matrix of an ensemble of homogenous        *
*  spheres. On entry, the following must be supplied :                 *
*     m            : complex index of refraction                       *
*     lambda       : wavelength                                        *
*     idis         : index of the size distribution                    *
*     nsub         : number of subintervals for integration over r     *
*     ngauss       : number of Gauss points used per subinterval       *
*     rmin         : lower bound for integration over r                *
*     rmax         : upper bound for integration over r                *
*     par1,2,3     : parameters of the size distribution               *
*     delta        : cutoff used in truncation of the Mie sum          *
*     thmin        : minimum scattering angle in degrees               *
*     thmax        : maximum scattering angle in degrees               *
*     step         : step in scattering angle in degrees               *
*     develop      : expansion in GSF (1) or not (0)                   *
*  On exit, the following results are returned :                       *
*     u            : cosines of scattering angles                      *
*     wth          : Gaussian weights associated with u                *
*     F            : scattering matrix for all cosines in u            *
*     miec         : array containing cross sections etc.              *
*     nangle       : the number of scattering angles                   *
*  When develop=1 u contains Gauss points, when develop=0 u contains   *
*  cosines of equidistant angles between thmin and thmax with step.    *
************************************************************************
      implicit double precision (a-h,o-z)
      parameter( NDn=10000000, NDr=1000, NDang=6000, NDdis=300,NDpart=4)
      double complex   m, ci, Splusf, Sminf, cSplusf
      double complex   cSminf, Splusb, Sminb, cSplusb, cSminb

      double complex, allocatable :: an(:), bn(:),D(:)
      double precision,allocatable :: pi(:), tau(:), fi(:), chi(:)
      double precision,allocatable :: facf(:), facb(:)

      double precision lambda, nwithr, miec, numpar, thmin, thmax, step
     +               , nwrdis
      integer     develop
      dimension   u(NDang), wth(NDang), F(4,NDang),rdis(NDpart, NDdis)
     	dimension   miec(13), nwrdis(NDpart, NDdis),r(NDr), w(NDr),
     +	nwithr(NDr)
      logical     symth
************************************************************************
*  Initialization                                                      *
************************************************************************
      do 11 j=1,NDang
          do 10 k=1,4
	     F(k,j)=0.D0
 10	  continue
 11	continue
      Csca  = 0.D0
      Cext  = 0.D0
      numpar= 0.D0
      G     = 0.D0
      reff  = 0.D0
      nfou  = 0
      fac90 = 1.D0
      ci    = dcmplx(0.D0,1.D0)
      call tstsym( develop, thmin, thmax, step, symth )
************************************************************************
*  Constants                                                           *
************************************************************************
      pie   = dacos(-1.d0)
      radfac= pie/180.D0
      rtox  = 2.D0*pie/lambda
      fakt  = lambda*lambda/(2.D0*pie)
* nfac is the number of precomputed factors (2n+1)/(n*(n+1))
      nfac  = 0
************************************************************************
*  distinguish between distribution or not                             *
************************************************************************
      if (idis.eq.0) then
         w(1)     = 1.D0
         r(1)     = rmin
         nwithr(1)= 1.D0
         nsub     = 1
         ngauss   = 1
         dr       = 0.D0
      else
         dr = (rmax-rmin)/dble(nsub)
*        call gausspt( ngauss, ngauss, (rmax-dr) , rmax ,r ,w )
         call gaulegCP(  ngauss, ngauss, (rmax-dr) , rmax ,r ,w )
         call sizedis( idis, par1, par2, par3, r, ngauss, nwithr
     +               , iparts, ndis, rdis, nwrdis )
      endif
************************************************************************
*  Start integration over radius r with largest radius !               *
************************************************************************
      do 60 l=nsub,1,-1
      do 50 i=ngauss,1,-1
*
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
************************************************************************
*  Precompute the factor (2n+1)/(n*(n+1)) needed in Mie sum over n     *
************************************************************************
      if (nmax .gt. nfac) then
          do 26 n=nfac+1, nmax
              facf(n) = dble(2*n+1)/dble(n*(n+1))
              facb(n) = facf(n)
              if (mod(n,2) .eq. 1) facb(n) = -facb(n)
   26     continue
          nfac = nmax
      endif
************************************************************************
*  Calculate extinction and scattering cross section                   *
*  Use the convergence criterion to determine the number of terms that *
*  will later be used in the mie sum for the scattering matrix itself  *
************************************************************************
      Cextsum = 0.D0
      Cscasum = 0.D0
      nstop = nmax
      do 52 n=1, nmax
          aux = (2.D0*dble(n)+1.D0) *
     +             dabs(dble(an(n)*conjg(an(n)) + bn(n)*conjg(bn(n))))
          Cscasum = Cscasum + aux
          Cextsum = Cextsum + (2.D0*n+1.D0)*dble(an(n)+bn(n))
          if (aux .lt. delta) then
              nstop = n
              goto 53
          endif
   52 continue
   53 nfou = nstop
      if (nfou .ge. nmax) then
          write(*,*) ' WARNING from scatmat : Mie sum not converged for'
     +           ,' scattering cross section'
          write(*,*) '   radius r = ',r(i),' sizeparameter x = ',x
     +           ,' sizedistribution nr. ',idis
          write(*,*) '   Re(m) = ',dble(m),' Im(m) = ',dimag(m)
          write(*,*) '   a priori estimate of number of Mie terms:',nmax
          write(*,*) '   term ',nmax,' for Csca was ',aux
          write(*,*) '   should have been less than ',delta
          write(*,*) '   the a priori estimate will be used'
      endif
************************************************************************
*  Only for the first run through the loop set points in u= dcos(th)   *
************************************************************************
      if ((l.eq.nsub) .and. (i.eq.ngauss)) then
************************************************************************
*  In case of expansion in GSF : set Gauss points for dcos(th)         *
************************************************************************
          if (develop .eq. 1) then
************************************************************************
*  Ensure accurate integrations: add two terms: nangle = 2*nfou+2      *
*  One should be sufficient, but total should be even !                *
************************************************************************
              nangle = 2*nfou+2
              if (nangle .gt. NDang) then
                 write(*,*) ' scatmat: need too many integration angles'
     +                  ,' nangle=',nangle
                 write(*,*) '       maximum array size NDang = ',NDang
                 stop 'too many integration angles'
              endif
*             call gausspt(nangle,nangle,-1.d0,1.D0,u,wth)
              call gaulegCP(nangle,nangle,-1.d0,1.D0,u,wth)
          else
************************************************************************
*  In case no expansion in GSF is desired : set u= dcos(th) for        *
*  for equidistant angles between thmin and thmax.                     *
************************************************************************
              nangle = idnint((thmax-thmin)/step) + 1
              if (nangle .gt. NDang) then
                 write(*,*) ' scatmat : need too many scattering angles'
     +                  ,' nangle=',nangle
                 write(*,*) '       maximum array size NDang = ',NDang
                 stop 'too many scattering angles'
              endif
              wfac = 2.d0/dble(nangle)
              do 260 iang=1, nangle
                  th = thmin + dble(iang-1)*step
                  u(nangle+1-iang) = dcos( radfac*th )
                  wth(iang) = wfac
  260         continue
          endif
      endif
************************************************************************
*  Integration for normalization of size distibution, geometrical      *
*  cross section and effective radius                                  *
************************************************************************
      numpar = numpar+sw
      G      = G     +sw*r(i)*r(i)
      reff   = reff  +sw*r(i)*r(i)*r(i)
      if (symth) then
************************************************************************
*  Start loop over scattering angles, do only half and use symmetry    *
*  between forward and backward scattering angles                      *
*  The factor fac90 will later be used to correct for the fact that    *
*  for a symmetrical set of scattering angles with an odd number of    *
*  angles the scattering matrix is a factor 2 too big at 90 degrees    *
*  due to the way we programmed the symmetry relations                 *
************************************************************************
        if (mod(nangle,2) .eq. 1) then
            nhalf = (nangle+1)/2
            fac90 = 0.5D0
        else
            nhalf = nangle/2
        endif
*
        do 40 j=1, nhalf
            call pitau( u(j), nmax, pi, tau )
            Splusf = dcmplx(0.D0,0.D0)
            Sminf  = dcmplx(0.D0,0.D0)
            Splusb = dcmplx(0.D0,0.D0)
            Sminb  = dcmplx(0.D0,0.D0)
*  THIS IS THE INNERMOST LOOP !! (Mie sum)
*  can be programmed more efficiently by taking the facf multiplication
*  outside the angle loop over index j
            do 20 n=1,nfou
                Splusf = Splusf + facf(n)*(an(n)+bn(n)) * (pi(n)+tau(n))
                Sminf  = Sminf  + facf(n)*(an(n)-bn(n)) * (pi(n)-tau(n))
                Splusb = Splusb + facb(n)*(an(n)+bn(n)) * (pi(n)-tau(n))
                Sminb  = Sminb  + facb(n)*(an(n)-bn(n)) * (pi(n)+tau(n))
   20       continue
            cSplusf = conjg(Splusf)
            cSminf  = conjg(Sminf )
            cSplusb = conjg(Splusb)
            cSminb  = conjg(Sminb )
            k = nangle-j+1
*  the forward scattering elements
            F(1,j) = F(1,j) +    sw*(Splusf*cSplusf + Sminf *cSminf)
            F(2,j) = F(2,j) -    sw*(Sminf *cSplusf + Splusf*cSminf)
            F(3,j) = F(3,j) +    sw*(Splusf*cSplusf - Sminf *cSminf)
            F(4,j) = F(4,j) + ci*sw*(Sminf *cSplusf - Splusf*cSminf)
*  the backward scattering elements
            F(1,k) = F(1,k) +    sw*(Splusb*cSplusb + Sminb *cSminb)
            F(2,k) = F(2,k) -    sw*(Sminb *cSplusb + Splusb*cSminb)
            F(3,k) = F(3,k) +    sw*(Splusb*cSplusb - Sminb *cSminb)
            F(4,k) = F(4,k) + ci*sw*(Sminb *cSplusb - Splusb*cSminb)
   40   continue
      else
************************************************************************
*  Start loop over scattering angles, do all angles                    *
************************************************************************
        do 400 j=1, nangle
            call pitau( u(j), nmax, pi, tau )
            Splusf = dcmplx(0.D0,0.D0)
            Sminf  = dcmplx(0.D0,0.D0)
*  THIS IS THE INNERMOST LOOP !! (Mie sum)
            do 200 n=1,nfou
                Splusf = Splusf + facf(n)*(an(n)+bn(n)) * (pi(n)+tau(n))
                Sminf  = Sminf  + facf(n)*(an(n)-bn(n)) * (pi(n)-tau(n))
  200       continue
            cSplusf = conjg(Splusf)
            cSminf  = conjg(Sminf )
            k = nangle-j+1
*  the forward scattering elements
            F(1,j) = F(1,j) +      sw*(Splusf*cSplusf + Sminf *cSminf)
            F(2,j) = F(2,j) -      sw*(Sminf *cSplusf + Splusf*cSminf)
            F(3,j) = F(3,j) +      sw*(Splusf*cSplusf - Sminf *cSminf)
            F(4,j) = F(4,j) + ci*sw*(Sminf *cSplusf - Splusf*cSminf)
  400   continue
      endif
************************************************************************
*  Integration for cross sections, shift radius to next subinterval    *
************************************************************************
      Csca = Csca + sw*Cscasum
      Cext = Cext + sw*Cextsum
      r(i) = r(i) - dr
   50 continue
      if (l .ne. 1)
     +         call sizedis( idis, par1, par2, par3, r, ngauss, nwithr
     +                     , iparts, ndis, rdis, nwrdis )
   60 continue
************************************************************************
*  End of integration over size distribution                           *
*  Some final corrections :                                            *
************************************************************************
      do 71 j=1,nangle
          do 70 k=1,4
	     F(k,j)= F(k,j)/(2.D0*Csca)
 70	  continue
 71	continue
	if (symth) then
	   do 80 k=1,4
	      F(k,nhalf) = fac90*F(k,nhalf)
 80	   continue
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
*
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
*

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
      end
      subroutine sizedis( idis, par1, par2, par3, r, numr, nwithr
     +                  , iparts, ndis, rdis, nwrdis )
************************************************************************
*  Calculate the size distribution n(r) for the numr radius values     *
*  contained in array r and return the results through the array nwithr*
*  The size distributions are normalized such that the integral over   *
*  all r is equal to one.                                              *
************************************************************************
      implicit double precision (a-h,o-z)
      parameter ( NDdis=300, NDpart = 4 )
      double precision nwithr, logC, logC1, logC2, nwrdis, nwrdisp
      dimension r(numr),nwithr(numr),nwrdis(NDpart,NDdis)
     +        , nwrdisp(NDdis), rdis(NDpart,NDdis), rdisp(NDdis)
     +        , y2(NDdis)
*
      pi     = dacos(-1.d0)
      root2p = dsqrt(pi+pi)
      if (idis .eq. 0) return
      goto (10, 20, 30, 40, 50, 60, 70, 80, 90 ) idis
      write(*,*) ' sizedis: illegal index : ',idis
      stop 'illegal index in sizedis'
************************************************************************
*  1 TWO PARAMETER GAMMA with alpha and b given                        *
************************************************************************
   10 alpha = par1
      b     = par2
      alpha1= alpha+1.D0
      logC  = alpha1*dlog(b)-gammlnCP(alpha1)
      do 11 i=1, numr
          nwithr(i) = dexp(logC+alpha*dlog(r(i))-b*r(i))
   11 continue
      goto 999
************************************************************************
*  2 TWO PARAMETER GAMMA with par1= reff and par2= veff given          *
************************************************************************
   20 alpha = 1.D0/par2 - 3.D0
      b     = 1.D0/(par1*par2)
      alpha1= alpha+1.D0
      logC  = alpha1*dlog(b)-gammlnCP(alpha1)
      do 21 i=1, numr
          nwithr(i) = dexp(logC+alpha*dlog(r(i))-b*r(i))
   21 continue
      goto 999
************************************************************************
*  3 BIMODAL GAMMA with equal mode weights                             *
************************************************************************
   30 alpha = 1.D0/par3 - 3.D0
      b1    = 1.D0/(par1*par3)
      b2    = 1.D0/(par2*par3)
      gamlna= gammlnCP(alpha+1.D0)
      logC1 = (alpha+1.D0)*dlog(b1)-gamlna
      logC2 = (alpha+1.D0)*dlog(b2)-gamlna
      do 31 i=1, numr
          nwithr(i) = 0.5D0*( dexp(logC1+alpha*dlog(r(i))-b1*r(i))
     +                      + dexp(logC2+alpha*dlog(r(i))-b2*r(i)) )
   31 continue
      goto 999
************************************************************************
*  4 LOG NORMAL with rg and sigma given                                *
************************************************************************
   40 flogrg = dlog(par1)
      flogsi = dabs(dlog(par2))
      C      = 1.D0/(root2p*flogsi)
      fac    = -0.5D0/(flogsi*flogsi)
      do 41 i=1, numr
          nwithr(i) = C * dexp( fac*(dlog(r(i))-flogrg)**2 ) / r(i)
   41 continue
      goto 999
************************************************************************
*  5 LOG NORMAL with reff and veff given                               *
************************************************************************
   50 rg     = par1/(1.D0+par2)**2.5D0
      flogrg = dlog(rg)
      flogsi = dsqrt(dlog(1.D0+par2))
      C      = 1.D0/(root2p*flogsi)
      fac    = -0.5D0/(flogsi*flogsi)
      do 51 i=1, numr
          nwithr(i) = C * dexp( fac*(dlog(r(i))-flogrg)**2 ) / r(i)
   51 continue
      goto 999
************************************************************************
*  6 POWER LAW                                                         *
************************************************************************
   60 alpha = par1
      rmin  = par2
      rmax  = par3
* DEBUG
*     write(7,*) ' sizedis : power law with alpha = ',alpha
*     write(7,*) '                          rmin  = ',rmin
*     write(7,*) '                          rmax  = ',rmax
* END DEBUG
      if (dabs(alpha+1.D0) .lt. 1.d-10) then
          C = 1.D0/dlog(rmax/rmin)
      else
          alpha1 = alpha-1.d0
          C = alpha1 * rmax**alpha1 / ((rmax/rmin)**alpha1-1.d0)
      endif
      do 61 i=1, numr
          if ((r(i) .lt. rmax) .and. (r(i) .gt. rmin)) then
              nwithr(i) = C*r(i)**(-alpha)
          else
              nwithr(i) = 0.D0
          endif
   61 continue
      goto 999
************************************************************************
*  7 MODIFIED GAMMA with alpha, rc and gamma given                     *
************************************************************************
   70 alpha = par1
      rc    = par2
      gamma = par3
      b     = alpha / (gamma*rc**gamma)
      aperg = (alpha+1.D0)/gamma
      logC  = dlog(gamma) + aperg*dlog(b) - gammlnCP(aperg)
      do 71 i=1, numr
          nwithr(i) = dexp( logC + alpha*dlog(r(i)) - b*r(i)**gamma )
   71 continue
      goto 999
************************************************************************
*  8 MODIFIED GAMMA with alpha, b and gamma given                      *
************************************************************************
   80 alpha = par1
      b     = par2
      gamma = par3
      aperg = (alpha+1.D0)/gamma
      logC  = dlog(gamma) + aperg*dlog(b) - gammlnCP(aperg)
      do 81 i=1, numr
          nwithr(i) = dexp( logC + alpha*dlog(r(i)) - b*r(i)**gamma )
   81 continue
      goto 999
************************************************************************
*  9 DISCRETE VALUED DISTRIBUTION                                      *
************************************************************************
   90 do 95 k=1, ndis
          rdisp(k)   = rdis(iparts,k)
          nwrdisp(k) = nwrdis(iparts,k)
   95 continue
      call splineCP(rdisp,nwrdisp,ndis,1.d30,1.d30,y2,NDdis,NDdis)
      do 91 j=1,numr
         call splintCP(rdisp,nwrdisp,y2,ndis,r(j),nwithr(j),NDdis,NDdis)
   91 continue
      goto 999
************************************************************************
  999 if (numr .le. 1) return
* DEBUG
*     write(7,*) ' sizedis:'
*     write(7,*) '   i             r(i)               n(r(i))'
*     write(7,*) ' ----------------------------------------------------'
*     do 1000 i=1, numr
*         write(7,1001) i,r(i),nwithr(i)
*1000 continue
*1001 format(i4,1pe24.14,1pe22.14)
* END DEBUG
      return
      end


      function gammlnCP(xarg)
************************************************************************
*  Return the value of the natural logarithm of the gamma function.    *
*  The argument xarg must be real and positive.                        *
*  This function is documented in :                                    *
*                                                                      *
*  W.H. Press et al. 1986, 'Numerical Recipes' Cambridge Univ. Pr.     *
*  page 157 (ISBN 0-521-30811)                                         *
*                                                                      *
*  When the argument xarg is between zero and one, the relation (6.1.4)*
*  on page 156 of the book by Press is used.                           *
*                                         V.L. Dolman April 18 1989    *
************************************************************************
      implicit double precision (a-h,o-z)
      parameter( eps = 1.d-7, one = 1.D0, two = 2.D0, half = 0.5D0
     +         , fpf = 5.5D0 )
      dimension cof(6)
      data cof,stp/76.18009173D0,-86.50532033D0, 24.01409822D0
     +   ,-1.231739516D0, 0.120858003D-2, -0.536382D-5, 2.50662827465D0/
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
      do 11 j=1, 6
          x = x+one
          ser = ser+cof(j)/x
   11 continue
      gtmp = tmp+dlog(stp*ser)
      if  (xarg .gt. one) then
          gammlnCP = gtmp
      else
          pix = pi*(one-xarg)
          gammlnCP = dlog(pix/dsin(pix))-gtmp
      endif
      return
      end


      subroutine splineCP(x, y, n, yp1, ypn, y2,NDmui,nn)
***********************************************************************
** Spline interpolation routine from Press et al. (1986, p.88).      **
**                                                                   **
** Given arrays x and y of length n containing a tabulated function, **
** i.e. y(i)=f(x(i)), with x(1)<x(2)<...<x(n), and given values yp1  **
** and ypn for the first derivative of the interpolating function at **
** points 1 and n respectively, this routine returns an array y2 of  **
** length n which contains the second derivatives of the interpola-  **
** ting function at the tabulated points x(i).                       **
** If yp1 and/or yp2 are equal to 1x10^30 or larger, the routine is  **
** signalled to set the corresponding boundary condition for a natu- **
** ral spline, with zero second derivative on that boundary.         **
***********************************************************************
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
      do 11 i=2,n-1
        sig=(x(i)-x(i-1))/(x(i+1)-x(i-1))
        p=sig*y2(i-1)+2.d0
        y2(i)=(sig-1.d0)/p
        u(i)=(6.d0*((y(i+1)-y(i))/(x(i+1)-x(i))-(y(i)-y(i-1))
     +       /(x(i)-x(i-1)))/(x(i+1)-x(i-1))-sig*u(i-1))/p
   11 continue
      if (ypn.gt..99d30) then
        qn=0.d0
        un=0.d0
      else
        qn=0.5d0
        un=(3.d0/(x(n)-x(n-1)))*(ypn-(y(n)-y(n-1))/(x(n)-x(n-1)))
      endif
      y2(n)=(un-qn*u(n-1))/(qn*y2(n-1)+1.d0)
      do 12 k=n-1,1,-1
        y2(k)=y2(k)*y2(k+1)+u(k)
   12 continue
      return
      end
      subroutine splintCP(xa, ya, y2a, n, x, y,NDmui,nn)
***********************************************************************
** Spline interpolation routine from Press et al. (1986, p.88).      **
**                                                                   **
** Given the arrays xa and ya of length n, which tabulate a function **
** (with the xa(i)'s in order), and given the array y2a, which is    **
** the output from splineCP above, and given a value of x, this        **
** routine returns a cubic-spline interpolated value y.              **
***********************************************************************
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
      y=a*ya(klo)+b*ya(khi)+
     +  ((a**3-a)*y2a(klo)+(b**3-b)*y2a(khi))*(h**2)/6.d0
      return
   10 format('ERROR in routine splintCP: Bad XA input.')
      end
      subroutine tstsym( develop, thmin, thmax, step, symth )
************************************************************************
*  Test if the set of theta points is symmetrical around 90 degrees    *
*  and return this information through logical symth                   *
*  In case of development in GSF we have a symmetrical Gauss set !!    *
************************************************************************
      implicit double precision(a-h,o-z)
      parameter( eps=1.d-6, heps=0.5d0*eps )
      integer develop
      logical symth
      symth = .false.
      if (develop .eq. 1) then
          symth = .true.
      else
          if ( (dabs( 180.d0 - thmin - thmax ) .lt. eps)  .and.
     +         (dmod( thmax-thmin+heps, step ) .lt. eps) )
     +              symth = .true.
      endif
* DEBUG
c      if (symth) then
c          write(7,*) ' tstsym: theta points symmetrical'
c      else
c          write(7,*) ' tstsym: theta points NOT symmetrical !!'
c      endif
* END DEBUG
      return
      end
