!------------------------------------------------------
! This routine calculates optical cross sections of
! fractal dust aggregates using a modified-mean-field 
! approximation.
!                               2017. 07.24. Ryo Tazaki
!------------------------------------------------------
! INPUT AND OUTPUT
!------------------------------------------------------
!
! Input parameters:
!
! lmd   : wavelength (um)
! R0    : monomer radius (um)
! Df    : fractal dimension
! PN    : number of monomers
! k0    : fractal prefactor PN = k0 (Rg/R0)^Df
! refrel: complex refractive index
! nang  : Number of angular grid between 0 ang 90 degrees.
!         If nang=91, then angular width is 1 degrees.
!
! iqsca : Switch for calculation method
!         iqsca = 1: Rayleigh-Gans-Debye theory
!                    Reference: Tazaki et al. 2016
!         iqsca = 2: Mean-field theory 
!                    Reference: Botet et al. 1997; Rannou et al. 1997
!         iqsca = 3: Modified mean-field theory
!                    Reference: Tazaki & Tanaka 2018
!
! iqcor : Switch for correlation function to describe aggregate
!         structure.
!         iqcor = 1: Gaussian cutoff (Tazaki et al. 2016)
!         iqcor = 2: Exponential cutoff (Berry & Percival 1986)
!         iqcor = 3: fractal dimension cutoff (Botet et al. 1995)
!         iqcor = 4: Homogeneous sphere 
!
! Output:
!
! Cext  : extinction cross section (um^2)
! Csca  : scattering cross section (um^2)
! Cabsp : absorption cross section (um^2)
! Gsca  : asymmetry parameter
! Smat  : Scattering matrix elements S_ij
!         The definition of S_ij obeys Bohren & Huffman's one.
!
!------------------------------------------------------
! Revision history
!------------------------------------------------------
!
! version 1.0 (July 24 2017)
! - Initial release 
!
! version 2.2 (Oct. 24 2017)
! - Implemented geometrical optics approximation (Modified Mean Field; MMF)
! - Implemented a stable algorithm for computations of the spherical Bessel function 
!       of the first kind at large order and/or large argument.
! - Fixed bug of Legendre polynomial calculations with large order.
! - Implemented dynamic allocation for some special function calculations.
! - Fixed bug about compiler dependence (ifort or gfortran)
!
! version 2.3 (Dec. 05 2017)
! - Implemented efficient computation of $a(\nu,n,p)$ and $b(\nu,n,p)$ using the parity arguments.
!
! version 2.5 (May. 31 2020)
! - Fixed a bug about initialization of Qsca and gsca in (renormalized) Lorentz-Mie routine.
!   This does not change the results, since both of them are not used to obtain final results.
! 
! version 2.5.1 (Oct. 19, 2020) 
!  - some minor changes to adjust OpTool
!    OpTool: https://github.com/cdominik/optool
!
!------------------------------------------------------
subroutine meanscatt(lmd,R0,PN,Df,k0,refrel,iqsca,iqcor,nang,&
                Cext,Csca,Cabsp,GSCA,Smat)
implicit none
double precision,parameter::pi=acos(-1.0d0)
double precision::PN
double precision::Rg,R0,df,k,xi,xg,x0,c,lmd,k0,xstop
integer::nstop,numax,nmax
integer::iqsca,iqcor
integer::nu,n,p,j,jm,pmin,pmax

complex(kind(0d0))::refrel
complex(kind(0d0)),allocatable,dimension(:)::an,bn
complex,allocatable,dimension(:)::y,r
complex,allocatable,dimension(:,:)::a,d,d4chk
complex(kind(0d0)),allocatable,dimension(:,:)::ad,dd
complex,allocatable,dimension(:,:)::S
complex,allocatable,dimension(:,:,:)::T

real,allocatable,dimension(:)::x,w
real::x1,x2,dx
complex(kind(0d0))::sumA,sumB
complex(kind(0d0))::Sp_tmp
complex(kind(0d0)),allocatable,dimension(:)::Sp
double precision::anunp,bnunp,itgr_a,itgr_b
double precision::gaunta,gauntb

double precision::Cext,Csca,Gsca,Smat(1:4,2*nang-1)

integer::nang
complex(kind(0d0)),allocatable,dimension(:)::S1,S2
double precision::Qext0,Qsca0,GSCA0

DOUBLE PRECISION::S11,S12,S33,S34
double precision::dang
double precision,allocatable,dimension(:)::ang,PF!,Isca
double precision,allocatable,dimension(:)::S11A,S12A,S33A,S34A
double precision::rad,q,AA,Sq
double precision,parameter::d2r=pi/180.0d0
double precision,parameter::r2d=1.0d0/d2r

double precision::hg
double precision::al,bb,xx,cc
double precision::CabspA,Cabsp,cn1,cn2
double precision::nrm,g_asym
double precision::GC,GS,tau

integer::mxnang
integer::jj
double precision::angip,S11ip,PFip

double precision::fr,fc,u,Rc,kk,xxg,xxc
integer::itest

Rg = R0 * (dble(PN)/k0) ** (1.0d0/df)
k   = 2.0d0*pi/lmd
xg = k * Rg

!nang = 90*3+1

allocate(S1(1:2*NANG-1),S2(1:2*NANG-1))
allocate(ang(1:2*NANG-1),PF(1:2*NANG-1))
allocate(S11A(1:2*NANG-1),S12A(1:2*NANG-1),S33A(1:2*NANG-1),S34A(1:2*NANG-1))

x0 = 2.0d0*pi*R0/lmd
xstop= x0 + 4.0d0 * x0 ** (1.0/3.0) + 2.0d0

nstop=nint(xstop)
numax = nstop
nmax  = nstop

!write(*,fmt='(A6,1PE12.5,A6,1PE12.5)',advance='no') " x0 = ",x0,", xg = ",xg
!
!safety check
!
if(numax + nmax .ge. 500) then
        write(*,*) "Monomer's size parameter...",x0
        write(*,*) "Monomer expansion order...",nstop
        write(*,*) "Maximum order of sph. Besssel : ",nstop*2
        write(*,*) "Warning: p = 2 * monomer expand. > 500"
        write(*,*) "This parameter space was NOT tested!"
        write(*,*) "Calculation is aborted!"
        stop
endif

if (xg .ge. 1.0d6) then
        write(*,*) "Safety Stop: "
        write(*,*) "Aggregate size parameter exceeds 1d6, where"
        write(*,*) "the accuracy of numerical integration of the s_p"
        write(*,*) "factor was not checked."
        write(*,*) "Calculation is aborted !"
        stop
endif


! output correlation function
!        open(50,file="gu_df3.0.out",status="unknown")
!        do itest=1,500
!         !u=2.0*R0*(2000.0d0/2.0d0)**(dble(itest-1)/dble(99))
!         kk=1.0d0
!         u=1.0e-3*(1.0e1/1.0e-3)**(dble(itest-1)/dble(499))
!         Rc=Rg*sqrt(5.0d0/3.0d0)
!         u=u*Rc
!         xxg=Rg*kk
!         
!         fr=1.0d0/(4.0*pi*Rg**3.0)*(u/Rg)**(df-3.0d0)
!         if(u .le. 2.0d0*Rc) then
!                 write(50,*) u/Rc,fr*fc(1,u,xxg,df),fr*fc(2,u,xxg,df),fr*fc(3,u,xxg,df),fr*fc(4,u,xxg,df)
!         else
!                write(50,*) u/Rc,fr*fc(1,u,xxg,df),fr*fc(2,u,xxg,df),fr*fc(3,u,xxg,df),0.0d0
!         endif
!
!        enddo
!        close(50)
!        stop


allocate(an(nstop),bn(nstop),a(2,nstop),d(2,nstop),d4chk(2,nstop))
allocate(ad(2,nstop),dd(2,nstop))

if(iqsca .ge. 2) then
allocate(T(2,numax,nmax))
allocate(S(2*numax,2*nmax))
allocate(y(2*nmax),r(2*nmax))
endif

!=============================================
! Exansion coefficients (d1n,d2n) of a monomer
! taking the multiple scattreing into account.
!============================================

!
! Lorentz-Mie coefficients for a single monomer
!

call lorentz_mie(x0,refrel,nstop,an,bn)
a = cmplx(0.0,0.0)
do n=1,nstop
        a(1,n) = cmplx(an(n)) 
        a(2,n) = cmplx(bn(n))
        ad(1,n) = an(n)
        ad(2,n) = bn(n)
enddo


if(iqsca .eq. 1) then


        !Single scattering (only external radiation field)
        do n=1,nmax
                d(1,n) = a(1,n)
                d(2,n) = a(2,n)
                dd(1,n)= ad(1,n)
                dd(2,n)= ad(2,n)
        enddo


elseif(iqsca .ge. 2) then


        !
        ! integral for a(nu,n,p) and b(nu,n,p)
        !
        
        jm = 400
        allocate(x(jm),w(jm))
        x = 0.0
        w = 0.0
        x1 = -1.0
        x2 = 1.0
        
                !Gauss-legendre quadrature
                call gauleg(x1,x2,x,w,jm) 

                !daikei
                !call daikei(x1,x2,x,jm) 

        !
        ! Storing the value of Sp(kRg)
        !
        !call spout(iqcor,df)

        allocate(Sp(0:numax+nmax))
        Sp = cmplx(0.0d0,0.0d0,kind(0d0))
        !write(*,*) numax+nmax
        !stop
        do p=0,numax+nmax
              Sp_tmp = cmplx(0.0d0,0.0d0,kind(0d0))
              
              if(iqcor .eq. 4) then
              call integration_of_Sp(iqcor,3.0d0,p,k,xg,Sp_tmp)
              else
              call integration_of_Sp(iqcor,df,p,k,xg,Sp_tmp)
              endif
                      !spherical bessel test:
                      !call integration_of_Sp(iqcor,df,200,xg,Sp_tmp)
                      !stop
              Sp(p) = Sp_tmp
              !write(*,*) p,real(Sp(p)),aimag(Sp(p))
        enddo
        !stop


        !
        ! Translation coefficients: A and B
        !
        do nu=1,numax
        !do nu=numax,numax
        !do nu=24,24
        !write(*,*) nu,"/",numax
                do n=1,nmax
                !do n=nmax,nmax
                !do n=81,81
                      !write(*,*) n,"/",nmax
                      pmin=abs(n-nu)
                      pmax=nu+n

                      sumA = 0.0d0
                      sumB = 0.0d0

                      do p=pmin,pmax 
                      !do p=105,150 
                      !write(*,*) "(nu,n,p)=",nu,n,p

                                anunp=0.0d0
                                bnunp=0.0d0

                                        !S Parity (ver. 2.2):
                                        !Considering parity of a and b for fast computation.
                                        !if p has same parity as n+nu, a(nu,n,p) .ne. 0 and b(nu,n,p) = 0
                                        !if p has not same parity as n+nu, b(nu,n,p) .ne. 0 and a(nu,n,p) = 0
                                        if(mod(pmax,2) .eq. mod(p,2)) then
                                                do j=1,jm
                                                   anunp = anunp + w(j) * itgr_a(nu,n,p,x(j))
                                                enddo
                                        else
                                                do j=1,jm
                                                   bnunp = bnunp + w(j) * itgr_b(nu,n,p,x(j))
                                                enddo
                                        endif
                                        !E: Parity

                                        !v2.1: p runs all integers between |n-nu| and n+nu
                                        !But, this calculation includes some unneccesary calc.
                                                !gauss-legendre quadrature
                                                !do j=1,jm
                                                !   anunp = anunp + w(j) * itgr_a(nu,n,p,x(j))
                                                !   bnunp = bnunp + w(j) * itgr_b(nu,n,p,x(j))
                                                !enddo
                                                !daikei
                                                !do j=1,jm-1
                                                !   anunp = anunp + 0.5d0 * (itgr_a(nu,n,p,x(j))+&
                                                !           itgr_a(nu,n,p,x(j+1)))*(x(j+1)-x(j))
                                                !   bnunp = bnunp + 0.5d0 * (itgr_b(nu,n,p,x(j))+&
                                                !           itgr_b(nu,n,p,x(j+1)))*(x(j+1)-x(j))
                                                !enddo

                                anunp = ((2.0d0*dble(p)+1.0d0)/2.0d0) * anunp
                                bnunp = ((2.0d0*dble(p)+1.0d0)/2.0d0) * bnunp
                               
                                !Sp = Sp_store(p)
                                !Sp = L/DBLE(PN-1)

                                sumA = sumA + (dble(n)*dble(n+1)+dble(nu)*dble(nu+1)-dble(p)*&
                                        & dble(p+1)) * anunp * Sp(p) 
                               
                                sumB = sumB + bnunp * Sp(p)
                      enddo

                        T(1,nu,n) = ((2.0*real(nu)+1.0)/(real(n)*real(n+1)*&
                                & real(nu)*real(nu+1))) * cmplx(sumA)

                        T(2,nu,n) = 2.0*((2.0*real(nu)+1.0)/(real(n)*real(n+1)*&
                                & real(nu)*real(nu+1))) * cmplx(sumB)

                enddo

        enddo

        deallocate(x,w)
        S = cmplx(0.0,0.0)

        do n=1,nmax
                do nu=1,numax
                S(2*n-1,2*nu-1) = real(PN-1) * a(1,n) * T(1,nu,n) !a*A
                S(2*n-1,2*nu)   = real(PN-1) * a(1,n) * T(2,nu,n) !a*B
                enddo

                do nu=1,numax
                S(2*n,2*nu-1) = real(PN-1) * a(2,n) * T(2,nu,n) !b*B
                S(2*n,2*nu)   = real(PN-1) * a(2,n) * T(1,nu,n) !b*A
                enddo
        enddo


        do n=1,2*nmax
               S(n,n) = 1.0 + S(n,n)
        enddo

        do n=1,nmax
                y(2*n-1) = a(1,n)
                y(2*n)   = a(2,n)
        enddo


        !***********************************************
        ! Solving for the mean field
        !
        ! This routine is checked.
        !
        ! Solving inverse matrix of 4*Nstop x 4*Nstop
        ! 4 = 2 modes x 2 complex variables of fields
        ! Stupid solution algorithm.
        !
        ! S * x = y 
        ! r = S^{-1} * y 
        !
        call complex_leqs_solver(2*nstop,4*nstop,S,y,r)
        !
        !***********************************************


        !Mean field coefficients (Renormalized Lorentz-Mie coefficients)
        do n=1,nmax
           d(1,n) = r(2*n-1) !d^(1)_n
           d(2,n) = r(2*n)   !d^(2)_n
        enddo


                !write(*,*) "-------------------------------------------------------------"
                !write(*,2000) "n","Re(d1n(1))","Im(d1n(1))","Re(d1n(2))","Im(d1n(2))"
                !do n=1,nstop
                !  write(*,1000) n,real(d(1,n)),aimag(d(1,n)),real(d(2,n)),aimag(d(2,n))
                !enddo
                !1000 format(' ',I3,1P4E15.5)
                !2000 format(' ',A3,4A15)


        !=================================================================
        ! Scattering by a fractal dust aggregate: Now (d1n,d2n) is 
        ! regarded as the new Lorentz-Mie coefficients of a single monomer.
        !===================================================================
        ! single complex --> double complex
        do n=1,nmax
           dd(1,n) = cmplx(dble(d(1,n)),aimag(d(1,n)),kind(0d0))
           dd(2,n) = cmplx(dble(d(2,n)),aimag(d(2,n)),kind(0d0))
        enddo

                !
                ! For geometrical cross-section test
                !
                !do n=1,nstop
                !        Cext = Cext + dble(2*n+1)*real(d(1,n)+d(2,n))
                !        !Monomer
                !        !Cabsp = Cabsp + dble(2*n+1)*&
                !        !       &real(cn1*cabs(a(1,n))**2.0+cn2*cabs(a(2,n))**2.0)
                !enddo
                !Cext = (2.0d0*pi*dble(PN)/k**2.0)*Cext
                !!Cabsp = (2.0d0*pi*dble(PN)/k**2.0)*Cabsp
                !return
                !!

endif !end iqsca==2

!
call renorm_mie(dd,x0,nstop,nang,S1,S2,QEXT0,QSCA0,GSCA0)
!

S11A = 0.0d0
S12A = 0.0d0
S33A = 0.0d0
S34A = 0.0d0

DANG=0.5d0*PI/DBLE(NANG-1)
DO J=1,2*NANG-1
        ANG(J)=DANG*DBLE(J-1)*(180.0d0/PI)
        
        rad=ANG(J)*d2r
        q=2.0d0*k*sin(rad/2.0d0)
 
        S11=0.5D0*ABS(S2(J))*ABS(S2(J))
        S11=S11+0.5D0*ABS(S1(J))*ABS(S1(J))
        S12=0.5D0*ABS(S2(J))*ABS(S2(J))
        S12=S12-0.5D0*ABS(S1(J))*ABS(S1(J))       
        S33=DBLE(S1(J)*CONJG(S2(J)))
        S34=AIMAG(S2(J)*CONJG(S1(J)))
        !write(*,*) ANG(J),S11,S12,S33,S34

        !---------------------------------
        !Calculate structure factor Sq
        !---------------------------------
        if(iqcor .eq. 1) then

                HG = 0.0d0
                AL = df/2.0d0
                BB = 1.5d0
                XX = -(q*Rg)**2.0d0/df
                call CHGM(AL,BB,XX,HG)
                if(HG*0.d0 /= 0.d0 .or. isnan(HG) .eqv. .true.) then
                CC=sqrt(pi)*Df**(Df/2.0)/(2.0*GAMMA((3.0-Df)/2.0))
                HG=CC*(q*RG)**(-Df)
                end if
                
                S11A(J) = dble(PN) * S11 * (1.0d0 + dble(PN-1) * HG)
                S12A(J) = dble(PN) * S12 * (1.0d0 + dble(PN-1) * HG)
                S33A(J) = dble(PN) * S33 * (1.0d0 + dble(PN-1) * HG)
                S34A(J) = dble(PN) * S34 * (1.0d0 + dble(PN-1) * HG)

         elseif(iqcor .eq. 2) then

                xi =  sqrt(2.0d0/(df*(df+1.0d0))) * Rg
               
                IF(J .eq. 1) then
                        Sq = 1.0d0
                else
                        Sq = sin((df-1.0d0)*atan(q*xi))/(q*xi*(df-1.0d0)&
                             *(1.0d0+(q*xi)**2.0)**((df-1.0d0)/2.0))
                endif
                
                S11A(J) = dble(PN) * S11 * (1.0d0 + dble(PN-1) * Sq)
                S12A(J) = dble(PN) * S12 * (1.0d0 + dble(PN-1) * Sq)
                S33A(J) = dble(PN) * S33 * (1.0d0 + dble(PN-1) * Sq)
                S34A(J) = dble(PN) * S34 * (1.0d0 + dble(PN-1) * Sq)

         elseif(iqcor .eq. 3) then

                c  = 0.5d0   ! Reference: Botet et al. (1995)
                AA = 4.0d0 * pi 
                AA = c*df/AA

                call structure_factor_integration(c,df,q,Rg,xg,Sq)
                
                if(j .eq. 1) then
                        Sq = 1.0d0
                else
                        Sq = (4.0d0*pi*aa/(q*rg)) * Sq
                endif

                S11A(J) = dble(PN) * S11 * (1.0d0 + dble(PN-1) * Sq)
                S12A(J) = dble(PN) * S12 * (1.0d0 + dble(PN-1) * Sq)
                S33A(J) = dble(PN) * S33 * (1.0d0 + dble(PN-1) * Sq)
                S34A(J) = dble(PN) * S34 * (1.0d0 + dble(PN-1) * Sq)

          elseif(iqcor .eq. 4) then !RG structure factor

                if(j .eq. 1) then
                        Sq = 1.0d0
                else
                       ! Sq = df*(sin(q*Rg)-q*Rg*cos(q*Rg))/(q*Rg)**3.0
                        Sq = (3.0d0*(sin(q*Rg)-q*Rg*cos(q*Rg))/(q*Rg)**3.0)**2.0d0
                       ! Sq = (3.0d0*(sin(q*sqrt(5.0/3.0)*Rg)-q*sqrt(5.0/3.0)*Rg&
                       !         *cos(q*sqrt(5.0/3.0)*Rg))/(q*sqrt(5.0/3.0)*Rg)**3.0)**2.0d0
                endif

                !N ---> infty.
                !S11A(J) = dble(PN) ** 2.0 * S11 * Sq
                !S12A(J) = dble(PN) ** 2.0 * S12 * Sq
                !S33A(J) = dble(PN) ** 2.0 * S33 * Sq
                !S34A(J) = dble(PN) ** 2.0 * S34 * Sq

                S11A(J) = dble(PN) * S11 * (1.0d0 + dble(PN-1) * Sq)
                S12A(J) = dble(PN) * S12 * (1.0d0 + dble(PN-1) * Sq)
                S33A(J) = dble(PN) * S33 * (1.0d0 + dble(PN-1) * Sq)
                S34A(J) = dble(PN) * S34 * (1.0d0 + dble(PN-1) * Sq)

          endif

enddo

!***********************************************
!   SCATTERING CROSS-SECTION VIA SIMPSON FORMULA
!***********************************************

DANG=0.5d0*pi/DBLE(NANG-1)
Csca=0.0d0
DO J=1,NANG-1
        Csca=Csca+DANG*(S11A(2*J-1)*sin(DANG*DBLE(2*J-2))&
        + 4.0d0*S11A(2*j)*sin(DANG*DBLE(2*J-1))&
        + S11A(2*j+1)*sin(DANG*DBLE(2*J)))/3.0d0 
END DO
Csca=2.0*pi*Csca/k**2.0


!***********************************************
!  PHASE FUNCTION
!  Definition : Bohren & Huffman textbook.
!***********************************************
PF=0.0d0
!open(10,file="PF_lmd4um.dat",status="unknown")
DO J=1,2*NANG-1
        PF(J) = S11A(J)/(CSCA*k**2.0)
!        write(10,*) ANG(J),S11A(J),PF(J)
END DO
!close(10)

5100    format(' ',1P6E15.5)

!***********************************************
!  ASYMMETRY PARAMETER: g = <cos\theta>
!***********************************************
g_asym=0.0d0 
nrm = 0.0d0

DO J=1,NANG-1
        g_asym=g_asym+DANG*(PF(2*J-1)*sin(DANG*DBLE(2*J-2))*cos(DANG*DBLE(2*J-2))&
        +4.0d0*PF(2*j)*sin(DANG*DBLE(2*J-1))*cos(DANG*DBLE(2*J-1))&
        +PF(2*j+1)*sin(DANG*DBLE(2*J))*cos(DANG*DBLE(2*J)))/3.0d0 

        nrm=nrm+DANG*(PF(2*J-1)*sin(DANG*DBLE(2*J-2))&
        +4.0d0*PF(2*j)*sin(DANG*DBLE(2*J-1))&
        +PF(2*j+1)*sin(DANG*DBLE(2*J)))/3.0d0 
END DO

nrm=2.0d0*pi*nrm
g_asym=2.0*pi*g_asym


!
!Check normalization of phase function
!
if ((nrm-1.0d0) .gt. 1.0d-2) then
        write(*,*) ABS(nrm-1.0)*100
        write(*,*) "ERROR : PHASE FUNCTION IS NOT NORMALIZED TO UNITY"
        write(*,*) "CALCULATION ABORTED !"
        stop
end if


GSCA=g_asym

!***********************************************
!  ABSORPTION AND EXTINCTION CROSS SECTIONS
!***********************************************
if(iqsca .eq. 1) then
        Cabsp= 0.0d0
        do n=1,nstop
                cn1=real(1.0/conjg(a(1,n))-1.0)
                cn2=real(1.0/conjg(a(2,n))-1.0)
                Cabsp = Cabsp + dble(2*n+1)*&
                       &real(cn1*cabs(d(1,n))**2.0+cn2*cabs(d(2,n))**2.0)
        enddo
        Cabsp= (2.0d0*pi*dble(PN)/k**2.0)*cabsp
        Cext = Cabsp + Csca
elseif(iqsca .eq. 2) then
        Cext = 0.0d0
        Cabsp= 0.0d0
        CabspA = 0.0d0
        do n=1,nstop
                Cext = Cext + dble(2*n+1)*real(d(1,n)+d(2,n))
                cn1=real(1.0/conjg(a(1,n))-1.0)
                cn2=real(1.0/conjg(a(2,n))-1.0)
        enddo
        Cext = (2.0d0*pi*dble(PN)/k**2.0)*Cext
        Cabsp = Cext - Csca        
elseif(iqsca .eq. 3) then
        Cext = 0.0d0
        Cabsp= 0.0d0
        do n=1,nstop
                Cext = Cext + dble(2*n+1)*real(d(1,n)+d(2,n))
                cn1=real(1.0/conjg(a(1,n))-1.0)
                cn2=real(1.0/conjg(a(2,n))-1.0)
                ! Monomer's Cabs here: 
                Cabsp = Cabsp + dble(2*n+1)*&
                       &real(cn1*cabs(a(1,n))**2.0+cn2*cabs(a(2,n))**2.0)
        enddo
        Cext = (2.0d0*pi*dble(PN)/k**2.0)*Cext
        Cabsp= (2.0d0*pi*dble(pn)/k**2.0)*cabsp
        RC=SQRT(5.0/3.0)*RG
        GC=pi*RC**2.0
        GS = GC 
        tau = (Cabsp/GS) !: Cabsp in RGD
        Cabsp = GS * (1.0d0 - dexp(-tau))
        Cabsp = max(Cabsp,Cext-Csca)
        Csca = Cext - Cabsp
endif

!Scatteirng matrix return
do j=1,2*nang-1
Smat(1,j)=S11A(j)
Smat(2,j)=S12A(j)
Smat(3,j)=S33A(j)
Smat(4,j)=S34A(j)
enddo

deallocate(an,bn,a,d,d4chk,ad,dd,&
        & S1,S2,ang,S11A,S12A,S33A,S34A,PF)

if(iqsca .ge. 2) then
        deallocate(T,S,y,r,Sp)
endif

return
end subroutine meanscatt


!
subroutine daikei(x1,x2,x,jm)
implicit none
integer::i,jm
real::dx,x1,x2
real,dimension(jm)::x

dx=(x2-x1)/dble(jm-1)
do i=1,jm
        x(i) = x1 + dble(i-1) * dx
enddo

return
end subroutine daikei

!
subroutine structure_factor_integration(c,D,q,Rg,xg,Sq)

implicit none
integer::m,i,n
double precision::x,c,D,q,Rg,s,xg,Sq,h
double precision::func,a,b

!
!Set integration range
!
a = 0.0d0 
b = 20.0

!a = log10(1.0d-3)
!b = log10(1.0d2)

! this is purely empiriacal.
!m=250
m = 1000 * int(1.0+xg**(1.0/2.0))
!

Sq = 0.0d0
s = 0.0d0 
n = 2 * m
h = (b-a)/dble(n)
x = a
do i=0,m-1
s = h*(func(x+2*i*h,c,D,q,Rg)+&
        & 4*func(x+(2*i+1)*h,c,D,q,Rg)&
        & +func(x+(2*i+2)*h,c,D,q,Rg))/3.0d0
Sq = Sq + s
end do


return
end subroutine structure_factor_integration


!
! Integrand of structure factor integration
!
function func(x,c,D,q,Rg)
implicit none
double precision::q,Rg,D,c,x
double precision::func

if(x .eq. 0.0d0) then
        func = 0.0d0
else
        func = dsin(q*Rg*x)*dexp(-c*x**D)*x**(D-2.0d0)
endif

return
end function func


!---------------------------------------------
! Calculate Lorentz-Mie coefficients: an, bn
! This routine is based on BHMIE code written 
! Bruce Draine, although this is in double precision.
!---------------------------------------------
subroutine lorentz_mie(x,refrel,nstop,a,b)
implicit none
integer,parameter::nmxx=1500000
!input variables
integer::nstop
double precision::x
complex(kind(0d0))::refrel

!Local variables
integer::nmx
complex(kind(0d0)),dimension(nmxx)::d
double precision::xstop,ymod
complex(kind(0d0))::y

integer::n,en
double precision::psi0,psi1,psi,chi0,chi1,chi
complex(kind(0d0))::xi1,xi
!complex(kind(0d0))::an,bn
complex(kind(0d0)),dimension(nstop)::a,b

!
y = refrel*x
ymod = abs(y)
xstop= x + 4.0d0 * x ** (1.0/3.0) + 2.0d0
!nstop=nint(xstop)
nmx = nint(max(xstop,ymod)) + 15

if(nmx .gt. nmxx) then
        write(*,*) "Error: nmx > nmxx=",nmxx,"for |m|x=",ymod
        stop
endif

!
!calculate D_n(mx) by downward recurrence.
!Initial value is set as D(mx)=0+0i at n=NMX.
!
d(nmx)=cmplx(0.0d0,0.0d0,kind(0d0))
do n=1,nmx-1
        en=nmx-n+1
        d(nmx-n) = (dble(en)/y) - (1.0d0/(d(en) + dble(en)/y))
enddo

!
psi0=cos(x)
psi1=sin(x)
chi0=-sin(x)
chi1=cos(x)
xi1=cmplx(psi1,-chi1,kind(0d0))

do n=1,nstop

        !Calculate psi and chi via upward recurrence:
        psi = (2.0d0*dble(n)-1.0d0)*psi1/x - psi0
        chi = (2.0d0*dble(n)-1.0d0)*chi1/x - chi0
        xi = cmplx(psi,-chi,kind(0d0))

        !Calculate Lorentz-Mie coefficients:
        !an=(d(n)/refrel+dble(n)/x)*psi-psi1
        !an=an/((d(n)/refrel+(dble(n)/x))*xi-xi1)
        !bn=(refrel*d(n)+dble(n)/x)*psi-psi1
        !bn=bn/((refrel*d(n)+(dble(n)/x))*xi-xi1)

        a(n)=(d(n)/refrel+dble(n)/x)*psi-psi1
        a(n)=a(n)/((d(n)/refrel+(dble(n)/x))*xi-xi1)
        b(n)=(refrel*d(n)+dble(n)/x)*psi-psi1
        b(n)=b(n)/((refrel*d(n)+(dble(n)/x))*xi-xi1)

        !write(*,*) n,real(an),aimag(an),real(bn),aimag(bn)

        !Preparation for next step
        psi0 = psi1
        psi1 = psi
        chi0 = chi1
        chi1 = chi
        xi1 = cmplx(psi1,-chi1,kind(0d0))

enddo

return
end subroutine lorentz_mie

!--------------------------------------------------
! Renormalizing mie code: 
! Modified from BHMIE code written by Bruce Draine.
!--------------------------------------------------
subroutine renorm_mie(d,x,nstop,nang,s1,s2,qext,qsca,gsca)
implicit none
integer::nang,j,jj,n,nn,nstop
double precision::x,qext,qsca,gsca
complex(kind(0d0)),dimension(1:2*nang-1)::S1,S2
double precision::dang,pii,theta
double precision,dimension(1:nang)::amu,pi0,pi1,pi,tau
complex(kind(0d0)),dimension(2,nstop)::d
!complex(kind(0d0)),dimension(2,nstop)::d
complex(kind(0d0))::an,bn,an1,bn1

double precision::en,fn,p

double complex dpcx
double precision realpart
double precision imagpart
realpart(dpcx)=(dble(dpcx))
imagpart(dpcx)=(dimag(dpcx))


!*** require nang.ge.1 in order to calculate scattering intensities

      pii=4.d0*atan(1.d0)
      dang=0.
      if(nang.gt.1)dang=.5*pii/dble(nang-1)
      do j=1,nang
         theta=dble(j-1)*dang
         amu(j)=cos(theta)
      enddo
      do j=1,nang
         pi0(j)=0.
         pi1(j)=1.
      enddo
      nn=2*nang-1
      !double complex (dcxs1,dcxs2)
      do j=1,nn
         !dcxs1(j)=(0.d0,0.d0)
         !dcxs2(j)=(0.d0,0.d0)
         S1(j)=cmplx(0.d0,0.d0,kind(0d0))
         S2(j)=cmplx(0.d0,0.d0,kind(0d0))
      enddo


    p=-1.0d0  

    an=cmplx(0.0d0,0.0d0,kind(0d0))
    an1=cmplx(0.0d0,0.0d0,kind(0d0))
    bn=cmplx(0.0d0,0.0d0,kind(0d0))
    bn1=cmplx(0.0d0,0.0d0,kind(0d0))

    ! added in v2.5
    qsca = 0.d0
    gsca = 0.d0
    !

    do n=1,nstop
         en = dble(n)
         fn = (2.0d0*en+1.0d0)/(en*(en+1.0d0))
         !en=n
         !fn=(2.e0*en+1.)/(en*(en+1.))

!*** store previous values of an and bn for use
!    in computation of g=<cos(theta)>

         if(n.gt.1)then
            an1=an
            bn1=bn
         endif

!*** compute an and bn:
     !renormalizing the lorentz-mie coeff.
         an = d(1,n)
         bn = d(2,n)

!*** augment sums for qsca and g=<cos(theta)>

         qsca=qsca+(2.0d0*en+1.0d0)*(abs(an)**2.0+abs(bn)**2.0)
         gsca=gsca+((2.0d0*en+1.0d0)/(en*(en+1.0d0)))*&
     &        (realpart(an)*realpart(bn)+imagpart(an)*imagpart(bn))
         if(n.gt.1)then
            gsca=gsca+((en-1.0d0)*(en+1.0d0)/en)*&
     &      (realpart(an1)*realpart(an)+imagpart(an1)*imagpart(an)+&
     &      realpart(bn1)*realpart(bn)+imagpart(bn1)*imagpart(bn))
         endif


!c*** now calculate scattering intensity pattern
!    first do angles from 0 to 90

         do j=1,nang
            jj=2*nang-j
            pi(j)=pi1(j)
            tau(j)=en*amu(j)*pi(j)-(en+1.0d0)*pi0(j)

            S1(j)=S1(j)+fn*(an*pi(j)+bn*tau(j))
            S2(j)=S2(j)+fn*(an*tau(j)+bn*pi(j))
         enddo

!*** now do angles greater than 90 using pi and tau from
!    angles less than 90.
!    p=1 for n=1,3,...; p=-1 for n=2,4,...

         p=-p
         do j=1,nang-1
            jj=2*nang-j
            S1(jj)=S1(jj)+fn*p*(an*pi(j)-bn*tau(j))
            S2(jj)=S2(jj)+fn*p*(bn*pi(j)-an*tau(j))
         enddo

!*** compute pi_n for next value of n
!    for each angle j, compute pi_n+1
!    from pi = pi_n , pi0 = pi_n-1

         do j=1,nang
            pi1(j)=((2.0d0*en+1.0d0)*amu(j)*pi(j)-(en+1.0d0)*pi0(j))/en
            pi0(j)=pi(j)
         enddo

      enddo

!*** have summed sufficient terms.
!    now compute qsca,qext,qback,and gsca
      gsca=2.d0*gsca/qsca
      qsca=(2.d0/(x*x))*qsca
      qext=(4.d0/(x*x))*realpart(S1(1))

return
end subroutine renorm_mie


!
!-----------------------------
! Integrand of a(nu,n,p)
!------------------------------
function itgr_a(nu,n,p,xr)
implicit none
integer::nu,n,p
real::xr
!local variables
double precision::x,P1nu,P1n,Pp
double precision,allocatable,dimension(:,:)::PMN,PMND
double precision,allocatable,dimension(:)::PN,PND

double precision::itgr_a

x = dble(xr)

!Calculate associated Legendre polynominals
!call lpmn(100,1,max(nu,n),X,PMN,PMND)

allocate(PMN(0:10,0:max(nu,n,10)),PMND(0:10,0:max(nu,n,10)))
allocate(PN(0:max(p,1)),PND(0:max(p,1)))

call lpmn(10,1,max(nu,n,10),X,PMN,PMND)
P1nu = PMN(1,nu)
P1n = PMN(1,n)

call lpn(max(p,1),X,PN,PND)
!call lpn(p,X,PN,PND)
Pp  = PN(p)

deallocate(PMN,PMND,PN,PND)

itgr_a = P1nu*P1n*Pp

return
end function itgr_a


!-----------------------------
! Integrand of b(nu,n,p)
!------------------------------
function itgr_b(nu,n,p,xr)
implicit none
integer::nu,n,p
real::xr

!local variables
double precision::x,P1nu,P1n,PDp
double precision,allocatable,dimension(:,:)::PMN,PMND
double precision,allocatable,dimension(:)::PN,PND

double precision::itgr_b

x = dble(xr)

allocate(PMN(0:10,0:max(nu,n,10)),PMND(0:10,0:max(nu,n,10)))
allocate(PN(0:max(p,1)),PND(0:max(p,1)))

!Calculate associated Legendre polynominals
call lpmn(10,1,max(nu,n,10),X,PMN,PMND)
!call lpmn(100,1,max(nu,n),X,PMN,PMND)
P1nu = PMN(1,nu)
P1n = PMN(1,n)

call lpn(max(p,1),X,PN,PND)
!call lpn(p,X,PN,PND)
PDp  = PND(p)

deallocate(PMN,PMND,PN,PND)

itgr_b = P1nu*P1n*PDp

return
end function itgr_b


subroutine complex_leqs_solver(n,np,T,p,r)
!Input : n, T, p 
!output: r
implicit none
integer::i,j,n,np
complex,dimension(n,n)::T
complex,dimension(n)::p,r
!local variables
real,dimension(np,np)::a
real,dimension(np)::d,b
integer,dimension(np)::indx

!---------------------------------------
! Create 2n x 2n set of "real" equations
! where n = 2 * nstop
!---------------------------------------
a = 0.0
b = 0.0
r = 0.0

do i=1,n
        do j=1,n
                a(i,j) = real(T(i,j))
                a(i,j+n) = -aimag(T(i,j))
                a(i+n,j) = aimag(T(i,j))
                a(i+n,j+n) = real(T(i,j))
        enddo
        b(i) = real(p(i))
        b(n+i) = aimag(p(i))
enddo

!do i=1,np
        !do j=1,np
        !        write(*,*) i,j,a(i,j)
        !enddo
        !write(*,*) b(i)
!enddo
!stop

call ludcmp(a,np,np,indx,d)
call lubksb(a,np,np,indx,b)

! Outputing complex solution
do i=1,n
        r(i) = cmplx(b(i),b(n+i))
enddo

return
end subroutine complex_leqs_solver



subroutine integration_of_Sp(iqcor,D,p,k,x,Sp)
integer::n,nn,p
integer::iqcor,isol,ll
double precision::umin,umax,du,k
double precision,allocatable,dimension(:)::u,xx
complex(kind(0d0)),allocatable,dimension(:)::intg
double precision,allocatable,dimension(:)::intg_unit
!double precision,dimension(0:75)::SJ,SY
!double precision,dimension(0:100)::SJ,SY
double precision,allocatable,dimension(:)::SJ,SY
double precision::jp,yp
double precision::x
double precision::JPH
complex(kind(0d0))::SH,SPH
complex(kind(0d0))::wa,Sp
double precision::D,fc
double precision,parameter::pi=acos(-1.0d0)
double precision::floorvalue
double precision::a1,a2,a3,a4,b1,b2,b3,b4,lnxa,lnxb
double precision::unitary,error,error2
integer::j

!---------------------------
!old boundary (Jablonski+94)
!---------------------------
!a1=-3.6693021d-8
!a2=-3.1745158d-5
!a3=2.1567720d-2
!a4=9.9123380d-1
!
!b1=-9.9921351d-8
!b2=8.7822303d-6
!b3=1.0238752d-2
!b4=3.7588265
!
!---------------------------
! new boundary defined by RT
!---------------------------
a1=1.69496268177237e-08 
a2=-2.43299782942114e-05
a3=0.0158750501131321 
a4=1.00672154148706

b1=2.49355951047228e-08 
b2=-2.9387731648675e-05 
b3=0.0135005554796179 
b4=3.72312019844119
!

floorvalue = 1.0e-30

!
if(iqcor .eq. 4) then
        umax = 2.0d0 * x * sqrt(5.0/3.0)
        !write(*,*) x,k,umax
        !stop
        umin = umax/1.d7
        nn = 10000
else
        umax = 1.0d1 * max(1.0d0,x)
        umin = umax/1.d10
        nn = 10000 * max(1,int((x/10.0**(2.5))))
endif

allocate(u(1:nn),xx(1:nn),intg(1:nn),intg_unit(1:nn))
du = (umax/umin) ** (1.0d0/dble(nn-1))

do n=1,nn
      u(n) = umin * du ** (dble(n-1))
enddo

!-------------------------------------------
! Integration using Trapezoid formula
!-------------------------------------------
intg = cmplx(0.0d0,0.0d0,kind(0d0))
intg_unit = 0.0d0
!open(1,file="sbessel_test.dat",status="unknown")
do n=1,nn

        lnxa=a1*dble(p)**3.0+a2*dble(p)**2.0+a3*dble(p)+a4
        lnxb=b1*dble(p)**3.0+b2*dble(p)**2.0+b3*dble(p)+b4

        if(log(u(n)) .lt. lnxa) then
                isol = 1
        elseif(log(u(n)) .gt. lnxb) then
                isol = 3
        else
                isol = 2
        endif

        ll=p+100
        allocate(SJ(0:ll),SY(0:ll))
                call sphbessel(ll,u(n),SJ,SY,isol)
                jp=SJ(p)
                yp=SY(p)
               
                !check
                !do j=0,p
                !        write(1,*) u(n),j,SJ(j),SY(j)
                !enddo
                !write(1,*)
                !check
        deallocate(SJ,SY)


        if(jp*0.d0 /= 0.d0 .or. isnan(jp) .eqv. .true.) then
               write(*,*) "--------------------------------------------------------"
               write(*,*) "Error: NaN is found at the spherical Bessel function of the first kind"
               write(*,*) "argument x : ",u(n)
               write(*,*) "order    p : ",p
               write(*,*) "j_p(x)     : ",jp
               if(isol .eq. 1) then
                       write(*,*) "scheme     : the series expansion"
               elseif(isol .eq. 2) then
                       write(*,*) "scheme     : downward recurrence"
               elseif(isol .eq. 3) then
                       write(*,*) "scheme     : upward recurrence"
               endif

               write(*,*) "Simulation is aborted !"
               write(*,*) "--------------------------------------------------------"
               stop

        elseif(yp*0.d0 /= 0.d0 .or. isnan(yp) .eqv. .true.) then
               write(*,*) "--------------------------------------------------------"
               write(*,*) "Error: NaN is found at the spherical Bessel function of the second kind"
               write(*,*) "argument x : ",u(n)
               write(*,*) "order    p : ",p
               write(*,*) "y_p(x)     : ",yp 
               write(*,*) "Simulation is aborted !"
               write(*,*) "--------------------------------------------------------"
               stop

        endif

        JPH = sqrt(2.0d0*u(n)/pi)*jp 
        SH = cmplx(jp,yp,kind(0d0))
        SPH = sqrt(2.0d0*u(n)/pi)*SH

        if(iqcor .eq. 4) then
                intg(n) = u(n) * JPH * SPH * fc(iqcor,u(n),x,D)
                intg_unit(n) = u(n) ** 2.0d0 * fc(iqcor,u(n),x,D)
                !write(*,*) intg_unit(n)
                !if(u(n)*k/x .ge. 2.0d0) exit
        else        
                !intg(n) = u(n) ** (D-2.0d0) * JPH * SPH * exp(-c*(u(n)/x)**D)
                intg(n) = u(n) ** (D-2.0d0) * JPH * SPH * fc(iqcor,u(n),x,D)
                intg_unit(n) = u(n) ** (D-1.0d0) * fc(iqcor,u(n),x,D)
                if(fc(iqcor,u(n),x,D) .le. floorvalue) exit
        endif

        JPH = 0.0d0
        SPH = 0.0d0

enddo

wa = cmplx(0.0d0,0.0d0,kind(0d0))

unitary = 0.0d0
do n=1,nn-1
        wa = wa + 0.5d0 * (intg(n)+intg(n+1)) * (u(n+1)-u(n))
        unitary = unitary + 0.5d0 * (intg_unit(n)+intg_unit(n+1))*(u(n+1)-u(n))
enddo

if(iqcor .eq. 4) then
        unitary = unitary / x ** 3.0d0
else
        unitary = unitary / x ** D
endif

error = abs(1.0d0-unitary)
        
        !do n=1,nn-1
        !write(*,*) u(n),intg_unit(n)
        !enddo
        !stop

        !write(*,*) x,error
        if(error .ge. 1.0e-3) then
                write(*,*) "--------------------------------------------------------"
                write(*,*) "Check unitary condition : two-points correlation function"
                write(*,*) "Numerical integration of g(u) with error of ",error*1.d2," (%)"
                write(*,*) "which exceed 0.1%. This means sp(kRg) integration"
                write(*,*) "may not converge."
                write(*,*) "Simulation is aborted !"
                write(*,*) "--------------------------------------------------------"
                stop
        endif

if(iqcor .eq. 4) then
        Sp = pi/(4.0d0*x**3.0d0)*wa
else
        Sp = pi/(4.0d0*x**D)*wa
endif

deallocate(u,xx,intg,intg_unit)

        !-------------------------------------------
        ! Analytic solution for Im(S_{p=0})
        !-------------------------------------------
        !hg = 0.0d0
        !al = 1.0d0/2.0d0
        !bb = 1.5d0
        !xxx = -2.0d0*xx(ix)**2.0!/2.0
        !call chgm(al,bb,xxx,hg)
        !if(hg*0.d0 /= 0.d0) then
        !hg=(sqrt(pi)/2.0d0)*(sqrt(2.0d0)*xx(ix))**(-1.0d0)
        !endif
        !ImS0 = -(sqrt(2.0d0*pi)/(4.0d0*xx(ix)))*hg
        !-------------------------------------------
        ! Output
        !-------------------------------------------
        !write(*,*) xx(ix),ix,dble(Sp),-aimag(Sp),dble(SpGL),-aimag(SpGL),-ImS0
        !write(*,*) xx(ix),dble(Sp),-aimag(Sp),-ImS0

return
end subroutine integration_of_Sp


!-----------------------------------------
!
! Cut off models for correlation function
!
!-----------------------------------------
function fc(iqcor,u,x,df)
implicit none
double precision::u,x,df
integer::iqcor
!local
double precision::c,fc,xc

if(iqcor .eq. 1) then

        !RT Gaussian cut-off model  
        c = (df/4.0d0)
        fc = (2.0d0*c**(df/2.0)/gamma(df/2.0))*exp(-c*(u/x)**2.0) 

elseif(iqcor .eq. 2) then

        c = sqrt(df*(df+1.0d0)/2.0d0)
        fc = c**(df)/gamma(df) * exp(-c*(u/x))

elseif(iqcor .eq. 3) then

        !Botet et al. 1995 & 1997's cut-off model
        c = 0.5d0
        fc = c*df*exp(-c*(u/x)**df)

elseif(iqcor .eq. 4) then
        
        !kRg --> kRc
        xc = sqrt(5.0d0/3.0d0) * x
        
        !homogeneous sphere: two-points correlation function
        fc = (3.0d0/16.0d0)*(u/xc)**3.0-(9.0d0/4.0d0)*(u/xc)+3.0d0

        !1/4piRg^3 --> 1/4piRc^3
        fc = (3.0d0/5.0d0)**(1.5) * fc

endif

return
end function fc


!---------------------------------------------
!
! Calculate spherical bessel function 
!     using recurrence relation
!
!---------------------------------------------
!subroutine sphbessel(M,x,SJ,DSJ,SY,DSY)
subroutine sphbessel(M,x,SJ,SY,isol)
implicit none
integer::M,n
double precision::K0,K1,x,s
double precision,dimension(0:M)::K,SJ,SY
double precision::Y0,Y1
double precision::floorvalue,ceilingvalue
double precision,parameter::pi=acos(-1.0d0)
integer::i,isol
double precision::wa,xi,FN
double precision::eps

eps=1.0d-30
floorvalue = 1.0d-300
ceilingvalue = 1.0d250

!------------------------------------------------
! Downward recurrence for 1st spherical Bessel 
! Note: Upward recurrence is numerically unstable.
!-----------------------------------------------
!K1 = 0.0d0
!K0 = 1.0d0
!K=0.0d0
SJ=0.0d0
SJ(0) = dsin(x)/x

!write(*,*) x
!stop

!if(x .gt. 3.0d-1 .and. x .le. 50.0d0) then

        if(isol .eq. 1) then
        !--------------------
        !series expansion
        !--------------------
        
        SJ(0) = dsin(x)/x
        !SJ(1) = dsin(x)/(x*x)-dcos(x)/x
        !FN = x / 3.0d0
        FN = 1.0d0
        do n=1,M
        FN = x / dble(2*n+1) * FN
                xi = 1.0d0
                wa = 0.0d0
                do i=1,100
                        xi = - x * x / (dble(2*i)*(2.0d0*dble(i+n)+1)) * xi
                        if(abs(xi) .le. eps) exit
                        !        write(*,*) "convergence order : ",i        
                        !        exit
                        !endif
                        wa = wa + xi
                enddo
        SJ(n) = FN * (1.0d0 + wa) 
        if(DABS(SJ(n)) .le. floorvalue) exit
        enddo

        elseif(isol .eq. 2) then

        !--------------------
        !Downward recurrence
        !--------------------
        K1 = 0.0d0
        K0 = 1.0d0
        K=0.0d0

        SJ=0.0d0
        SJ(0) = dsin(x)/x
        do n=M,0,-1
                K(n) = -K1 + ((2.0d0*dble(n+1)+1)/x)*K0
                K1 = K0
                K0 = K(n)
        enddo

        s = SJ(0)/K(0)

        do n=1,M
                SJ(n)=s*K(n)
                !if(DABS(SJ(n)) .le. floorvalue) SJ(n) = 0.0d0
                if(DABS(SJ(n)) .le. floorvalue) exit
                !write(*,*) n,K(n),SJ(n)
        enddo

        elseif(isol .eq. 3) then

        !--------------------
        !Upward recurrence
        !--------------------
        SJ(0) = dsin(x)/x
        SJ(1) = dsin(x)/(x*x)-dcos(x)/x
        K0 = SJ(0)
        K1 = SJ(1)
        do n=1,M-1
         SJ(n+1) = ((2.0d0*dble(n)+1)/x)*K1 - K0
         K0 = K1
         K1 = SJ(n+1)
        enddo


        endif

!elseif(x .le. 3.0d-1) then
!
!        !write(*,*) "Small x limit; x = ",x
!        SJ(0) = 1.d0
!        do i=1,M
!                SJ(i) = SJ(i-1) * x / dble(2*i+1)
!
!                if(SJ(i) .le. floorvalue) exit
!                !then
!                !        SJ(i) = 0.0d0
!                !endif
!        enddo
!        !stop
!
!elseif(x .gt. 50.0d0) then
!!        !
!!        !For large argument, I use asymptotic relation.
!!        !
!        !write(*,*) "Large"
!        do n=1,M
!                SJ(n)=dcos(x-dble(n+1)*pi/2.0d0)/x
!                if(DABS(SJ(n)) .le. floorvalue) exit
!        enddo
!
!endif


!--------------------------------------------
! Upward recurrence for 2nd spherical Bessel
!--------------------------------------------
SY=0.0d0
SY(0) = -dcos(x)/x
SY(1) = -dcos(x)/(x*x) - dsin(x)/x 

Y0 = SY(0)
Y1 = SY(1)

!if(x .lt. 100.0d0) then
        do n=2,M
                SY(n) = ((2.0d0*dble(n-1)+1.0d0)/x)*Y1-Y0
                Y0 = Y1
                Y1 = SY(n)
                
                !write(*,*) n,M,SY(n) 
                if(DABS(SY(n)) .ge. ceilingvalue) exit
                !SY(n) = 0.0d0
        enddo
!elseif(x .ge. 100.0d0) then
!        do n=1,M
!                SY(n)=dsin(x-dble(n+1)*pi/2.0d0)/x
!        enddo
!endif

!!----------------------------------------------
!! Upward recurrence for derivatives of 1st and 2nd
!! spherical Bessel function
!!----------------------------------------------
!DSJ(0) = (dcos(x)-dsin(x)/x)/x
!DSY(0) = (dsin(x)+dcos(x)/x)/x
!do n=1,M
!        DSJ(n) = dble(n)/dble(2*n+1)*SJ(n-1)-dble(n+1)/dble(2*n+1)*SJ(n+1)
!        DSY(n) = dble(n)/dble(2*n+1)*SY(n-1)-dble(n+1)/dble(2*n+1)*SY(n+1)
!
!        if(DABS(DSJ(n)) .le. floorvalue) DSJ(n) = 0.0d0
!        !if(DABS(DSY(n)) .ge. ceilingvalue) DSY(n) = 0.0d0
!enddo

return
end subroutine sphbessel
!******************************************************************
!*      Purpose: This program computes the confluent              *
!*               hypergeometric function M(a,b,x) using           *
!*               subroutine CHGM                                  *
!*      Input  : a  --- Parameter                                 *
!*               b  --- Parameter ( b <> 0,-1,-2,... )            *
!*               x  --- Argument                                  *
!*      Output:  HG --- M(a,b,x)                                  *
!*      Example:                                                  *
!*                  a       b       x          M(a,b,x)           *
!*                -----------------------------------------       *
!*                 1.5     2.0    20.0     .1208527185D+09        *
!*                 4.5     2.0    20.0     .1103561117D+12        *
!*                -1.5     2.0    20.0     .1004836854D+05        *
!*                -4.5     2.0    20.0    -.3936045244D+03        *
!*                 1.5     2.0    50.0     .8231906643D+21        *
!*                 4.5     2.0    50.0     .9310512715D+25        *
!*                -1.5     2.0    50.0     .2998660728D+16        *
!*                -4.5     2.0    50.0    -.1807604439D+13        *
!* -------------------------------------------------------------- *
!* REFERENCE: "Fortran Routines for Computation of Special        *
!*             Functions jin.ece.uiuc.edu/routines/routines.html" *
!*                                                                *
!*                              F90 Release By J-P Moreau, Paris. *
!*                                     (www.jpmoreau.fr)          *
!******************************************************************
        SUBROUTINE CHGM(A,B,X,HG)
!       ===================================================
!       Purpose: Compute confluent hypergeometric function
!                M(a,b,x)
!       Input  : a  --- Parameter
!                b  --- Parameter ( b <> 0,-1,-2,... )
!                x  --- Argument
!       Output:  HG --- M(a,b,x)
!       Routine called: GAMMA for computing ג(x)
!       ===================================================
        IMPLICIT DOUBLE PRECISION (A-H,O-Z)
        PI=3.141592653589793D0
        A0=A
        A1=A
        X0=X
        HG=0.0D0
        IF (B.EQ.0.0D0.OR.B.EQ.-ABS(INT(B))) THEN
           HG=1.0D+300
        ELSE IF (A.EQ.0.0D0.OR.X.EQ.0.0D0) THEN
           HG=1.0D0
        ELSE IF (A.EQ.-1.0D0) THEN
           HG=1.0D0-X/B
        ELSE IF (A.EQ.B) THEN
           HG=DEXP(X)
        ELSE IF (A-B.EQ.1.0D0) THEN
           HG=(1.0D0+X/B)*DEXP(X)
        ELSE IF (A.EQ.1.0D0.AND.B.EQ.2.0D0) THEN
           HG=(DEXP(X)-1.0D0)/X
        ELSE IF (A.EQ.INT(A).AND.A.LT.0.0D0) THEN
           M=INT(-A)
           R=1.0D0
           HG=1.0D0
           DO 10 K=1,M
              R=R*(A+K-1.0D0)/K/(B+K-1.0D0)*X
10            HG=HG+R
        ENDIF
        IF (HG.NE.0.0D0) RETURN
        IF (X.LT.0.0D0) THEN
           A=B-A
           A0=A
           X=DABS(X)
        ENDIF
        IF (A.LT.2.0D0) NL=0
        IF (A.GE.2.0D0) THEN
           NL=1
           LA=INT(A)
           A=A-LA-1.0D0
        ENDIF
        DO 30 N=0,NL
           IF (A0.GE.2.0D0) A=A+1.0D0
           IF (X.LE.30.0D0+DABS(B).OR.A.LT.0.0D0) THEN
              HG=1.0D0
              RG=1.0D0
              DO 15 J=1,500
                 RG=RG*(A+J-1.0D0)/(J*(B+J-1.0D0))*X
                 HG=HG+RG
                 IF (DABS(RG/HG).LT.1.0D-15) GO TO 25
15            CONTINUE
           ELSE
              CALL GAMMA(A0,TA)
              CALL GAMMA(B,TB)
              XG=B-A
              CALL GAMMA(XG,TBA)
              SUM1=1.0D0
              SUM2=1.0D0
              R1=1.0D0
              R2=1.0D0
              DO 20 I=1,8

                 R1=-R1*(A+I-1.0D0)*(A-B+I)/(X*I)
                 R2=-R2*(B-A+I-1.0D0)*(A-I)/(X*I)

                 !S: RT 
!                 R1=-R1*(A+I-1.0D0)*(B-A+I-1.0D0)/(X*I)
!                 R2=-R2*(B-A+I-1.0D0)*(A-I)/(X*I)
                 !E: RT

                 SUM1=SUM1+R1
20               SUM2=SUM2+R2
              HG1=TB/TBA*(X)**(-A)*DCOS(PI*A)*SUM1
              HG2=TB/TA*DEXP(X)*X**(A-B)*SUM2

              !S: RT
!              HG1=(TB/TBA)*(X)**(-A)!*SUM1
!              HG2=TB/TA*DEXP(-X)*(-X)**(A-B)*SUM2
!              if(HG2 .le. 1.0e-30) then
!              HG2 = 0.0d0
!              end if
!11               write(*,*) HG1,HG2 
              !E: RT

               HG=HG1+HG2
           ENDIF
25         IF (N.EQ.0) Y0=HG
           IF (N.EQ.1) Y1=HG
30      CONTINUE
        IF (A0.GE.2.0D0) THEN
           DO 35 I=1,LA-1
              HG=((2.0D0*A-B+X)*Y1+(B-A)*Y0)/A
              Y0=Y1
              Y1=HG
35            A=A+1.0D0
        ENDIF

        IF (X0.LT.0.0D0) HG=HG*DEXP(X0)
        A=A1
        X=X0
        RETURN
        END SUBROUTINE CHGM


        SUBROUTINE GAMMA(X,GA)
!       ==================================================
!       Purpose: Compute gamma function ג(x)
!       Input :  x  --- Argument of ג(x)
!                       ( x is not equal to 0,-1,-2,תתת)
!       Output:  GA --- ג(x)
!       ==================================================
        IMPLICIT DOUBLE PRECISION (A-H,O-Z)
        DIMENSION G(26)
        PI=3.141592653589793D0
        IF (X.EQ.INT(X)) THEN
           IF (X.GT.0.0D0) THEN
              GA=1.0D0
              M1=X-1
              DO 10 K=2,M1
10               GA=GA*K
           ELSE
              GA=1.0D+300
           ENDIF
        ELSE
           IF (DABS(X).GT.1.0D0) THEN
              Z=DABS(X)
              M=INT(Z)
              R=1.0D0
              DO 15 K=1,M
15               R=R*(Z-K)
              Z=Z-M
           ELSE
              Z=X
           ENDIF
           DATA G/1.0D0,0.5772156649015329D0,  &
                -0.6558780715202538D0, -0.420026350340952D-1, &
                0.1665386113822915D0,-.421977345555443D-1,    &
                -.96219715278770D-2, .72189432466630D-2,      &
                -.11651675918591D-2, -.2152416741149D-3,      &
                .1280502823882D-3, -.201348547807D-4,         &
                -.12504934821D-5, .11330272320D-5,            &
                -.2056338417D-6, .61160950D-8,                &
                .50020075D-8, -.11812746D-8,                  &
                .1043427D-9, .77823D-11,                      & 
                -.36968D-11, .51D-12,                         & 
                -.206D-13, -.54D-14, .14D-14, .1D-15/
           GR=G(26)
           DO 20 K=25,1,-1
20            GR=GR*Z+G(K)
           GA=1.0D0/(GR*Z)
           IF (DABS(X).GT.1.0D0) THEN
              GA=GA*R
              IF (X.LT.0.0D0) GA=-PI/(X*GA*DSIN(PI*X))
           ENDIF
        ENDIF
        RETURN
        END SUBROUTINE GAMMA

!end of file mchgm.f90
subroutine spout(iqcor,df)
implicit none
integer::iqcor,p,imax,i
double precision::df,k,xg
double precision::xmin,xmax,dx
complex(kind(0d0))::sp

k = 1.0d0
!xmin=1.0d-2
xmin=1.0d-2
xmax=1.0d4
imax=250
!imax=10
dx=(xmax/xmin)**(1.0d0/dble(imax-1))

write(*,*) "Writing Sp ..."
open(21,file="sp.out",status="unknown")
write(21,6600) "exp. order p","k*R_g","Re(sp)","Im(sp)","Im(s_{p=0})"
!do p=0,0
do p=0,25
        write(*,*) "p = ",p
        do i=1,imax
        write(*,*) i," / ",imax," x = ",xg
        xg = xmin * dx ** (dble(i-1))
        call integration_of_Sp(iqcor,df,p,k,xg,Sp) 
        write(21,6000) p,xg,dble(Sp),aimag(Sp),-(acos(-1.0d0)/8.0d0)*xg**(-2.0d0)
        enddo

write(21,*)
enddo
close(21)
write(*,*) "sp.out is produced"

6000 format(' ',I15,1P4E15.5)
6600 format(' ',4A15)

stop

return
end subroutine spout
