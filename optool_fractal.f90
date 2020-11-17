!--------------------------------------------------------------------------------
!
!       OPACFRACTAL version 3.0
!
!--------------------------------------------------------------------------------
!
! OPACFRACTAL v3.0 computes light scattering properties of randomly oriented 
! fractal dust aggregates by means of the modified mean field theory developed 
! in Tazaki & Tanaka (2018). This code is also capable of computing the light 
! scattering solution based on the Rayleigh-Gans-Debye theory and the Mean field 
! theory.
! 
! If you wish to publish a paper that contains results obtained by this code, 
! please acknowledge this code and cite relevant papers. 
!
! Disclaimer:
! I reject all responsibility for the use of this code. Although it has been 
! tested, the code may still contain a bug. I am not responsible for any damages
! caused by the use of the code. If you find a bug, please let me know.
!
!                                                   Ryo Tazaki (2020/Nov/11)                
!                                                   email: r.tazaki -at- uva.nl 
!
!--------------------------------------------------------------------------------
!
! LIMITATION OF THE CODE
!
!--------------------------------------------------------------------------------
!
! Here, I summarize some limitations of the code for each approximation (iqsca).
! The limitation can be simply judged by the phase shift induced by an aggregate
! (Equation 9 in Tazaki & Tanaka 2018; see also Section 3.2 in this paper).
!
! iqsca=1: All outputs would be physically reasonable for the phase shift <~ 1. 
! 
! iqsca=2: The extinction cross section would be calculaed without limitation.
!          However, the other outputs would be reliable for the phase shift < ~1.
!
! iqsca=3: The extinction cross section would be calculaed without limitation.
!          Scattering and absorption cross sections could be calculated
!          for the phase shift > 1, however too large phase shift may
!          cause some problem.
!          The asymmetry parameter and the sattering matrix elements would be
!          reliable for the phase shift < ~1.
!
! For the two-point correlation function, I suggest to use iqcor=1 because
! it is numerically stable and is likely to reproduce optical properties
! of fractal aggregates.
!
!--------------------------------------------------------------------------------
!
! IMPORTANT REFERENCES
!
!--------------------------------------------------------------------------------
!
! If you wish to publish a paper by using this code, following references should 
! be cited depending on the code options (iqsca and iqcor) you used.
!
! Light scattering solver:
! iqsca = 1     : Tazaki et al. 2016, ApJ, 823, 70
! iqsca = 2     : Botet et al. 1997, ApOpt, 36, 8791
! iqsca = 3     : Tazaki & Tanaka 2018, ApJ, 860, 79
!
! The two-point correlation function of monomers:
! iqcor = 1     : Tazaki et al. 2016, ApJ, 823, 70 
! iqcor = 2     : Berry & Percival 1986, AcOpt, 33, 577
! iqcor = 3     : Botet et al. 1995, JPhA, 28, 297
!
! The geometric cross section of aggregates (to be used when iqsca=3):
! iqgeo = 2     : Okuzumi et al. 2009, ApJ, 707, 1247
! iqgeo = 3     : Empirical formulae obtained by mean-field extinction
!
!--------------------------------------------------------------------------------
!
! LIST OF SUBROUTINES BY OTHER AUTHORS
!
!--------------------------------------------------------------------------------
!
! OPACFRACTAL v3.0 contains following subroutines provided by other authors.
! I thank the authors for the availability of the subroutines.
!
! lorenz_mie  : BHMIE written by Bruce T. Draine with modifications 
! renorm_mie  : BHMIE written by Bruce T. Draine with modifications
! gamma       : Zhang, S. and Jin, J. (1996) "Computation of Special Functions." 
! chgm        : Zhang, S. and Jin, J. (1996) "Computation of Special Functions." 
! lpmns       : Zhang, S. and Jin, J. (1996) "Computation of Special Functions."
! lpn         : Zhang, S. and Jin, J. (1996) "Computation of Special Functions."
! ludcmp      : Press, W. H. et al. (1997), "Numerical Recipes in Fortran 77"
! lubksb      : Press, W. H. et al. (1997), "Numerical Recipes in Fortran 77"
! gauleg      : Press, W. H. et al. (1997), "Numerical Recipes in Fortran 77"
!
!--------------------------------------------------------------------------------
!
! INPUT PARAMETERS
!
!--------------------------------------------------------------------------------
!
! lmd           : wavelength     (micron)
! R0            : monomer radius (micron)
! Df            : fractal dimension 
! PN            : number of monomers
! k0            : fractal prefactor; PN = k0 (Rg/R0)^Df
! refrel        : complex refractive index of a monomer particle
! nang          : Number of angular grid between 0 ang 90 degrees.
!                 If nang=91, then angular width is 1 degree.
!                 <TIPS>
!                 For iqsca=3, I suggest the users to check the convergence
!                 of the results by changing nang.
! iqsca         : Switch for the light scattering solver
!                 iqsca = 1  : Rayleigh-Gans-Debye theory
!                 iqsca = 2  : Mean-field theory 
!                 iqsca = 3  : Modified mean-field theory
! iqcor         : Switch for the two-point correlation function 
!                 iqcor = 1  : Gaussian cutoff 
!                 iqcor = 2  : Exponential cutoff 
!                 iqcor = 3  : Fractal dimension cutoff
! iqgeo         : Switch for the two-point correlation function 
!                 iqgeo = 1  : pi * rc^2
!                 iqgeo = 2  : Okuzumi et al. (2009)
!                 iqgeo = 3  : Empirical formulae from mean-field extinction
! iquiet        : Switch for standard output during a calculation
!                 iquiet = 0 : show stdout 
!                 iquiet = 1 : suppress stdout (including warnings)
!
!--------------------------------------------------------------------------------
!
! OUTPUT PARAMETERS
!
!--------------------------------------------------------------------------------
!
! Cext               : extinction cross section (micron^2)
! Csca               : scattering cross section (micron^2)
! Cabsp              : absorption cross section (micron^2)
! Gsca               : asymmetry parameter
! Smat(4,2*nang-1)   : Scattering matrix elements S_ij
!                      where S_ij obeys Bohren & Huffman's definition.
!                      The first column contains each element of scat. matrix
!                        Smat(1,:) : S_11,               
!                        Smat(2,:) : S_12,               
!                        Smat(3,:) : S_33,              
!                        Smat(4,:) : S_34,                
!                      and the second column contains values at each scattering
!                      angle from 0 deg to 180 deg with 2*nang-1 grid points.
! dphi               : Phase shift induced by the aggregates
!
!--------------------------------------------------------------------------------
!
! HISTORY
!
!--------------------------------------------------------------------------------
!
! version 1.0 (July 24 2017)
! - Initial release 
!
! version 2.2 (Oct. 24 2017)
! - Implemented geometrical optics approximation (Modified Mean Field; MMF)
! - Implemented a stable algorithm for computations of the spherical Bessel 
!   function of the first kind at large order and/or large argument.
! - Fixed bug of Legendre polynomial calculations with large order.
! - Implemented dynamic allocation for some special function calculations.
! - Fixed bug about compiler dependence (ifort or gfortran)
!
! version 2.3 (Dec. 05 2017)
! - Implemented efficient computation of $a(\nu,n,p)$ and $b(\nu,n,p)$ 
!   using the parity arguments.
!
! version 2.5 (May. 31 2020)
! - Fixed a bug about initialization of Qsca and gsca in (renormalized) 
!   Lorentz-Mie routine. This does not change the results, since both of 
!   them are not used to obtain final results.
! 
! version 3.0 (Nov. 05, 2020) 
! - BUG fixed: Add a simple prescription for negative S(q) for iqcor=3.
! - Fixed a code crush for OpenMP mode.
! - Fixed some over/underflow signals.
! - Fixed a type mismatch of "d" for the subroutine 'ludcmp'.
! - Fixed a type mismatch of cn1, cn2 in the main subroutine.
! - Unified notation for Cabs in iqsca=1 for safety.
! - Updated f77 subroutines to f90.
! - Implemented parameterized real precision
! - Implemented more efficient loop structure at calculations of 
!   gaunt coefficients : a(nu,n,p) and b(nu,n,p).
! - Removed the iqcor=4 option because it is practically unimportant.
! - Implemented Okuzumi et al. (2009) formula for geometrical cross-section.
! - Added a return variable "phase shift"
! - Added a new switch for geometric cross section
!
!--------------------------------------------------------------------------------
!       Some constants
!--------------------------------------------------------------------------------
module types
implicit none
integer, parameter      :: dp = selected_real_kind(P=15)
end module types

module const
use types
implicit none
real(kind=dp),parameter :: pi = 3.1415926535897932384_dp
real(kind=dp),parameter :: twopi  = 2.0_dp * pi    
real(kind=dp),parameter :: halfpi = 0.5_dp * pi    
real(kind=dp),parameter :: d2r    = pi/180.0_dp  ! degree to radian
real(kind=dp),parameter :: r2d    = 1.0_dp/d2r   ! radian to degree
integer      ,parameter :: jm     = 400          ! Number of grid points of Gauss-Ledengre
                                                 ! for integration of Gaunt coefficients.
integer,parameter       :: nn     = 10000        ! Number of grid points in the numerical
                                                 ! integration for Sp(kRg)
end module const

!--------------------------------------------------------------------------------
!
!
!       MAIN ROUTINE
!
!
!--------------------------------------------------------------------------------

subroutine meanscatt(lmd,R0,PN,Df,k0,refrel,iqsca,iqcor,iqgeo,iquiet,nang,&
                Cext,Csca,Cabsp,g_asym,Smat,dphi)
use types; use const
implicit none

!--------------------------------------------------------------------------------
!       Input variables
!--------------------------------------------------------------------------------

integer                 :: iqsca              ! switch of method
integer                 :: iqcor              ! switch of correlation function
integer                 :: iqgeo              ! switch of geometric cross section
integer                 :: iquiet             ! switch of suppress warning messages
integer                 :: nang               ! Number of angle mesh between 0 to 90 deg.
real(kind=dp)           :: PN                 ! Number of monomers
real(kind=dp)           :: R0                 ! Monomer radius
real(kind=dp)           :: df, k0             ! fractal dimension and prefactor
real(kind=dp)           :: lmd                ! wavelength
complex(kind=dp)        :: refrel             ! complex refractive index

!--------------------------------------------------------------------------------
!       Output variables
!--------------------------------------------------------------------------------

real(kind=dp)           :: Cext               ! Extinction cross section
real(kind=dp)           :: Csca               ! Scattering
real(kind=dp)           :: Cabsp              ! Absorption
real(kind=dp)           :: g_asym             ! asymmetry parameter
real(kind=dp)           :: dphi               ! Phase shift induced by aggregate
real(kind=dp)           :: ang(1:2*nang-1)    ! angle grids (from 0.0 to 180.0 deg)
real(kind=dp)           :: Smat(1:4,2*nang-1) ! Scatteing matrix elements
real(kind=dp)           :: PF(1:2*nang-1)     ! Phase function (BH def.)

!--------------------------------------------------------------------------------
!       Local variables
!--------------------------------------------------------------------------------
real(kind=dp),parameter :: eta=25.0_dp         
real(kind=dp),parameter :: qRg_crit=26.0_dp 
real(kind=dp)           :: Rg                 ! radius of gyration
real(kind=dp)           :: Rc                 ! characteristic radius
real(kind=dp)           :: k                  ! wavenumber
real(kind=dp)           :: x0,xg              ! size parameter of monomer and agg.
real(kind=dp)           :: xstop              ! used in Mie truncation 
real(kind=dp)           :: xi                 ! used in xi (Berry corr)
integer::nstop,numax,nmax
integer::nu,n,p,j,pmin,pmax
integer::degmax,order
real(kind=dp),allocatable,dimension(:)::x,w
real(kind=dp)::x1,x2,anunp,bnunp
real(kind=dp)::dang,S11,S12,S33,S34
real(kind=dp)::q,Sq,al,bb,xx,cc 
real(kind=dp)::nrm
real(kind=dp)::GC,tau
real(kind=dp),allocatable,dimension(:)::PMN,PMND,LP,DLP
real(kind=dp),allocatable,dimension(:,:)::AL1N,LN,DLN
real(kind=dp)::ff,dphic,dphi0
complex(kind=dp)::sumA,sumB,Sp_tmp
complex(kind=dp),allocatable,dimension(:)::an,bn,y,r,Sp
complex(kind=dp),allocatable,dimension(:,:)::ad,dd,S
complex(kind=dp),allocatable,dimension(:,:,:)::T
complex(kind=dp),dimension(1:2*nang-1)::S1,S2
complex(kind=dp)::cn1,cn2
complex(kind=dp):: mgmref

!--------------------------------------------------------------------------------
k       =  twopi/lmd                    ! wavenumber
Rg      =  R0 * (PN/k0) ** (1.0_dp/df)  ! Radius of gyration of the aggregate
Rc      =  sqrt(5.0_dp/3.0_dp)*Rg       ! Characteristic radius of the aggregate
xg      =  k * Rg                       ! size parameter of aggregates
x0      =  k * R0                       ! size parameter of a monomer


!--------------------------------------------------------------------------------
!
! Truncation order of the monomer's scattering field 
!
!--------------------------------------------------------------------------------
xstop = x0 + 4.0_dp * x0 ** (1.0_dp/3.0_dp) + 2.0_dp
nstop = nint(xstop)
numax = nstop
nmax  = nstop

!--------------------------------------------------------------------------------
!
! Check input parameters
! 
!--------------------------------------------------------------------------------
if(iquiet .eq. 0) then
write(*,fmt='(1PE15.5,A)') lmd,   " =  wavelength (um)"
write(*,fmt='(1PE15.5,A)') R0 ,   " =  monomer radius (um)"
write(*,fmt='(1PE15.5,A)') Rg ,   " =  gyration radius (um)"
write(*,fmt='(1PE15.5,A)') Rc ,   " =  characteristic radius (um)"
write(*,fmt='(1PE15.5,A)') x0,    " =  size parameter of a monomer x0 (um)"
write(*,fmt='(1PE15.5,A)') xg,    " =  size parameter of an aggregate xg (um)"
write(*,fmt='(I15,A)'    ) nstop, " =  expansion order of scat. field."
endif

if (iqsca .ne. 1 .and. iqsca .ne. 2 .and. iqsca .ne. 3) then
        print *, 'ERROR: Inappropriate iqsca value.'
        print *, '       STOP.'
        stop
endif

if (iqcor .ne. 1 .and. iqcor .ne. 2 .and. iqcor .ne. 3) then
        print *, 'ERROR: Inappropriate iqcor value.'
        print *, '       STOP.'
        stop
endif

if (iqgeo .ne. 1 .and. iqgeo .ne. 2 .and. iqgeo .ne. 3) then
        print *, 'ERROR: Inappropriate iqgeo value.'
        print *, '       STOP.'
        stop
endif

if (nang .le. 1) then
        print *, 'ERROR: Inappropriate nang value.    '
        print *, '       nang should be larger than 1.'
        print *, '       STOP.'
        stop
endif

if (iquiet .ne. 0 .and. iquiet .ne. 1) then
        print *, 'ERROR: Inappropriate iquiet value.'
        print *, '       STOP.'
        stop
endif

if (PN .lt. 1.0) then
        print *, 'ERROR: Number of monomer is less than unity'
        print *, '       STOP.'
        stop
endif

if (df .gt. 3.0) then
        print *, 'ERROR: Fractal dimension should be less than 3.0.'
        print *, '       STOP.'
        stop
endif

if(numax + nmax .ge. 500 .and. iquiet .eq. 0) then
        print *, "WARNING: the truncation order of monomer's scattered light"
        print *, "         field exceeds the maximum value (=500).          "
        print *, "         This may cause a code crush at computations of   "
        print *, "         the spherical Bessel function of the 1st kind.   "
endif

!
! compute the phase shift (Equation 9 in Tazaki & Tanaka 2018)
!
ff = PN*(R0/RC)**3.0_dp
call MGmixing(refrel,mgmref,ff)
dphic = 2.0_dp*k*Rc*abs(mgmref-1.0_dp)
dphi0 = 2.0_dp*x0*abs(refrel-1.0_dp)
dphi  = max(dphic,dphi0)

if (dphi .ge. 1.0_dp .and. iquiet .eq. 0) then
        print *, 'WARNING: The phase shift by an aggregate exceeds unity.'
        print *, '         Output of scattering matrix elements are not  ' 
        print *, '         physically guaranteed.'
endif


!--------------------------------------------------------------------------------
!
! Computation part
!
!--------------------------------------------------------------------------------
!
! First, compute Lorenz-Mie scattering coefficients of a monomer
!
allocate(an(nstop),bn(nstop),ad(2,nstop),dd(2,nstop))

! initilizing arrays
ad = cmplx(0.0_dp,0.0_dp,kind=dp)
dd = cmplx(0.0_dp,0.0_dp,kind=dp)
!
call lorenz_mie(x0,refrel,nstop,an,bn)
!
ad(1,:) = an(:)
ad(2,:) = bn(:)
!
! Solve multiple scattering? 
!        iqsca = 1 : No
!        iqsca = 2 : Yes
!
if(iqsca .eq. 1) then
        !
        ! Scattering coefficients of the monomer is 
        ! equal to Lorenz-Mie solution.
        !
        dd(1,:)= ad(1,:)
        dd(2,:)= ad(2,:)

elseif(iqsca .ge. 2) then

        !--------------------------------------------------------------------------------
        !
        !  Calculate a(nu,n,p) and b(nu,n,p):
        !
        !               2p + 1  /+1
        !  a(nu,n,p) = -------- | dx P_nu^1(x) * P_n^1(x) * P_p(x),
        !                 2     /-1
        !
        !               2p + 1  /+1                         dP_p(x)
        !  b(nu,n,p) = -------- | dx P_nu^1(x) * P_n^1(x) * -------,
        !                 2     /-1                            dx
        !
        !  where P_n^m is the associated Legendre function (n: degree, m: order), 
        !  P_n is the Legendre polynominal function 
        !  (see Equations (29, 30) in Tazaki & Tanaka 2018).
        !  The integration is performed with the Gauss-Legendre quadrature.
        !
        !--------------------------------------------------------------------------------
        !
        ! Preparing for gauss-legendre quadrature
        !
        allocate(x(jm),w(jm))

        ! initializing arrays
        x  =   0.0_dp
        w  =   0.0_dp

        ! range of integration 
        x1 =  -1.0_dp
        x2 =   1.0_dp
        call gauleg(x1,x2,x,w,jm) 

        !--------------------------------------------------------------------------------
        !
        ! Storing values of the associated Legendgre function and 
        ! the Legendre polynominals and its derivative  at each Gauss-Ledengre point.
        !
        ! The values are stored in 
        !       AL1N(n,j) = P_n^1(x_j) : Associated Legendre function with order 1.
        !       LN  (n,j) = P_n(x_j)   : Ledengre polynominals
        !       DLN (n,j) = P_n(x_j)'  : Derivative of Legendre polynominals
        !
        !--------------------------------------------------------------------------------
        order  = 1
        degmax = nstop
        pmax   = 2*nstop
        allocate(PMN(0:degmax),PMND(0:degmax),LP(0:pmax),DLP(0:pmax))
        allocate(AL1N(0:degmax,1:jm),LN(0:pmax,1:jm),DLN(0:pmax,1:jm))
        do j=1,jm
                call lpmns(order,degmax,x(j),PMN,PMND)
                call lpn(pmax,x(j),LP,DLP)
                AL1N(:,j) = PMN(:)
                LN (:,j)  = LP(:)
                DLN(:,j)  = DLP(:)
        enddo
        deallocate(PMN,PMND,LP,DLP)
        
        !
        ! Preparing arrays for the structure integration Sp(kRg),
        ! and the translation matrix.
        !
        allocate(Sp(0:numax+nmax),T(2,numax,nmax))

        ! initializing arrays
        Sp = cmplx(0.0_dp,0.0_dp,kind=dp)
        T  = cmplx(0.0_dp,0.0_dp,kind=dp)
       
        !
        ! The structure integration Sp 
        !
                !call spout(iqcor,df) !for output 
        do p=0,numax+nmax
              Sp_tmp = cmplx(0.0_dp,0.0_dp,kind=dp)
              call integration_of_Sp(iqcor,df,p,xg,Sp_tmp)
              Sp(p) = Sp_tmp
        enddo
        
        !--------------------------------------------------------------------------------
        !
        ! Caltulate translation matrix coefficients
        !
        ! Translation coefficients: A and B (Equation # from Tazaki & Tanaka 2018)
        ! T(1,nu,n) : \bar{A}_{1,n}^{1,nu} defined in Equation (14) 
        ! T(2,nu,n) : \bar{A}_{1,n}^{1,nu} defined in Equation (15)
        ! anunp     : a(nu,n,p)            defined in Equation (29)
        ! bnunp     : b(nu,n,p)            defined in Equation (30)
        !
        ! Note about the rule of rum with respect to the index p
        ! (RT and Botet, in private communication, 2017).
        !
        ! Just before Equation (4) in Botet et al. (1997), it is written that 
        ! p should have the same parity as n+nu. However, this is only true
        ! true for a(nu,n,p), and not for b(nu,n,p) 
        ! Thus, in this code, index p runs ALL integers between |nu-n| and nu+n.
        !
        ! Tips for efficient computation. 
        ! The parity of the integrand of a and b leads following properties:
        ! - a(nu,n,p) .ne. 0 and b(nu,n,p) = 0 when p has the same parity as n+nu 
        ! - b(nu,n,p) .ne. 0 and a(nu,n,p) = 0 when p has the opposite parity to n+nu 
        !
        !--------------------------------------------------------------------------------
        do nu=1,numax
                do n=1,nmax
                      pmin=abs(n-nu)
                      pmax=nu+n
                      sumA = 0.0_dp
                      sumB = 0.0_dp
                      do p=pmin,pmax 
                                anunp=0.0_dp
                                bnunp=0.0_dp
                                if(mod(pmax,2) .eq. mod(p,2)) then
                                        do j=1,jm
                                            anunp = anunp + w(j) * AL1N(nu,j) * AL1N(n,j) * LN(p,j) 
                                        enddo
                                else
                                        do j=1,jm
                                            bnunp = bnunp + w(j) * AL1N(nu,j) * AL1N(n,j) * DLN(p,j) 
                                        enddo
                                endif
                                anunp = 0.5_dp * real(2*p+1,kind=dp) * anunp
                                bnunp = 0.5_dp * real(2*p+1,kind=dp) * bnunp
                                sumA = sumA + real(n*(n+1)+nu*(nu+1)-p*(p+1),kind=dp) * anunp * Sp(p) 
                                sumB = sumB + bnunp * Sp(p)
                      enddo
                      T(1,nu,n) = sumA * real(2*nu+1,kind=dp) / real(n*(n+1)*nu*(nu+1),kind=dp)
                      T(2,nu,n) = sumB * 2.0_dp * real(2*nu+1,kind=dp) / real(n*(n+1)*nu*(nu+1),kind=dp)
                enddo
        enddo
        deallocate(x,w,Sp,AL1N,LN,DLN)
        allocate(S(2*numax,2*nmax),y(2*nmax),r(2*nmax))
        
        ! initializing arrays
        S = cmplx(0.0_dp,0.0_dp,kind=dp)
        y = cmplx(0.0_dp,0.0_dp,kind=dp)
        r = cmplx(0.0_dp,0.0_dp,kind=dp)
        
        do n=1,nmax
                do nu=1,numax
                S(2*n-1,2*nu-1) = (PN-1.0_dp) * ad(1,n) * T(1,nu,n)  !a*A
                S(2*n-1,2*nu)   = (PN-1.0_dp) * ad(1,n) * T(2,nu,n)  !a*B
                S(2*n,2*nu-1)   = (PN-1.0_dp) * ad(2,n) * T(2,nu,n)  !b*B
                S(2*n,2*nu)     = (PN-1.0_dp) * ad(2,n) * T(1,nu,n)  !b*A
                enddo
                y(2*n-1) = ad(1,n)
                y(2*n)   = ad(2,n)
        enddo
        do n=1,2*nmax
               S(n,n) = 1.0_dp + S(n,n)
        enddo
        !
        ! Inverting the translation matrix elements (S), 
        ! and find the mean-field vector (r).
        !       S * r = y 
        !       ==> r = S^{-1} * y 
        call complex_leqs_solver(2*nstop,4*nstop,S,y,r)
        ! 
        ! Finally, I obtain mean-field scattering coefficients
        ! dd(1,n) : d^(1)_n
        ! dd(2,n) : d^(2)_n
        !
        do n=1,nmax
           dd(1,n) = r(2*n-1) 
           dd(2,n) = r(2*n)  
        enddo
        deallocate(T,S,y,r)
endif 

if(iquiet .eq. 0) then
print *, '----- Lorenz-Mie coefficients -----'
write(*,fmt='(A3,4A15)') "n","Re(an)","Im(an)","Re(bn)","Im(bn)"
do n=1,nstop
  write(*,fmt='(I3,1P4E15.5)') n,real(ad(1,n)),aimag(ad(1,n)),real(ad(2,n)),aimag(ad(2,n))
enddo
print *, '----- Mean field coefficients -----'
write(*,fmt='(A3,4A15)') "n","Re(d1n(1))","Im(d1n(1))","Re(d1n(2))","Im(d1n(2))"
do n=1,nstop
  write(*,fmt='(I3,1P4E15.5)') n,real(dd(1,n)),aimag(dd(1,n)),real(dd(2,n)),aimag(dd(2,n))
enddo
endif

!
! Replacing the lorenz-mie scattering coefficients by the mean-field 
! scattering coefficients, and find scatternig amplitude of each
! monomer.
!
dang=halfpi/real(nang-1,kind=dp)
call renorm_mie(dd,nstop,nang,dang,S1,S2)

!
! Compute scattering matrix elements based on scattering amplitude
! of each monomer and the static structure factor of an aggregate
!
Smat = 0.0_dp
do j=1,2*nang-1
        ! scattering angle in radian
        ang(j)  = dang * real(j-1,kind=dp)        
        ! magnitude of scattering vector
        q       = 2.0_dp*k*sin(0.5_dp*ang(j))     
        ! scattering matrix elements of each monomer
        S11     =       0.5_dp*ABS(S2(j))*ABS(S2(j))
        S11     = S11 + 0.5_dp*ABS(S1(j))*ABS(S1(j))
        S12     =       0.5_dp*ABS(S2(j))*ABS(S2(j))
        S12     = S12 - 0.5_dp*ABS(S1(j))*ABS(S1(j))       
        S33     =       DBLE( S1(j)*CONJG(S2(j)) )
        S34     =       AIMAG( S2(j)*CONJG(S1(j)) )
        !
        ! Compute interference of scattered waves based on a statistical 
        ! distribution of monomers specified by the two-point correlation
        ! function. The inteference pattern is characterized by the quantity
        ! called the static structure factor (Sq) as defined in Eq (17) in
        ! Tazaki & Tanaka (2018). 
        !
        ! The analytical expression of the static structure factor is available 
        ! for the case of iqcor = 1 and 2, but not for iqcor =3. 
        ! Hence, iqcor=3 requires additional numerical integration.
        !
        if(iqcor .eq. 1) then
                ! compute confluent hypergeometric function
                if(df .eq. 3.0_dp) then
                        if(q*q*Rg*Rg/3.0 .ge. eta) then
                                Sq = 0.0_dp
                        else
                                Sq = exp(-q*q*Rg*Rg/3.0_dp)
                        endif
                else
                        AL = 0.5_dp * df
                        BB = 1.5_dp
                        XX = -q*q*Rg*Rg/df
                        if(q*Rg .lt. qRg_crit) then
                                call chgm(AL,BB,XX,Sq)
                        else
                                CC=sqrt(pi)*Df**(0.5_dp*Df)/(2.0*gamma(0.5_dp*(3.0-Df)))
                                Sq=CC*(q*RG)**(-Df)
                        endif
                endif
         elseif(iqcor .eq. 2) then
                xi =  sqrt(2.0_dp/(df*(df+1.0_dp))) * Rg
                IF(J .eq. 1) then
                        Sq = 1.0_dp
                else
                        Sq = sin((df-1.0_dp)*atan(q*xi))/(q*xi*(df-1.0_dp)&
                             *(1.0_dp+(q*xi)**2.0_dp)**((df-1.0_dp)/2.0_dp))
                endif
          elseif(iqcor .eq. 3) then
                call structure_factor_integration(df,q,Rg,Sq)
          endif
          
          Smat(1,j) = PN * S11 * (1.0_dp + (PN-1.0_dp) * Sq)
          Smat(2,j) = PN * S12 * (1.0_dp + (PN-1.0_dp) * Sq)
          Smat(3,j) = PN * S33 * (1.0_dp + (PN-1.0_dp) * Sq)
          Smat(4,j) = PN * S34 * (1.0_dp + (PN-1.0_dp) * Sq)

          ! safety check.
          if(Smat(1,j) .lt. 0.d0) then
                  print *, 'ERROR: S11 for aggregte becomes negative. Strange!'
                  print *, '       STOP.'
                  stop
          endif
enddo

!
! Compute scattering cross section of aggregates
! by integrating S11 via Simpson's rule.
!
Csca = 0.0_dp
do j=1,nang-1
        Csca  = Csca + dang * &
                       ( Smat(1,2*j-1) * sin(ang(2*j-1))   + &
                  4.0_dp*Smat(1,2*j  ) * sin(ang(2*j  ))   + &
                         Smat(1,2*j+1) * sin(ang(2*j+1)) ) / 3.0_dp
enddo
Csca  = twopi     * Csca / k / k
PF(:) = Smat(1,:) / Csca / k / k
!
! asymmetry parameter by integrating
! phase function with Simpson's rule.
! 
g_asym =  0.0_dp
nrm    =  0.0_dp
DO j=1,nang-1
        g_asym = g_asym + dang * &
                       ( PF(2*j-1) * sin(ang(2*j-1)) * cos(ang(2*j-1)) + &
                  4.0_dp*PF(2*j  ) * sin(ang(2*j  )) * cos(ang(2*j  )) + &
                         PF(2*j+1) * sin(ang(2*j+1)) * cos(ang(2*j+1)) ) / 3.0d0

        nrm    =  nrm   + dang * &
                       ( PF(2*j-1) * sin(ang(2*j-1))   + &
                  4.0_dp*PF(2*j  ) * sin(ang(2*j  ))   + &
                         PF(2*j+1) * sin(ang(2*j+1)) ) / 3.0_dp
END DO
nrm     = twopi * nrm
g_asym  = twopi * g_asym

!
! Check normalization of phase function
!
if (abs(nrm-1.0_dp) .gt. 1.0e-3_dp) then
        write(*,*) ABS(nrm-1.0_dp)*1.0e2_dp
        print *, 'ERROR: Phase function is not normalized to unity. Strange!'
        stop
end if

!
! Calculate optical cross sections.
!
! For the case of an isolcated individual monomer, 
! each cross section can be written by 
! (Equation (4.61; 4.62) in Bohren & Huffman textbook.) 
!
!                   nstop
!            2*pi   ---
!  C_ext =  ------  \    (2n+1) Re ( an + bn )
!            k^2    /--
!                   n=1
!
!                   nstop
!            2*pi   ---
!  C_sca =  ------  \    (2n+1) ( |an|^2 + |bn|^2 )
!            k^2    /--
!                   n=1
!
!                   nstop
!            2*pi   ---                      1                    1
!  C_abs =  ------  \    (2n+1) Re{|an|^2 ( --- - 1 ) + |bn|^2 ( --- - 1 )}
!            k^2    /--                     an*                  bn*
!                   n=1
! 
! where (an,bn) is scattering coefficients of each monomer and * denotes 
! the complex conjugate.
!
if(iqsca .eq. 1) then
        Cabsp = 0.0_dp
        do n=1,nstop
                cn1   = 1.0_dp/conjg(ad(1,n))-1.0_dp
                cn2   = 1.0_dp/conjg(ad(2,n))-1.0_dp
                Cabsp = Cabsp + real(2*n+1,kind=dp)*&
                       &real(cn1*abs(ad(1,n))*abs(ad(1,n))+cn2*abs(ad(2,n))*abs(ad(2,n)))
        enddo
        Cabsp = Cabsp * twopi * PN / k / k
        Cext  = Cabsp + Csca
elseif(iqsca .eq. 2) then
        Cext = 0.0_dp
        do n=1,nstop
                Cext = Cext + real(2*n+1,kind=dp)*real(dd(1,n)+dd(2,n))
        enddo
        Cext  = Cext * twopi * PN / k / k
        Cabsp = Cext - Csca        
elseif(iqsca .eq. 3) then
        Cext  = 0.0_dp
        Cabsp = 0.0_dp
        ! Cext  : Mean-field extinction
        ! Cabsp : RGD theory
        do n=1,nstop
                Cext = Cext + real(2*n+1,kind=dp)*real(dd(1,n)+dd(2,n))
                cn1   = 1.0_dp/conjg(ad(1,n))-1.0_dp
                cn2   = 1.0_dp/conjg(ad(2,n))-1.0_dp
                Cabsp = Cabsp + real(2*n+1,kind=dp)*&
                       &real(cn1*abs(ad(1,n))*abs(ad(1,n))+cn2*abs(ad(2,n))*abs(ad(2,n)))
        enddo
        Cext  = Cext  * twopi * PN / k / k 
        Cabsp = Cabsp * twopi * PN / k / k 
        !
        ! Geometrical cross section of an aggregate
        !
        if(iqgeo .eq. 1) then
                GC = pi * RC * RC
        elseif(iqgeo .eq. 2) then
                call geocross(PN,R0,RC,GC)
        elseif(iqgeo .eq. 3) then
                call geocross_mft(PN,df,GC)
                GC = GC * PN * pi * R0 * R0
        endif
        !
        ! Compute "tau" of an aggregate, where Cabsp is RGD theory's one
        !
        tau   = (Cabsp/GC)                      
        !
        ! RGD theory with geometrical optics approximation.
        ! TIP: To avoid underflow, tau>eta=10, Cabsp = GC.
        !
        if(tau .ge. 10.0_dp) then 
                Cabsp = GC
        else
                Cabsp = GC * (1.0_dp - exp(-tau))      
        endif
        ! Compare MFT vs. RGD+GeoOpt, take the larger one.
        Cabsp = max(Cabsp,Cext-Csca)            
        ! Scattrering is computed from Cext(MFT)-Cabs.
        Csca  = Cext - Cabsp                    
endif

deallocate(an,bn,ad,dd)

return
end subroutine meanscatt

!--------------------------------------------------------------------------------
!
!  This subroutine finds effective refractive index of an aggregates
!  by using Maxwell--Garnett mixing rule.
!
!  eps1: dielectric function of material 1
!  eps2: dielectric function of material 2 (vacuum)
!  f1  : volume filling factor of material 1
!  f2  : volume filling factor of material 2 ( = Porosity )
!  
!--------------------------------------------------------------------------------
subroutine MGmixing(refrel,mgav,f1)
use types
implicit none
real(kind=dp)::f1
complex(kind=dp)::eps1,eps2,refrel
complex(kind=dp)::mg,mgav
eps1 = refrel*refrel
eps2 = cmplx(1.0_dp,0.0_dp,kind=dp)
mg=eps2*(2.0_dp*f1*(eps1-eps2)+eps1+2.0_dp*eps2)/&
         (eps1+2.0_dp*eps2-f1*(eps1-eps2))
mgav=sqrt(mg)
return
end subroutine MGmixing

!--------------------------------------------------------------------------------
!  
!  This subroutine performs integration of the static structure factor by 
!  Botet et al. (1995) based on  df-dependent cut-off function for two-point 
!  correlation function.
!  The statis structure factor S(q) is given by
!
!
!          c*df     /x_max
!  S(q) =  ----  *  |  dx x^{df-2}*sin(q*Rg*x)*exp[-c*x^df],
!          q*Rg     /0
!
! 
!  where x = u/Rg; where Rg is the radius of gyration and u is the 
!  distance between two monomers; q is the magnitude of scattering
!  vector, df is the fractal dimension, c=0.5 is normalization factor.
!
!  Since we have a steep cut-off function for x>1, it would be enough 
!  to take x_max slightly larger than unity. 
!  Here we determine xmax so that eta=c*x^df becomes 25.
!  
!  We use Simpson's rule for the numerical integration, and then,
!  the integration range [0,xmax] will be divided by 2*m points.
!
!  Caution:
!       For large df(>2) and large qRg, the numerical integration sometimes 
!       returns a negative value of S(q). Since S(q) is a power spectrum,
!       it must be positive. I also noticed that the same problem arises with 
!       the code written by Pascal Rannou. Because of this ill behavior, I 
!       prefer to use Gausian-cut-off model, which is more stable.
!    
!       A not-too-bad tentative prescription is simply set S(q)=0 when S(q)<0, 
!       which means scattered waves from each monomer is incoherent. This is 
!       what I did in if-sentence at the end of this subroutine.
!
!--------------------------------------------------------------------------------
subroutine structure_factor_integration(D,q,Rg,Sq)
use types
implicit none
integer::i
integer      , parameter :: m    = 100000    ! # of integration grids 
real(kind=dp), parameter :: c    = 0.5_dp    ! normalization factor by Botet+97
real(kind=dp)::xmax,D,q,Rg,Sq,h
real(kind=dp)::func,x0,x1,x2,eta

! estimate upper limit of integration: xmax
eta = 25
xmax = (eta/c)**(1.0_dp/D)

if(q .eq. 0.0_dp) then
        Sq = 1.0_dp
else
        Sq = 0.0_dp
        h = xmax/real(2*m,kind=dp)
        do i=0,m-1
                x0 = real(2*i,kind=dp)*h
                x1 = real(2*i+1,kind=dp)*h
                x2 = real(2*i+2,kind=dp)*h
                Sq = Sq + h * ( func(x0,c,D,q,Rg) +&
                         4.0_dp*func(x1,c,D,q,Rg) +&
                                func(x2,c,D,q,Rg) )/3.0d0
        end do
        Sq  = Sq * c * D / q / Rg
endif

! switch-off the monomer-monomer correlation when S(q)<0.
if(Sq .lt. 0.0_dp) then
        Sq = 0.0_dp
endif

return
end subroutine structure_factor_integration

function func(x,c,D,q,Rg)
use types
implicit none
real(kind=dp), parameter :: floorval = 1.e-30_dp
real(kind=dp)::q,Rg,D,c,x
real(kind=dp)::func
if(x .eq. 0.0_dp) then
        func = 0.0_dp
else
        func = sin(q*Rg*x)*exp(-c*x**D)*x**(D-2.0_dp)
        if(abs(func) .le. floorval) then
                func = 0.0_dp
        endif
endif
return
end function func

!--------------------------------------------------------------------------------
!
! Calculate Lorenz-Mie scattering coefficients (an,bn) for a monomer particle.
! 
! Since monomer's size parameter is supposed not to be very large, 
! I use simple Bohren & Huffman Mie algorithm is used. 
!
! The original BHMIE code is taken from Bruce Draine's HP:
!       https://www.astro.princeton.edu/~draine/scattering.html
! although we have slightly modified it.
!
!--------------------------------------------------------------------------------
subroutine lorenz_mie(x,refrel,nstop,a,b)
use types
implicit none
integer,parameter                 :: nmxx=150000
integer                           :: nstop
real(kind=dp)                     :: x
complex(kind=dp)                  :: refrel
!--------------------------------------------------------------------------------
integer                           :: n,nmx,en
real(kind=dp)                     :: xstop,ymod,enr,nr
real(kind=dp)                     :: psi0,psi1,psi,chi0,chi1,chi
complex(kind=dp)                  :: xi1,xi,y
!complex(kind=dp),dimension(nmxx)  :: d
complex(kind=dp),allocatable,dimension(:)  :: d
complex(kind=dp),dimension(nstop) :: a,b

y     = refrel * x
ymod  = abs(y)
xstop = x + 4.0_dp * x ** (1.0_dp/3.0_dp) + 2.0_dp
nmx   = nint(max(xstop,ymod)) + 15
if(nmx .gt. nmxx) then
        write(*,*) "Error: nmx > nmxx=",nmxx,"for |m|x=",ymod
        stop
endif

! calculate logarithmic derivative D_n(mx) 
! by downward recurrence.
! Initial value is set as D(mx)=0+0i at n=NMX.
allocate(d(1:nmx))
d(nmx)=cmplx(0.0_dp,0.0_dp,kind=dp)
do n=1,nmx-1
        en  = nmx-n+1
        enr = real(nmx-n+1,kind=dp)
        d(nmx-n) = (enr/y) - (1.0_dp/(d(en) + enr/y))
enddo

!Initial values
psi0 =  cos(x)
psi1 =  sin(x)
chi0 = -sin(x)
chi1 =  cos(x)
xi1  =  cmplx(psi1,-chi1,kind=dp)

do n=1,nstop
        nr = real(n,kind=dp)
        !Calculate psi and chi via upward recurrence:
        psi = (2.0_dp*nr-1.0_dp)*psi1/x - psi0
        chi = (2.0_dp*nr-1.0_dp)*chi1/x - chi0
        xi  = cmplx(psi,-chi,kind=dp)
        !Calculate Lorentz-Mie coefficients:
        a(n) = (d(n)/refrel+nr/x)*psi-psi1
        a(n) = a(n)/((d(n)/refrel+(nr/x))*xi-xi1)
        b(n) = (refrel*d(n)+nr/x)*psi-psi1
        b(n) = b(n)/((refrel*d(n)+(nr/x))*xi-xi1)
        !Preparation for next step
        psi0 = psi1
        psi1 = psi
        chi0 = chi1
        chi1 = chi
        xi1 = cmplx(psi1,-chi1,kind=dp)
enddo

deallocate(d)

return
end subroutine lorenz_mie

!--------------------------------------------------------------------------------
!
! Calculate scattering amplitude S1, S2 from given scattering coefficients d.
!
! This subroutine is based on the BHMIE code is taken from Bruce Draine's HP:
!       https://www.astro.princeton.edu/~draine/scattering.html
! although we have slightly modified it.
!
!--------------------------------------------------------------------------------
subroutine renorm_mie(d,nstop,nang,dang,s1,s2)
use types
implicit none
integer                                 :: nang,j,jj,n,nn,nstop
real(kind=dp)                           :: dang,theta 
real(kind=dp)                           :: en,fn,p
real(kind=dp),dimension(1:nang)         :: amu,pi0,pi1,pi,tau
complex(kind=dp)                        :: an,bn,an1,bn1
complex(kind=dp),dimension(1:2*nang-1)  :: S1,S2
complex(kind=dp),dimension(2,nstop)     :: d

!*** require nang.ge.1 in order to calculate scattering intensities
do j=1,nang
        theta=dble(j-1)*dang
        amu(j)=cos(theta)
enddo
do j=1,nang
        pi0(j)=0.0_dp
        pi1(j)=1.0_dp
enddo
nn=2*nang-1
do j=1,nn
        S1(j)=cmplx(0.0_dp,0.0_dp,kind=dp)
        S2(j)=cmplx(0.0_dp,0.0_dp,kind=dp)
enddo

p   = -1.0_dp 
an  = cmplx(0.0_dp,0.0_dp,kind=dp)
an1 = cmplx(0.0_dp,0.0_dp,kind=dp)
bn  = cmplx(0.0_dp,0.0_dp,kind=dp)
bn1 = cmplx(0.0_dp,0.0_dp,kind=dp)

do n=1,nstop
        en = real(n,kind=dp)
        fn = (2.0_dp*en+1.0_dp)/(en*(en+1.0_dp))

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

        !c*** now calculate scattering intensity pattern
        !    first do angles from 0 to 90
        do j=1,nang
                !jj=2*nang-j
                pi(j)=pi1(j)
                tau(j)=en*amu(j)*pi(j)-(en+1.0_dp)*pi0(j)
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
                pi1(j)=((2.0_dp*en+1.0_dp)*amu(j)*pi(j)-(en+1.0_dp)*pi0(j))/en
                pi0(j)=pi(j)
        enddo
enddo

return
end subroutine renorm_mie

!--------------------------------------------------------------------------------
!
! LPMNS computes associated Legendre functions Pmn(X) and derivatives P'mn(x).
!
!  Licensing:
!
!    This routine is copyrighted by Shanjie Zhang and Jianming Jin.  However, 
!    they give permission to incorporate this routine into a user program 
!    provided that the copyright is acknowledged.
!
!  Modified:
!
!    18 July 2012
!
!  Author:
!
!    Shanjie Zhang, Jianming Jin
!
!  Reference:
!
!    Shanjie Zhang, Jianming Jin,
!    Computation of Special Functions,
!    Wiley, 1996,
!    ISBN: 0-471-11963-6,
!    LC: QA351.C45.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, the order of Pmn(x).
!
!    Input, integer ( kind = 4 ) N, the degree of Pmn(x).
!
!    Input, real ( kind = 8 ) X, the argument.
!
!    Output, real ( kind = 8 ) PM(0:N), PD(0:N), the values and derivatives
!    of the function from degree 0 to N.
!
!--------------------------------------------------------------------------------
!
!   Note by R.T.
!   .f90 version was downloaded from 
!   https://people.sc.fsu.edu/~jburkardt/f_src/special_functions/special_functions.html
!   
!   Revision History 
!       2020. Oct 31:   - Minor change so that a real type is specified by "dp".
!                       - D+00 --> _dp
!
!--------------------------------------------------------------------------------
subroutine lpmns ( m, n, x, pm, pd )
use types
implicit none
integer ( kind = 4 ) n
integer ( kind = 4 ) k
integer ( kind = 4 ) m
real ( kind = dp ) pm(0:n)
real ( kind = dp ) pm0
real ( kind = dp ) pm1
real ( kind = dp ) pm2
real ( kind = dp ) pmk
real ( kind = dp ) pd(0:n)
real ( kind = dp ) x
real ( kind = dp ) x0

do k = 0, n
pm(k) = 0.0_dp
pd(k) = 0.0_dp
end do

if ( abs ( x ) == 1.0_dp ) then
        do k = 0, n
        if ( m == 0 ) then
                pm(k) = 1.0_dp
                pd(k) = 0.5_dp * k * ( k + 1.0_dp )
                if ( x < 0.0_dp ) then
                        pm(k) = ( -1.0_dp ) ** k * pm(k)
                        pd(k) = ( -1.0_dp ) ** ( k + 1 ) * pd(k)
                end if
        else if ( m == 1 ) then
                pd(k) = 1.0D+300
        else if ( m == 2 ) then
                pd(k) = -0.25_dp * ( k + 2.0_dp ) * ( k + 1.0_dp ) &
                  * k * ( k - 1.0_dp )
                if ( x < 0.0_dp ) then
                        pd(k) = ( -1.0_dp ) ** ( k + 1 ) * pd(k)
                end if
        end if
        end do
        return
end if

x0 = abs ( 1.0_dp - x * x )
pm0 = 1.0_dp
pmk = pm0
do k = 1, m
        pmk = ( 2.0_dp * k - 1.0_dp ) * sqrt ( x0 ) * pm0
        pm0 = pmk
end do
pm1 = ( 2.0_dp * m + 1.0_dp ) * x * pm0
pm(m) = pmk
pm(m+1) = pm1
do k = m + 2, n
        pm2 = ( ( 2.0_dp * k - 1.0_dp ) * x * pm1 &
        - ( k + m - 1.0_dp ) * pmk ) / ( k - m )
        pm(k) = pm2
        pmk = pm1
        pm1 = pm2
end do

pd(0) = ( ( 1.0_dp - m ) * pm(1) - x * pm(0) ) &
/ ( x * x - 1.0_dp )  
do k = 1, n
        pd(k) = ( k * x * pm(k) - ( k + m ) * pm(k-1) ) &
        / ( x * x - 1.0_dp )
end do

return
end subroutine lpmns

!--------------------------------------------------------------------------------
!
! LPN computes Legendre polynomials Pn(x) and derivatives Pn'(x).
!
!  Licensing:
!
!    This routine is copyrighted by Shanjie Zhang and Jianming Jin.  However, 
!    they give permission to incorporate this routine into a user program 
!    provided that the copyright is acknowledged.
!
!  Modified:
!
!    07 July 2012
!
!  Author:
!
!    Shanjie Zhang, Jianming Jin
!
!  Reference:
!
!    Shanjie Zhang, Jianming Jin,
!    Computation of Special Functions,
!    Wiley, 1996,
!    ISBN: 0-471-11963-6,
!    LC: QA351.C45.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the maximum degree.
!
!    Input, real ( kind = 8 ) X, the argument.
!
!    Output, real ( kind = 8 ) PN(0:N), PD(0:N), the values and derivatives
!    of the polyomials of degrees 0 to N at X.
!
!--------------------------------------------------------------------------------
!
!   Note by R.T.
!   .f90 version was downloaded from 
!   https://people.sc.fsu.edu/~jburkardt/f_src/special_functions/special_functions.html
!   
!   Revision History 
!       2020. Oct 31:   - Minor change so that a real type is specified by "dp".
!                       - D+00 --> _dp
!
!--------------------------------------------------------------------------------
subroutine lpn ( n, x, pn, pd )
use types
implicit none
integer ( kind = 4 ) n
integer ( kind = 4 ) k
real ( kind = dp ) p0
real ( kind = dp ) p1
real ( kind = dp ) pd(0:n)
real ( kind = dp ) pf
real ( kind = dp ) pn(0:n)
real ( kind = dp ) x
pn(0) = 1.0_dp
pn(1) = x
pd(0) = 0.0_dp
pd(1) = 1.0_dp
p0 = 1.0_dp
p1 = x

do k = 2, n

        pf = ( 2.0_dp * k - 1.0_dp ) / k * x * p1 &
        - ( k - 1.0_dp ) / k * p0
        pn(k) = pf

        if ( abs ( x ) == 1.0_dp ) then
                pd(k) = 0.5_dp * x ** ( k + 1 ) * k * ( k + 1.0_dp )
        else
                pd(k) = k * ( p1 - x * pf ) / ( 1.0_dp - x * x )
        end if

        p0 = p1
        p1 = pf

end do

return
end subroutine lpn

!--------------------------------------------------------------------------------
!
!  This subroutine finds the inverse matrix of ( n , n ) complex matrix.
!  I adopt a straightforward way to do this, although the method is 
!  inefficient in terms of memory (Press et al. 1992).
!
!  Consider a set of complex n linear equations:
!
!        T * r = p,
!
!  where T is (n,n) matrix with complex elements, r and p are the 
!  complex column vector of size n, and * denotes the inner product.
!  This complex systems of equations can be rewritten in the form:
!
!       ( Re(T) + i Im(T) ) * ( Re(r) + i Im(r) ) = ( Re(p) + i Im(p) ),
!
!  where A, C are the (n x n) matrix with real elements, x,y,b,d are
!  column vectors of size n. This linear equations can be reduced to  
!
!        Re(T) * Re(r) - Im(T) * Im(r) = Re(p),
!        Im(T) * Re(r) + Re(T) * Im(r) = Im(p),
!
!  As a result, we obtain (2n,2n) set of linear equations with real elements.
!
!        /               \       /       \       /       \
!       |  Re(T)  -Im(T)  |     |  Re(r)  |     |  Re(p)  |
!       |                 |  *  |         |  =  |         |
!       |  Im(T)   Re(T)  |     |  Im(r)  |     |  Im(p)  |
!        \               /       \       /       \       /
!   
!  This (2n,2n) real matrix can be inverted by a common method
!  with LU decomposition and LU backsubstitution.
!
!--------------------------------------------------------------------------------
subroutine complex_leqs_solver(n,np,T,p,r)
use types
implicit none
integer::i,j,n,np
integer,dimension(np)::indx
complex(kind=dp),dimension(n,n)::T
complex(kind=dp),dimension(n)::p,r
real(kind=dp),dimension(np,np)::a
real(kind=dp),dimension(np)::b
real(kind=dp)::d 

!---------------------------------------
! Create np x np = 2n x 2n set of "real" equations
! where n = 2 * nstop
!---------------------------------------
a = 0.0_dp
b = 0.0_dp
r = 0.0_dp
do i=1,n
        do j=1,n
                a(i,j)     =   real(T(i,j))
                a(i,j+n)   = -aimag(T(i,j))
                a(i+n,j)   =  aimag(T(i,j))
                a(i+n,j+n) =   real(T(i,j))
        enddo
        b(i)   = real(p(i))
        b(n+i) = aimag(p(i))
enddo
call ludcmp(a,np,np,indx,d)
call lubksb(a,np,np,indx,b)
do i=1,n
        r(i) = cmplx(b(i),b(n+i),kind=dp)
enddo

return
end subroutine complex_leqs_solver

!--------------------------------------------------------------------------------
!
!  This subroutine performs integration of S_p(kRg) (Equation 31):
!
!               pi^2     /Infinity
!  S_p(k*Rg) =  ----  *  |  du  u * J_{p+1/2}(u) * H_{p+1/2}^(1)(u) * g(u/k), 
!                k^3     /0
!
!  where g(u) is the two-point correlation function (Equations 18 and 19): 
!
!               1       / u  \ ^{d_f-3}       / u  \
!  g(u) =   ---------  |----- |       *  fc  | ---- |,
!           4*pi*Rg^3   \ Rg /                \ Rg /
!
!  where fc is the cut-off function. By substituting g(u) into S_p, we have
!
!                 pi   /u_max
!  S_p(k*Rg) = ------ *|  du  u^{df-2} * J_{p+1/2}(u) * H_{p+1/2}^(1)(u) * fc(u/xg), 
!               4Xg^df /u_min
!
!  where the integration range is approximated by the range [u_min,u_max]. 
!  By using the spherical Bessel j_p(u) and Hankel functions of 1st kind h_p^(1)(u),
!  the Bessel function and the Hankel function are rewritten by
!
!               J_{p+1/2}    (u) = sqrt(2u/pi) * j_p(u) 
!               H_{p+1/2}^(1)(u) = sqrt(2u/pi) * h_p^(1)(u) 
!
!  we have
!
!                 1        /u_max
!  S_p(k*Rg) =  ------  *  |  du  u^{df-1} * j_{p}(u) * h_{p}^(1)(u) * fc(u/xg), 
!               2*Xg^df    /u_min
!
!  For the unitary condition of the two-point correlation function is
!
!                /Infinity              
!  1         =   |  dw  4*pi*w^2 g(w),
!                /0                     
!
!  If we take the integration variable as w = u/k, then we obtain 
!
!                1     /u_max
!  1         = ------  |  du u^{df-1} fc(u/xg),    .... Eq. (*)
!               Xg^df  /u_min
!
!
!  The integration range [umin,umax] is determined as follows.
!  The integrand of Equation (*) is
!       
!         (u/xg)^{df-1}fc(u/xg) du = (u/xg) ^{d_f}fc(u/xg) dlnu
! 
!  u_max is chosen so that fc(u/xg) ~ exp[-eta1].
!       For iqcor=1 (Gauss)
!               u_max ~ 2 * xg * sqrt(eta1 / d_f)
!       For iqcor=1 (Exponential)
!               u_max ~ xg * eta1 * sqrt(2.0 /(d_f*(d_f+1))
!       For iqcor=1 (FLDIM)
!               u_max ~ xg * sqrt( 2.0 * eta1 ) ** (1.0/d_f)
!
!  u_min is chosen so that (u_min/xg)^{d_f} ~ exp[-eta2], thus, 
!
!               umin ~ xg * exp(-eta2/d_f)
!
!  I adopt eta1 = 25.0, and eta2 = 40.0, respectively.
!
!--------------------------------------------------------------------------------
subroutine integration_of_Sp(iqcor,D,p,xg,Sp)
use types; use const
!integer,parameter      :: nn =  10000                    !  integration grid.
!--------------------------------------------------------------------------------
! Old boundary given by Jablonski et al. 1994 : !!! NOT FAVORED !!!
!--------------------------------------------------------------------------------
!real(kind=dp),parameter::a1=-3.6693021d-8
!real(kind=dp),parameter::a2=-3.1745158d-5
!real(kind=dp),parameter::a3=2.1567720d-2
!real(kind=dp),parameter::a4=9.9123380d-1
!real(kind=dp),parameter::b1=-9.9921351d-8
!real(kind=dp),parameter::b2=8.7822303d-6
!real(kind=dp),parameter::b3=1.0238752d-2
!real(kind=dp),parameter::b4=3.7588265
!--------------------------------------------------------------------------------
! New boundary used in by Tazaki & Tanaka 2018
!--------------------------------------------------------------------------------
real(kind=dp),parameter:: a1 =  1.69496268177237e-08 
real(kind=dp),parameter:: a2 = -2.43299782942114e-05
real(kind=dp),parameter:: a3 =  0.0158750501131321 
real(kind=dp),parameter:: a4 =  1.00672154148706
real(kind=dp),parameter:: b1 =  2.49355951047228e-08 
real(kind=dp),parameter:: b2 = -2.9387731648675e-05 
real(kind=dp),parameter:: b3 =  0.0135005554796179 
real(kind=dp),parameter:: b4 =  3.72312019844119
!--------------------------------------------------------------------------------
real(kind=dp),parameter:: floorvalue=1.0e-30_dp
real(kind=dp),parameter:: eta1 = 25.0_dp
real(kind=dp),parameter:: eta2 = 40.0_dp
integer                :: n,p,iqcor,isol
real(kind=dp)::umin,umax,du,xg,D,fc
real(kind=dp)::lnxa,lnxb,jp,yp,unitary,error
real(kind=dp),allocatable,dimension(:)::u
real(kind=dp),allocatable,dimension(:)::intg_unit
real(kind=dp),allocatable,dimension(:)::SJ,SY
complex(kind=dp),allocatable,dimension(:)::intg
complex(kind=dp)::hp,wa,Sp

! estimate umax, umin from cut-off function.
umax = 0.0_dp
if(iqcor .eq. 1) then
        umax = 2.0_dp * xg * sqrt(eta1/D) 
elseif(iqcor .eq. 2) then
        umax = xg * eta1 * sqrt(2.0_dp/(D*(D+1.0_dp)))
elseif(iqcor .eq. 3) then
        umax = xg * (2.0_dp*eta1)**(1.0_dp/D)
endif
umin = xg * exp(-eta2/D)
du   = (umax/umin) ** (1.0_dp/real(nn-1,kind=dp))



allocate(u(1:nn),intg(1:nn),intg_unit(1:nn))
do n=1,nn
        u(n) = umin * du ** real(n-1,kind=dp)
enddo

! calculate boundary for Bessel solver.
lnxa=a1*dble(p)**3.0_dp+a2*dble(p)**2.0_dp+a3*dble(p)+a4
lnxb=b1*dble(p)**3.0_dp+b2*dble(p)**2.0_dp+b3*dble(p)+b4

intg      = cmplx(0.0_dp,0.0_dp,kind=dp)
intg_unit = 0.0_dp
do n=1,nn
        if(log(u(n)) .lt. lnxa) then
                isol = 1
        elseif(log(u(n)) .gt. lnxb) then
                isol = 3
        else
                isol = 2
        endif
        allocate(SJ(0:p),SY(0:p))
        call sphbessel(p,u(n),SJ,SY,isol)
        jp=SJ(p)
        yp=SY(p)
        hp = cmplx(jp,yp,kind=dp)
        deallocate(SJ,SY)

        if(jp*0.0_dp /= 0.0_dp .or. isnan(jp) .eqv. .true.) then
               write(*,*) "--------------------------------------------------------"
               write(*,*) "Error: NaN is found at the spherical Bessel function of the first kind"
               write(*,*) "argument u : ",u(n)
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
        elseif(yp*0.0_dp /= 0.0_dp .or. isnan(yp) .eqv. .true.) then
               write(*,*) "--------------------------------------------------------"
               write(*,*) "Error: NaN is found at the spherical Bessel function of the second kind"
               write(*,*) "argument u : ",u(n)
               write(*,*) "order    p : ",p
               write(*,*) "y_p(x)     : ",yp 
               write(*,*) "Simulation is aborted !"
               write(*,*) "--------------------------------------------------------"
               stop
        endif
        intg(n)      = u(n) ** (D-1.0_dp) * jp * hp * fc(iqcor,u(n),xg,D)
        intg_unit(n) = u(n) ** (D-1.0_dp) * fc(iqcor,u(n),xg,D)
enddo
wa      = cmplx(0.0_dp,0.0_dp,kind=dp)
unitary = 0.0_dp
do n=1,nn-1
        wa      = wa      + 0.5_dp * (intg(n)+intg(n+1)) * (u(n+1)-u(n))
        unitary = unitary + 0.5_dp * (intg_unit(n)+intg_unit(n+1))*(u(n+1)-u(n))
enddo

unitary = unitary / xg ** D
Sp      = 0.5_dp * wa / xg ** D
error   = abs(1.0_dp-unitary)
        
!
! Set the floorvalue:
! For large p and small xg, the real part of Sp will be very small. 
!
if(abs(real(Sp)) .lt. floorvalue) then
        Sp = cmplx(floorvalue,aimag(Sp))
endif
if(abs(aimag(Sp)) .lt. floorvalue) then
        Sp = cmplx(real(Sp),-floorvalue)
endif
!

if(error .ge. 1.0e-3_dp) then
        write(*,*) "--------------------------------------------------------"
        write(*,*) "Check unitary condition : the two-points correlation function"
        write(*,*) "Numerical integration of g(u) with error of ",error*1.d2," (%)"
        write(*,*) "which exceed 0.1%. This means sp(kRg) integration"
        write(*,*) "may not converge."
        write(*,*) "Simulation is aborted !"
        write(*,*) "--------------------------------------------------------"
        stop
endif

deallocate(u,intg,intg_unit)


return
end subroutine integration_of_Sp

!--------------------------------------------------------------------------------
!
! The cut-off models for the two-point correlation function
! (see Equation (19) of Tazaki & Tanaka 2018).
!
! For the Gaussian cut-off model (iqcor = 1):
! 
!        / u  \       2c^{df/2}       /     / u  \^2 \
!     fc |----|   = ------------ exp | - c |-----|   |, with c = df/4,
!        \ xg /      GAMMA(df/2)      \     \ Xg /   /
!
! where, GAMMA is the gamma function.
! For the exponential cut-off model (iqcor =2):
!
!        / u  \       c^{df}          /     / u \  \
!     fc |----|   = ------------ exp | - c |-----| |,  with c = [df*(df+1)/2]^{1/2},
!        \ xg /      GAMMA(df)        \     \ Xg / /
!
!
! for the fractal-dimension cut-off model (iqcor=3):
!
!        / u  \                  /      / u  \^{df} \
!     fc |----|   =  c * df exp | - c  |-----|      |,  with c = 0.5
!        \ xg /                  \      \ Xg /      /
!
!--------------------------------------------------------------------------------
function fc(iqcor,u,x,df)
use types
implicit none
integer::iqcor
real(kind=dp)::u,x,df
real(kind=dp)::c,fc
! initialize
fc = 0.0_dp
if(iqcor .eq. 1) then
        c  = 0.25_dp * df
        fc = (2.0d0*c**(0.5_dp*df)/gamma(0.5_dp*df))*exp(-c*u*u/x/x)
elseif(iqcor .eq. 2) then
        c  = sqrt(0.5_dp*df*(df+1.0d0))
        fc = c**df/gamma(df) * exp(-c*u/x)
elseif(iqcor .eq. 3) then
        c  = 0.5d0
        fc = c*df*exp(-c*(u/x)**df)
endif
return
end function fc

!--------------------------------------------------------------------------------
!
! This subroutine computes spherical bessel function of the first kind j_p(x)
! and the second kind y_p(x) at a given order p and the argument x.
!
! Calculations of the spherical Bessel function of the first kind j_p(x) for a 
! wide range of parameters is not straightforward because of numerical instability.
! However, Jablonski (1994) reported that a proper combination of (i) series 
! expansion,(ii) downward recurrence, and (ii) upward recuurence can yield
! correct results for a wide range of parameters. 
!
! The spherical Bessel function of the second kind y_p(x) can be computed by
! upward recurrence relation.
!
! The 0th and 1st order of these functions are
!
!            sin(x)              sin(x)       cos(x)
! j_0(x) =  -------- , j_1(x) = --------  -  --------
!              x                  x^2           x
!
!             cos(x)               cos(x)     sin(x)
! y_0(x) = - ------- , y_1(x) = - -------  - --------
!               x                   x^2         x
!
!--------------------------------------------------------------------------------
! isol = 1 : Series expansion 
!--------------------------------------------------------------------------------
!
! The series expansion (see Equation 2 of Jablonski 1994) can be rewritten 
! in the form (n>=1):
!
!                --       imax    --
!                |        ---      |
!  j_n(x) =  A_n |  1  +  \    B_i |
!                |        /--      |
!                --       i=1     --
!
!  where, imax determines the truncation order of the series expansion.
!  A_n and B_i are computed by using recurrence relations:
!
!           x                        - x^2
!  A_n = -------- A_{n-1},  B_i = ------------- B_{i-1}
!         2n + 1                  2i*(2i+2n+1)
!
!  with K_0 = 1.0 and A_0 = 1.0
!
!--------------------------------------------------------------------------------
! isol = 2:  Downward recurrence relation 
!--------------------------------------------------------------------------------
! 
!  First of all, we accumulate 'pseudo' series k_n in downward direction:
!
!                       2n + 3
!  k_n(x) = -k_{n+2} + -------- k_{n+1}, k_L=0.0, K_L-1=1.0
!                          x
!
!  where L = l + n_warm-up, where n_warm-up=100 is warm-up. 
!  After that calibrate the pseudo series to be the same as the spherical Bessel
!  function:
!
!       j_n(x) = s * k_n(x) ( n = 1, ..., l )
!  
!  where s = j_0(x)/k_0(x).
!
!--------------------------------------------------------------------------------
! isol = 3 : Upward recurrence 
!--------------------------------------------------------------------------------
!
!                  2n + 1
!    z_{n+1}(x) = ---------  z_n(x) - z_{n-1}(x)
!                     x
!
!  where, z_n = j_n or y_n.
!
!--------------------------------------------------------------------------------
!
!  The spherical Bessel functions of higher order and small argument 
!  can be very small/large number, which easily causes under/over flow. 
!  Thus, I adopt floor and ceilling values for j_p(x) and y_p(x).
!
!--------------------------------------------------------------------------------
subroutine sphbessel(M,x,SJ,SY,isol)
use types
implicit none
integer::M,n,i,isol
integer,parameter        :: imax = 100    ! Truncation order of series expansion
integer,parameter        :: nwarmup = 100 ! Warming-up 
real(kind=dp),parameter  :: eps=1.0e-30_dp
!real(kind=dp),parameter  :: floorvalue=1.0e-140_dp
!real(kind=dp),parameter  :: ceilingvalue=1.0e140_dp
real(kind=dp),parameter  :: floorvalue=1.0e-70_dp
real(kind=dp),parameter  :: ceilingvalue=1.0e70_dp
real(kind=dp)            :: K0,K1,x,s,Y0,Y1
real(kind=dp)            :: wa,xi,FN
real(kind=dp),dimension(0:M)::SJ,SY
real(kind=dp),dimension(0:M+nwarmup)::K
!
! The spherical Bessel function of the first kind j_p(x)
!
SJ    = 0.0_dp
SJ(0) = sin(x)/x
SY=0.0_dp
SY(0) = -cos(x)/x

if(M .eq. 0) return

if(isol .eq. 1) then
        ! series expansion
        SJ(0) = sin(x)/x
        FN = 1.0_dp
        do n=1,M
        FN = x / real(2*n+1,kind=dp) * FN
        xi = 1.0_dp
        wa = 0.0_dp
        do i=1,imax
                xi = - x * x * xi / real(2*i*(2*i+2*n+1),kind=dp)
                wa = wa + xi
                if(abs(xi/wa) .le. floorvalue) exit
        enddo
        SJ(n) = FN * (1.0_dp + wa) 
        if(abs(SJ(n)) .le. floorvalue) exit
        enddo

elseif(isol .eq. 2) then
        ! downward recurrence
        K1 = 0.0_dp
        K0 = 1.0_dp
        K  = 0.0_dp
        SJ = 0.0_dp
        SJ(0) = sin(x)/x
        do n=M+nwarmup,0,-1
                K(n) = -K1 + (real(2*n+3,kind=dp)/x)*K0
                K1 = K0
                K0 = K(n)
        enddo
        s = SJ(0)/K(0)
        do n=1,M
                SJ(n)=s*K(n)
                if(abs(SJ(n)) .le. floorvalue) exit
        enddo

elseif(isol .eq. 3) then
        ! upward recurrence
        SJ(0) = sin(x)/x
        SJ(1) = sin(x)/(x*x)-cos(x)/x
        K0 = SJ(0)
        K1 = SJ(1)
        do n=1,M-1
                SJ(n+1) = ((2.0_dp*real(n,kind=dp)+1)/x)*K1 - K0
                K0 = K1
                K1 = SJ(n+1)
        enddo
endif

!
! The spherical Bessel function of the second kind y_p(x)
!
SY(1) = -cos(x)/(x*x) - sin(x)/x 
Y0 = SY(0)
Y1 = SY(1)
do n=2,M
        SY(n) = ((2.0_dp*real(n-1,kind=dp)+1.0_dp)/x)*Y1-Y0
        Y0 = Y1
        Y1 = SY(n)
        if(ABS(SY(n)) .ge. ceilingvalue) exit
enddo

return
end subroutine sphbessel

!--------------------------------------------------------------------------------
!
! Output structure integration s_p(kRg)
!
!--------------------------------------------------------------------------------
subroutine spout(iqcor,df)
use types; use const
implicit none
integer          :: iqcor,p,imax,i
real(kind=dp)    :: df,k,xg,xmin,xmax,dx
real(kind=dp)    :: ImS0_xlt1,ImS0_xgt1,ReSp
complex(kind=dp) :: sp
k    = 1.0_dp
xmin = 1.0e-8_dp
xmax = 1.0e8_dp
imax = 250
dx=(xmax/xmin)**(1.0_dp/real(imax-1,kind=dp))
write(*,*) "Writing Sp ..."
open(21,file="sp.out",status="unknown")
write(21,6600) "order p","k*R_g","Re(Sp)","Im(Sp)"
do p=0,20
        write(*,*) "p = ",p
        do i=1,imax
        xg = xmin * dx ** real(i-1,kind=dp)
        call integration_of_Sp(iqcor,df,p,xg,Sp)
        ! asymptotic solution for Im(S_p) with p=0.
        ImS0_xlt1 = -sqrt(2.0_dp*pi)/4.0_dp/xg
        ImS0_xgt1 = -pi/8.0_dp/xg/xg
        ! asymptotic solution for Re(s_p) for all p.
        ! only valid for 2 < df <= 3.
        ReSp = (df/(16.0*xg*xg))&
                *Gamma(0.5_dp+0.5_dp*(df-3.0_dp))/Gamma(0.5_dp*df)
        write(21,6000) p,xg,real(Sp),aimag(Sp),ReSp
        enddo
        write(21,*)
enddo
close(21)

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

write(*,*) "sp.out is produced"
6000 format(' ',I15,1P8E15.5)
6600 format(' ',8A15)
return
end subroutine spout

!*****************************************************************************80
!
!! CHGM computes the confluent hypergeometric function M(a,b,x).
!
!  Licensing:
!
!    This routine is copyrighted by Shanjie Zhang and Jianming Jin.  However, 
!    they give permission to incorporate this routine into a user program 
!    provided that the copyright is acknowledged.
!
!  Modified:
!
!    27 July 2012
!
!  Author:
!
!    Shanjie Zhang, Jianming Jin
!
!  Reference:
!
!    Shanjie Zhang, Jianming Jin,
!    Computation of Special Functions,
!    Wiley, 1996,
!    ISBN: 0-471-11963-6,
!    LC: QA351.C45.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) A, B, parameters.
!
!    Input, real ( kind = 8 ) X, the argument.
!
!    Output, real ( kind = 8 ) HG, the value of M(a,b,x).
!
!--------------------------------------------------------------------------------
!
!   Note by R.T.
!   .f90 version was downloaded from 
!   https://people.sc.fsu.edu/~jburkardt/f_src/special_functions/special_functions.html
!   
!   Revision History 
!       2020. Oct 31:   - Minor change so that a real type is specified by "dp".
!                       - D+00 --> _dp
!
!--------------------------------------------------------------------------------
  subroutine chgm ( a, b, x, hg )
  use types; use const
  implicit none
  real ( kind = dp ) a
  real ( kind = dp ) a0
  real ( kind = dp ) a1
  real ( kind = dp ) b
  real ( kind = dp ) hg
  real ( kind = dp ) hg1
  real ( kind = dp ) hg2
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  integer ( kind = 4 ) la
  integer ( kind = 4 ) m
  integer ( kind = 4 ) n
  integer ( kind = 4 ) nl
!  real ( kind = dp ) pi
  real ( kind = dp ) r
  real ( kind = dp ) r1
  real ( kind = dp ) r2
  real ( kind = dp ) rg
  real ( kind = dp ) sum1
  real ( kind = dp ) sum2
  real ( kind = dp ) ta
  real ( kind = dp ) tb
  real ( kind = dp ) tba
  real ( kind = dp ) x
  real ( kind = dp ) x0
  real ( kind = dp ) xg
  real ( kind = dp ) y0
  real ( kind = dp ) y1

  !pi = 3.141592653589793D+00
  a0 = a
  a1 = a
  x0 = x
  hg = 0.0_dp
  !S: RT added initialization of y1
  y1 = 0.0_dp
  !E: RT

  if ( b == 0.0_dp .or. b == - abs ( int ( b ) ) ) then
    hg = 1.0D+300
  else if ( a == 0.0_dp .or. x == 0.0_dp ) then
    hg = 1.0_dp
  else if ( a == -1.0_dp ) then
    hg = 1.0_dp - x / b
  else if ( a == b ) then
    hg = exp ( x )
  else if ( a - b == 1.0_dp ) then
    hg = ( 1.0_dp + x / b ) * exp ( x )
  else if ( a == 1.0_dp .and. b == 2.0_dp ) then
    hg = ( exp ( x ) - 1.0_dp ) / x
  else if ( a == int ( a ) .and. a < 0.0_dp ) then
    m = int ( - a )
    r = 1.0_dp
    hg = 1.0_dp
    do k = 1, m
      r = r * ( a + k - 1.0_dp ) / k / ( b + k - 1.0_dp ) * x
      hg = hg + r
    end do
  end if

  if ( hg /= 0.0_dp ) then
    return
  end if

  if ( x < 0.0_dp ) then
    a = b - a
    a0 = a
    x = abs ( x )
  end if

  if ( a < 2.0_dp ) then
    nl = 0
  end if

  if ( 2.0_dp <= a ) then
    nl = 1
    la = int ( a )
    a = a - la - 1.0_dp
  end if

  do n = 0, nl

    if ( 2.0_dp <= a0 ) then
      a = a + 1.0_dp
    end if

    if ( x <= 30.0_dp + abs ( b ) .or. a < 0.0_dp ) then

      hg = 1.0_dp
      rg = 1.0_dp
      do j = 1, 500
        rg = rg * ( a + j - 1.0_dp ) &
          / ( j * ( b + j - 1.0_dp ) ) * x
        hg = hg + rg
        if ( abs ( rg / hg ) < 1.0D-15 ) then
          exit
        end if
      end do

    else

      call gamma ( a, ta )
      call gamma ( b, tb )
      xg = b - a
      call gamma ( xg, tba )
      sum1 = 1.0_dp
      sum2 = 1.0_dp
      r1 = 1.0_dp
      r2 = 1.0_dp
      do i = 1, 8
        r1 = - r1 * ( a + i - 1.0_dp ) * ( a - b + i ) / ( x * i )
        r2 = - r2 * ( b - a + i - 1.0_dp ) * ( a - i ) / ( x * i )
        sum1 = sum1 + r1
        sum2 = sum2 + r2
      end do
      hg1 = tb / tba * x ** ( - a ) * cos ( pi * a ) * sum1
      hg2 = tb / ta * exp ( x ) * x ** ( a - b ) * sum2
      hg = hg1 + hg2

    end if

    if ( n == 0 ) then
      y0 = hg
    else if ( n == 1 ) then
      y1 = hg
    end if

  end do

  if ( 2.0_dp <= a0 ) then
    do i = 1, la - 1
      hg = ( ( 2.0_dp * a - b + x ) * y1 + ( b - a ) * y0 ) / a
      y0 = y1
      y1 = hg
      a = a + 1.0_dp
    end do
  end if

  if ( x0 < 0.0_dp ) then
    hg = hg * exp ( x0 )
  end if

  a = a1
  x = x0

  return
end subroutine chgm


!*****************************************************************************80
!
!  GAMMA evaluates the Gamma function.
!
!  Licensing:
!
!    The original FORTRAN77 version of this routine is copyrighted by 
!    Shanjie Zhang and Jianming Jin.  However, they give permission to 
!    incorporate this routine into a user program that the copyright 
!    is acknowledged.
!
!  Modified:
!
!    08 September 2007
!
!  Author:
!
!    Original FORTRAN77 version by Shanjie Zhang, Jianming Jin.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Shanjie Zhang, Jianming Jin,
!    Computation of Special Functions,
!    Wiley, 1996,
!    ISBN: 0-471-11963-6,
!    LC: QA351.C45
!
!  Parameters:
!
!    Input, real ( kind = 8 ) X, the argument.
!    X must not be 0, or any negative integer.
!
!    Output, real ( kind = 8 ) GA, the value of the Gamma function.
!
!--------------------------------------------------------------------------------
!
!   Note by R.T.
!   .f90 version was downloaded from 
!   https://people.sc.fsu.edu/~jburkardt/f_src/special_functions/special_functions.html
!   
!   Revision History 
!       2020. Oct 31:   - Minor change so that a real type is specified by "dp".
!                       - D+00 --> _dp
!
!--------------------------------------------------------------------------------

  subroutine gamma ( x, ga )
  use types; use const
  implicit none
  real ( kind = dp ), dimension ( 26 ) :: g = (/ &
    1.0_dp, &
    0.5772156649015329_dp, &
   -0.6558780715202538_dp, &
   -0.420026350340952e-01_dp, &
    0.1665386113822915_dp, &
   -0.421977345555443e-01_dp, &
   -0.96219715278770e-02_dp, &
    0.72189432466630e-02_dp, &
   -0.11651675918591e-02_dp, &
   -0.2152416741149e-03_dp, &
    0.1280502823882e-03_dp, & 
   -0.201348547807e-04_dp, &
   -0.12504934821e-05_dp, &
    0.11330272320e-05_dp, &
   -0.2056338417e-06_dp, & 
    0.61160950e-08_dp, &
    0.50020075e-08_dp, &
   -0.11812746e-08_dp, &
    0.1043427e-09_dp, & 
    0.77823e-11_dp, &
   -0.36968e-11_dp, &
    0.51e-12_dp, &
   -0.206e-13_dp, &
   -0.54e-14_dp, &
    0.14e-14_dp, &
    0.1e-15_dp /)
  real ( kind = dp ) ga
  real ( kind = dp ) gr
  integer ( kind = 4 ) k
  integer ( kind = 4 ) m
  integer ( kind = 4 ) m1
  !real ( kind = dp ), parameter :: pi = 3.141592653589793D+00
  real ( kind = dp ) r
  real ( kind = dp ) x
  real ( kind = dp ) z

  if ( x == aint ( x ) ) then

    if ( 0.0_dp < x ) then
      ga = 1.0_dp
      m1 = int ( x ) - 1
      do k = 2, m1
        ga = ga * k
      end do
    else
      ga = 1.0D+300
    end if

  else

    if ( 1.0_dp < abs ( x ) ) then
      z = abs ( x )
      m = int ( z )
      r = 1.0_dp
      do k = 1, m
        r = r * ( z - real ( k, kind = dp ) )
      end do
      z = z - real ( m, kind = dp )
    else
      z = x
    end if

    gr = g(26)
    do k = 25, 1, -1
      gr = gr * z + g(k)
    end do

    ga = 1.0_dp / ( gr * z )

    if ( 1.0_dp < abs ( x ) ) then
      ga = ga * r
      if ( x < 0.0_dp ) then
        ga = - pi / ( x* ga * sin ( pi * x ) )
      end if
    end if

  end if

  return
end subroutine gamma

!--------------------------------------------------------------------------------
!
! This subroutine computes the geometrical cross section of an aggregate
! based on Equation (47) of Okuzumi et al. (2009).
!
! In Tazaki & Tanaka (2018), the geometrical cross section is computed by 
! using the fitting formulae for BCCA and BPCA giben by Minato et al. (2006).
! 
! However, this code is still applicable for an aggregate having arbitrary
! fractal dimension. To find the geometrical cross section of such aggregates,
! I use Equation (47) in Okuzumi et al. (2009), where they found the empirical
! formula by interpolating the Minato et al.'s formula for such cases.
!
!--------------------------------------------------------------------------------
subroutine geocross(PN,a0,ac,GC)
use types; use const
implicit none
real(kind=dp),parameter :: k0_bcca=1.04_dp ! Tazaki+2016; BCCA with oblique collisions
real(kind=dp),parameter :: df_bcca=1.90_dp ! Tazaki+2016; BCCA with oblique collisions
real(kind=dp) :: PN,a0,ac,ac_bcca,GC0,GC_BCCA,GC_BCCA2,GCC,GC
ac_bcca = sqrt(5.0_dp/3.0_dp) * a0 * (PN/k0_bcca)**(1.0_dp/df_bcca)
GC0  = pi * a0 * a0
GCC  = pi * ac * ac
GC_BCCA2  = pi * ac_bcca * ac_bcca
if(PN .le. 16.d0) then
      GC_BCCA = 12.5_dp * PN ** 0.685_dp * exp(-2.53_dp/PN**0.0920_dp)
else
      GC_BCCA = 0.352_dp * PN + 0.566_dp * PN ** 0.862_dp 
endif
GC_BCCA = GC0 * GC_BCCA
!
GC = (1.0_dp/GC_BCCA+1.0_dp/GCC-1.0_dp/GC_BCCA2)
GC = 1.0_dp / GC
! If GC > N*pi*a0^2, then GC=N*pi*a0^2
GC = min(GC,PN*pi*a0*a0)
!
return
end subroutine geocross

!--------------------------------------------------------------------------------
!
! Fitting formulae of geometrical cross section derived from 
! the mean-field extinction cross section at short-wavelength limit.
!
! CAUTION:
! The fractal prefactor is assumed to have a value obtained by
! linear interpolation (or extrapolation) of those values of 
! BCCA and BPCA.
!
!--------------------------------------------------------------------------------
subroutine geocross_mft(PN,df,Gratio)
use types
implicit none
integer::i
integer,parameter::ndim=15
real(kind=dp),dimension(0:ndim)::a,b,c,d,e,f,g,h,fdim
real(kind=dp)                ::aa,bb,cc,dd,ee,ff,gg,hh
real(kind=dp)                ::PN,df,Gratio
DATA (fdim(i),i=0,15) / 1.50000E+00,    1.60000E+00,    1.70000E+00,    1.80000E+00,&
                        1.90000E+00,    2.00000E+00,    2.10000E+00,    2.20000E+00,&
                        2.30000E+00,    2.40000E+00,    2.50000E+00,    2.60000E+00,&
                        2.70000E+00,    2.80000E+00,    2.90000E+00,    3.00000E+00 /
DATA (a(i),i=0,15) /    7.13973E+02,    1.16423E+03,    3.33395E+02,    3.56949E+02,&
                        2.74600E+00,    3.04834E+02,    3.22197E+02,    4.03725E+02,&
                        4.39512E+02,    4.26389E+02,    4.53031E+02,    4.52352E+02,&
                        4.43447E+02,    4.35774E+02,    4.40590E+02,    4.50909E+02 /
DATA (b(i),i=0,15) /   -1.08072E-01,   -1.05630E-01,   -1.15695E-01,   -1.24080E-01,&
                       -2.52337E-03,   -8.12464E-02,   -3.39455E-02,    7.18452E-02,&
                        6.30843E-02,   -1.98444E-02,   -1.24053E-01,   -1.74670E-01,&
                       -1.95718E-01,   -2.04959E-01,   -2.07496E-01,   -2.04854E-01 /
DATA (c(i),i=0,15) /    6.96352E+00,    7.44311E+00,    6.16310E+00,    6.16324E+00,&
                        1.18113E+00,    5.90229E+00,    5.93429E+00,    6.13033E+00,&
                        6.21396E+00,    6.37749E+00,    6.34347E+00,    6.30390E+00,&
                        6.25537E+00,    6.21575E+00,    6.20834E+00,    6.21497E+00 /
DATA (d(i),i=0,15) /    2.89145E-02,    2.28560E-02,    2.83661E-02,    3.18105E-02,&
                        9.26620E-02,    2.80617E-02,    1.93156E-02,    2.41211E-05,&
                       -9.34488E-05,   -8.23762E-05,    1.91956E-02,    2.93221E-02,&
                        3.44086E-02,    3.71157E-02,    3.81358E-02,    3.78502E-02 /
DATA (e(i),i=0,15) /    1.52260E+00,    1.29338E+00,    1.10944E+00,    9.16614E-01,&
                        6.77090E-01,    6.26109E-01,    5.38993E-01,    4.52138E-01,&
                        3.94309E-01,    6.19632E-01,    4.83615E-01,    4.00043E-01,&
                        3.32825E-01,    2.78355E-01,    2.32785E-01,    1.93256E-01 /
DATA (f(i),i=0,15) /   -5.05108E-02,   -1.51396E-02,    9.11257E-03,    2.25686E-02,&
                        1.34713E-02,    5.34920E-02,    7.53823E-02,    1.08882E-01,&
                        1.53220E-01,    1.65246E-01,    1.98459E-01,    2.27960E-01,&
                        2.55767E-01,    2.81834E-01,    3.06280E-01,    3.29231E-01 /
DATA (g(i),i=0,15) /    1.55847E+00,    1.41197E+00,    1.32293E+00,    1.31261E+00,&
                        1.46648E+00,    1.32681E+00,    1.30102E+00,    1.31291E+00,&
                        1.19121E+00,    8.13122E-01,    8.65730E-01,    8.40165E-01,&
                        8.05872E-01,    7.64332E-01,    7.18741E-01,    6.72135E-01 /
DATA (h(i),i=0,15) /    8.90935E-01,    8.75141E-01,    8.65479E-01,    8.80248E-01,&
                        9.46148E-01,    9.44240E-01,    9.74293E-01,    1.04852E+00,&
                        1.11467E+00,    8.78236E-01,    9.61048E-01,    1.00913E+00,&
                        1.05517E+00,    1.10051E+00,    1.14855E+00,    1.20307E+00 /
do i=0,ndim-1
      if(fdim(i) .le. df .and. fdim(i+1) .ge. df) then
              aa = (a(i+1)-a(i))/(fdim(i+1)-fdim(i))*(df-fdim(i))+a(i)
              bb = (b(i+1)-b(i))/(fdim(i+1)-fdim(i))*(df-fdim(i))+b(i)
              cc = (c(i+1)-c(i))/(fdim(i+1)-fdim(i))*(df-fdim(i))+c(i)
              dd = (d(i+1)-d(i))/(fdim(i+1)-fdim(i))*(df-fdim(i))+d(i)
              ee = (e(i+1)-e(i))/(fdim(i+1)-fdim(i))*(df-fdim(i))+e(i)
              ff = (f(i+1)-f(i))/(fdim(i+1)-fdim(i))*(df-fdim(i))+f(i)
              gg = (g(i+1)-g(i))/(fdim(i+1)-fdim(i))*(df-fdim(i))+g(i)
              hh = (h(i+1)-h(i))/(fdim(i+1)-fdim(i))*(df-fdim(i))+h(i)
              exit
      endif
enddo
Gratio  = aa*PN**bb*exp(-cc/PN**dd)+ee*PN**ff*exp(-gg/PN**hh)
Gratio  = 1.0_dp / Gratio
return
end subroutine geocross_mft

subroutine opticslimit(iqcor,lmd,refrel,df,k0,PN,R0,Cext)
use types; use const
implicit none
integer::iqcor,i,nstop,n
real(kind=dp)::lmd,R0,Rg,x0,sp,xstop,df,k0,PN,k,xg,Cext,sumA,sumB
complex(kind=dp),allocatable,dimension(:)::an,bn
complex(kind=dp),allocatable,dimension(:,:)::dd
complex(kind=dp)::refrel

if(iqcor .ne. 1) then
        print *, 'ERROR'
        stop
endif
k       =  twopi/lmd                   ! wavenumber
!k       =  twopi/0.1_dp                 ! wavenumber
Rg      =  R0 * (PN/k0) ** (1.0_dp/df)  ! Radius of gyration of the aggregate
xg      =  k * Rg                       ! size parameter of aggregates
x0      =  k * R0                       ! size parameter of a monomer
!refrel = cmplx(2.0_dp,1.0_dp,kind=dp)
sp = (df/(16.0*xg*xg))&
    *Gamma(0.5_dp*(df-2.0_dp))/Gamma(0.5_dp*df)
xstop = x0 + 4.0_dp * x0 ** (1.0_dp/3.0_dp) + 2.0_dp
nstop = nint(xstop)
allocate(an(nstop),bn(nstop),dd(2,nstop))
call lorenz_mie(x0,refrel,nstop,an,bn)

! difflimit
!an=0.5_dp
!bn=0.5_dp

sumA=0.0_dp
sumB=0.0_dp
do n=1,nstop
        sumA = sumA + real(2*n+1,kind=dp)*an(n)
        sumB = sumB + real(2*n+1,kind=dp)*bn(n)
enddo
do n=1,nstop
        dd(1,n) = an(n)/(1.0_dp+(PN-1.0_dp)*Sp*sumA)
        dd(2,n) = bn(n)/(1.0_dp+(PN-1.0_dp)*Sp*sumB)
enddo
Cext = 0.0_dp
do n=1,nstop
        Cext = Cext + real(2*n+1,kind=dp)*real(dd(1,n)+dd(2,n))
enddo
Cext  = Cext * twopi * PN / k / k
deallocate(an,bn,dd)
return
end subroutine opticslimit


!--------------------------------------------------------------------------------
!
! ludcmp      : Press, W. H. et al. (1997), "Numerical Recipes in Fortran 77"
!
!--------------------------------------------------------------------------------
!
!   Note by R.T.
!       Original code is written in F77.
!
!   Revision History 
!       2020. Oct 31:   
!              - Avoid using old grammers; PAUSE, CONTINUE...
!              - Minor change so that a real type is specified by "dp".
!              - D+00 --> _dp
!--------------------------------------------------------------------------------
SUBROUTINE ludcmp(a,n,np,indx,d)
use types
implicit none
INTEGER n,np,indx(n),NMAX
real(kind=dp)::d,a(np,np),TINY
PARAMETER (NMAX=1000,TINY=1.0e-20)
INTEGER::i,imax,j,k
real(kind=dp)::aamax,dum,sum,vv(NMAX)
d=1.
imax=1
do i=1,n
aamax=0.
do j=1,n
  if (abs(a(i,j)).gt.aamax) aamax=abs(a(i,j))
enddo
!if (aamax.eq.0.) pause 'singular matrix in ludcmp' ! pause is now deleted feature.
!if (aamax.eq.0.) 'singular matrix' stop
if (aamax.eq.0.) then
        write(*,*) 'singular matrix in ludcmp. STOP'
        stop
endif
vv(i)=1./aamax
enddo
do j=1,n
do i=1,j-1
  sum=a(i,j)
  do k=1,i-1
    sum=sum-a(i,k)*a(k,j)
  enddo
  a(i,j)=sum
enddo
aamax=0.
do i=j,n
  sum=a(i,j)
  do k=1,j-1
    sum=sum-a(i,k)*a(k,j)
  enddo
  a(i,j)=sum
  dum=vv(i)*abs(sum)
  if (dum.ge.aamax) then
    imax=i
    aamax=dum
  endif
enddo
if (j.ne.imax)then
  do k=1,n
    dum=a(imax,k)
    a(imax,k)=a(j,k)
    a(j,k)=dum
  enddo
  d=-d
  vv(imax)=vv(j)
endif
indx(j)=imax
if(a(j,j).eq.0.)a(j,j)=TINY
if(j.ne.n)then
  dum=1./a(j,j)
  do i=j+1,n
    a(i,j)=a(i,j)*dum
  enddo
endif
enddo
return
END SUBROUTINE ludcmp
!--------------------------------------------------------------------------------
!
! lubksb      : Press, W. H. et al. (1997), "Numerical Recipes in Fortran 77"
!
!--------------------------------------------------------------------------------
!
!   Note by R.T.
!       Original code is written in F77.
!
!   Revision History 
!       2020. Oct 31:   
!              - Avoid using old grammers; PAUSE, CONTINUE...
!              - Minor change so that a real type is specified by "dp".
!              - D+00 --> _dp
!--------------------------------------------------------------------------------

SUBROUTINE lubksb(a,n,np,indx,b)
use types
implicit none
INTEGER n,np,indx(n)
!REAL a(np,np),b(n)
real(kind=dp)::a(np,np),b(n)
INTEGER i,ii,j,ll
!REAL sum
real(kind=dp)::sum
ii=0
do i=1,n
ll=indx(i)
sum=b(ll)
b(ll)=b(i)
if (ii.ne.0)then
  do j=ii,i-1
    sum=sum-a(i,j)*b(j)
  enddo
else if (sum.ne.0.) then
  ii=i
endif
b(i)=sum
enddo
do i=n,1,-1
sum=b(i)
do j=i+1,n
  sum=sum-a(i,j)*b(j)
enddo
b(i)=sum/a(i,i)
enddo
return
END SUBROUTINE lubksb
!--------------------------------------------------------------------------------
!
! gauleg      : Press, W. H. et al. (1997), "Numerical Recipes in Fortran 77"
!
!--------------------------------------------------------------------------------
!
!   Note by R.T.
!       Original code is written in F77.
!
!   Revision History 
!       2020. Oct 31:   
!              - Avoid using old grammers; PAUSE, CONTINUE...
!              - Minor change so that a real type is specified by "dp".
!              - D+00 --> _dp
!              - EPS value is changed from 3.e-14 to 1.e-14 according to 
!                the precision change (Numerical Recipes, Appendix C1, p.1363).
!
!--------------------------------------------------------------------------------
SUBROUTINE gauleg(x1,x2,x,w,n)
use types
implicit none
INTEGER n
real(kind=dp)::x1,x2,x(n),w(n)
real(kind=dp)::EPS
PARAMETER (EPS=1.0e-14_dp) !as suggested by NumericalRecipeF90
INTEGER i,j,m
real(kind=dp)::p1,p2,p3,pp,xl,xm,z,z1
!S: RT added initialization of pp 
! to suppress warning from [-Wmaybe-uninitialized]
pp = 0.0_dp
!E: RT
m=(n+1)/2
xm=0.5d0*(x2+x1)
xl=0.5d0*(x2-x1)
do i=1,m
z=cos(3.141592654d0*(i-.25d0)/(n+.5d0))
!S: RT 
! Remove goto sentence, and add do while.
! To enter the first do-while loop, 
! I put arbitrary some z1, which should differ 
! by more than EPS.
Z1=2.0*z 
!E: RT
do while(abs(z-z1).gt.EPS)
  p1=1.d0
  p2=0.d0
  do j=1,n
    p3=p2
    p2=p1
    p1=((2.d0*j-1.d0)*z*p2-(j-1.d0)*p3)/j
  enddo
  pp=n*(z*p1-p2)/(z*z-1.d0)
  z1=z
  z=z1-p1/pp
enddo
x(i)=xm-xl*z
x(n+1-i)=xm+xl*z
w(i)=2.d0*xl/((1.d0-z*z)*pp*pp)
w(n+1-i)=w(i)
enddo
return
END SUBROUTINE gauleg
