!--------------------------------------------------------------------------------
!
!       GEOFRACTAL version 1.0
!
!--------------------------------------------------------------------------------
!
! GEOFRACTAL v1.0 computes the geometrical cross section of randomly orientated 
! fractal dust aggregates based on a statistical distribution model of monomer 
! particles developed in Tazaki (in prep.).
! 
! If you wish to publish a paper that contains results obtained by this code, 
! please acknowledge this code and cite relevant papers. 
!
!
! Disclaimer:
! I reject all responsibility for the use of this code. Although it has been
! tested, the code may still contain a bug. I am not responsible for any damages
! caused by the use of the code. If you find a bug, please let me know.
!
!                                                   Ryo Tazaki (2021/Jan/1)
!                                                   email: r.tazaki -at- uva.nl 
!
!--------------------------------------------------------------------------------
! INPUT PARAMETERS
!--------------------------------------------------------------------------------
!
! Df            : fractal dimension
! PN            : number of monomers
! k0            : fractal prefactor defined by PN = k0*(Rg/R0)^Df
!
! iqapp         : Select a method for solving radial and angular integration 
!                 when computing the mean overlapping efficiency.
!                --------------------------------------------------------
!                          |    radial  (x)        |  angular (u)        |
!                ---------------------------------------------------------
!                iqapp = 1 |   numerical           |  numerical          |  
!                iqapp = 3 |   numerical  ( D<=2 ) |  analytical         |
!                          |   analytical ( D> 2 ) |                     |
!                --------------------------------------------------------
!
! iqcon         : Switch to specify the numerical factor A.
!                 iqcon  = 1 : A = 1 
!                 iqcon  = 2 : A is chosen to connect two regimes.
!
! iqcor         : Switch to specify a cut-off model of the correlation function
!                 iqcon  = 1 : The Gaussian cut-off function 
!                 iqcon  = 2 : The exponential cut-off function
!                 iqcon  = 3 : The fractal dimension cut-off function
!
! Prefered set of options is iqcor=3, iqapp=3 and iqcon=2.
! Note that iqcor=1 is also possible, but avoid using iqcor=2.
!
!--------------------------------------------------------------------------------
! OUTPUT PARAMETER
!--------------------------------------------------------------------------------
!
! G  : The geocross section of the aggregate normalized by a factor PN*pi*R0^2. 
!      Thus, "G" is non-dimensional quantity. 
!
!--------------------------------------------------------------------------------

!module types
!implicit none
!integer, parameter      :: dp = selected_real_kind(P=15)
!end module types
!
!module const
!use types
!implicit none
!real(kind=dp),parameter :: pi = 3.1415926535897932384_dp
!end module const

subroutine geofrac(iqapp,iqcon,iqcor,PN,k0,df,G)
use types
implicit none
integer                 :: iqapp,iqcon,iqcor
real(kind=dp)           :: PN,k0,df,GS,G
real(kind=dp)           :: pnth,gth,sigmath,A,sigma
logical,parameter::debug=.false.

!--------------------------------------------------------------------------------
! some checks
!--------------------------------------------------------------------------------
if(PN .le. 0.9999_dp) then
        print *, 'number of monomers =',PN
        print *, 'error: The number of monomers should be larger than 1.'
        print *, 'stop'
        stop
endif
if(df .lt. 0.9999_dp .or. df .gt. 3.00001) then
        print *, 'fractal dimension df =',df
        print *, 'error: Fractal dimension df should be 1<=df<=3.'
        print *, 'stop'
        stop
endif
if(iqcon .ne. 1 .and. iqcon .ne. 2) then
        print *, 'error: iqcon should be either 1 or 2.'
        print *, 'stop'
        stop
endif
if(iqcor .ne. 1 .and. iqcor .ne. 2 .and. iqcor .ne. 3) then
        print *, 'error: iqcor should be either 1, 2, or 3.'
        print *, 'stop'
        stop
endif
if(iqapp .ne. 1 .and. iqapp .ne. 2 .and. iqapp .ne. 3) then
        print *, 'error: iqapp should be either 1 or 3.'
        print *, 'stop'
        stop
endif
if(iqapp .eq. 2 .and. (debug .eqv. .false.)) then
        print *, 'error: iqapp = 2 is only allowed for a debug mode.'
        print *, 'stop'
        stop
endif

!--------------------------------------------------------------------------------
! start calculations
!--------------------------------------------------------------------------------
A = 1.0_dp
if(iqcon .eq. 2) then
        PNTH = 11.0_dp*df-8.5_dp
        PNTH = min(PNTH,8.0_dp)
        if(PN .lt. PNTH) then
                G=GS(PN)
                return
        endif
        call mean_overlap_efficiency(iqapp,iqcor,k0,df,PNTH,sigmath)
        A = (1.0_dp+(PNTH-1.0_dp)*sigmath)*GS(PNTH)
endif
call mean_overlap_efficiency(iqapp,iqcor,k0,df,PN,sigma)
G = A / (1.0_dp + (PN-1.0_dp)*sigma)

return
end subroutine geofrac

!--------------------------------------------------------------------------------
!
! Minato et al. (2006)'s fitting formula for N < 16
!
!--------------------------------------------------------------------------------
function GS(PN)
use types
implicit none
real(kind=dp)::GS,PN
        GS = 12.5_dp*PN**(-0.315_dp)*exp(-2.53_dp/PN**0.0920_dp)
return
end function GS

!--------------------------------------------------------------------------------
!
! This subroutine calculations the mean overlapping efficiency: sigma
! (In Tazaki in prep., this quantity corresponds to \tilde{sigma} )
!
! For iqcor=1 [Gaussian cut-off model]
!
!               xmin       /ln(xmax)  
! sigma  = --------------  |   dlnx frad(x)*S(rho),
!          16*Gamma(df/2)  /ln(xmin)
!
!       xmin = df*(k0/PN)^(2/df); a=(df-2)/2; rho = sqrt(x/xmin)
!
! For iqcor=2 [Exponential cut-off model]
!
!              xmin^2      /ln(xmax)  
! sigma  = --------------  |   dlnx frad(x)*S(rho),
!           16*Gamma(df)   /ln(xmin)
!
!       xmin = 2*sqrt(df(df+1)/2)*(k0/PN)^(1.0/df); a=df-2; rho = x/xmin
!
! For iqcor=3 [Fractal dimension cut-off model]
!
!          xmin^{2/df} /ln(xmax)  
! sigma  = ----------  |  dlnx frad(x)*S(rho),
!             16      /ln(xmin)
!
!       xmin = 2^(df-1)*k0/PN; a=(df-2)/df; rho = (x/xmin)^(1/df)
!
! where df is fractal dimension, R0 is the monomer radius, Rg is 
! the radius of gyration, frad(x) is the radial integrand function, 
! and S(rho) is the angular integral function.
!
! frad(x)dlnx = x^{a}*exp(-x) dlnx = x^{a-1}exp(-x) dx
!
!           16*rho^2  / (1/rho)   
! S(rho)  = --------  |  du fang(rho,u),
!              pi     /0      
! 
! fang(rho,u) = [Arcsin{sqrt(1-rho^2u^2)} - rho*u*sqrt{1-rho^2u^2}]*u/sqrt(1-u^2)
!
! Note that the angular integral function S(rho) ≈1 when rho >> 1.
!
!--------------------------------------------------------------------------------
subroutine mean_overlap_efficiency(iqapp,iqcor,k0,df,PN,sigma)
use types; use const
implicit none
integer,parameter       :: nx = 1000
real(kind=dp),parameter :: xmax = 25.0_dp
integer                 :: i,iqapp,iqcor
real(kind=dp)           :: xmin,dlnx,stmp,k0,df,PN,sigma,frad
real(kind=dp)           :: factor,aicgm,gammq
real(kind=dp),allocatable::x(:),sang(:)

if(iqcor .eq. 1) then                           ! Gaussian
        aicgm  = 0.5_dp*(df-2.0_dp)
        xmin   = df * (k0/PN) ** (2.0/df)
        factor = xmin/16.0_dp/gamma(0.5_dp*df)
elseif(iqcor .eq. 2) then                       ! Exponential
        aicgm  = df - 2.0_dp
        xmin   = 2.0_dp * sqrt(0.5_dp*df*(df+1.0_dp))*(k0/PN)**(1.0/df)
        factor = xmin*xmin/16.0_dp/gamma(df)
elseif(iqcor .eq. 3) then                       ! Fractal dimension
        aicgm  = (df - 2.0_dp) / df
        xmin   = 2.0_dp**(df-1.0_dp) * k0 / PN
        factor = xmin ** (2.0_dp/df) / 16.0_dp
endif

if(df .gt. 2.0_dp .and. iqapp .eq. 3) then
        sigma  = factor*gamma(aicgm)*gammq(aicgm,xmin)
else
        dlnx = log(xmax/xmin)/real(nx-1,kind=dp)
        allocate(x(1:nx),sang(1:nx))
        do i=1,nx
                x(i) = exp(log(xmin)+real(i-1,kind=dp)*dlnx)
                if(iqapp .eq. 1) then
                        call angular_integration(iqcor,x(i),xmin,df,stmp)
                        sang(i) = stmp
                else
                        sang(i) = 1.0_dp
                endif
        enddo
        sigma = 0.0_dp
        do i=1,nx-1
        sigma = sigma + 0.5_dp * (frad(aicgm,x(i))*sang(i)+&
                frad(aicgm,x(i+1))*sang(i+1))
        enddo
        sigma = sigma * dlnx * factor
        deallocate(x,sang)
endif
return
end subroutine mean_overlap_efficiency

!--------------------------------------------------------------------------------
! The radial integrand function in ln(x) space
!--------------------------------------------------------------------------------
function frad(aicgm,x)
use types
implicit none
real(kind=dp)::aicgm,x,frad
        frad = x ** aicgm * exp(-x)
return
end function frad

!--------------------------------------------------------------------------------
! This subroutine computes the angular integral function S(rho) 
!--------------------------------------------------------------------------------
subroutine angular_integration(iqcor,x,xmin,df,sang)
use types; use const
implicit none
integer            :: i,j,iqcor
integer,parameter  :: nmax_u = 500
real(kind=dp)      :: x,xmin,df,rho,umin,umax,du,fang,sang
real(kind=dp)      :: u(1:nmax_u)

if(iqcor .eq. 1) then
        rho = sqrt(x/xmin)
elseif(iqcor .eq. 2) then
        rho = x/xmin
elseif(iqcor .eq. 3) then
        rho = (x/xmin) ** (1.0_dp/df)
endif

sang = 0.0_dp
umin = 0.0_dp
umax = 1.0_dp / rho
du   = (umax - umin) / real(nmax_u-1,kind=dp)
do j=1,nmax_u
        u(j) = umin + du * real(j-1,kind=dp)
enddo
do j=1,nmax_u-1
        sang = sang + 0.5_dp*(fang(rho,u(j))+fang(rho,u(j+1)))*du
enddo
sang = 16.0_dp * rho * rho * sang / pi
return
end subroutine angular_integration

!--------------------------------------------------------------------------------
!
! This function returns the integrand of the angular integral.
!  fang (rho,u) =  [ Arcsin(sqrt(1-rho^2u^2)) - & 
!               rho*u*sqrt(1-rho^2u^2)] * u / sqrt(1-u^2)
!
!--------------------------------------------------------------------------------
function fang(rho,u)
use types
implicit none
real(kind=dp)::rho,u,fang

if(abs(1.0_dp-rho*u) .le. 1.0e-10_dp) then
        fang = 0.0_dp
else
        fang = (asin(sqrt(1.0_dp-rho*rho*u*u))-&
                rho*u*sqrt(1.0_dp-rho*rho*u*u))*u/sqrt(1.0_dp-u*u)
endif

return
end function fang

!--------------------------------------------------------------------------------
!
! Following subroutines (gammq, gser, gcf, gammln) compute 
! the incomplete Gamma function Q(a,x)
!     
!             1     /infinity
! Q(a,x) = -------- | x^{a-1} exp(-x) dx
!          Gamma(a) /x
!
! where a>0. The definition of Q(a,x) obeys the definition in
!
! "Numerical Recipes in Fortran 77"
! Press, W. H., Teukolsky, S. A., Vetterling, W. T., et al.1992, 
! Cambridge: University Press, |c1992, 2nd ed
!
! Following subroutines are taken from Numerical Recippes in Fortran 77
!
!--------------------------------------------------------------------------------
!
! History (R.T.)
!
! The original subroutine was slightly modified so that the old F77 grammers
! are rewritten by modern grammers.
!
! Convergence accuracy eps is changed from 3e-7 to 3e-14.
!
!--------------------------------------------------------------------------------
function gammq(a,x)
use types
real(kind=dp) a,gammq,x
real(kind=dp) gammcf,gamser,gln
if(x.lt.0..or.a.le.0.) then
      write(*,*) 'error' 
      stop
endif
if(x.lt.a+1.)then
        call gser(gamser,a,x,gln)
        gammq=1.-gamser
else
        call gcf(gammcf,a,x,gln)
        gammq=gammcf
endif
return
end function gammq

subroutine gser(gamser,a,x,gln)
use types
integer itmax
real(kind=dp) a,gamser,gln,x,eps
parameter (itmax=100,eps=3.e-14)
integer n
real(kind=dp) ap,del,sum,gammln
gln=gammln(a)
if(x.le.0.)then
        if(x.lt.0.) then
                write(*,*) "error"
                stop 
        endif
        gamser=0.
        return
endif
ap=a
sum=1./a
del=sum
do n=1,itmax
        ap=ap+1.
        del=del*x/ap
        sum=sum+del
        if(abs(del).lt.abs(sum)*eps) exit
        if(n .eq. itmax) then
                print *, 'a too large, itmax too small in gser'
                stop
        endif
enddo
gamser=sum*exp(-x+a*log(x)-gln)
return
end subroutine gser
subroutine gcf(gammcf,a,x,gln)
use types
integer itmax
real(kind=dp) a,gammcf,gln,x,eps,fpmin
parameter (itmax=100,eps=3.e-14,fpmin=1.e-30)
integer i
real(kind=dp) an,b,c,d,del,h,gammln
gln=gammln(a)
b=x+1.-a
c=1./fpmin
d=1./b
h=d
do i=1,itmax
        an=-i*(i-a)
        b=b+2.
        d=an*d+b
        if(abs(d).lt.fpmin)d=fpmin
        c=b+an/c
        if(abs(c).lt.fpmin)c=fpmin
        d=1./d
        del=d*c
        h=h*del
        if(abs(del-1.).lt.eps) exit
        if(n .eq. itmax) then
                print *, 'a too large, itmax too small in gcf'
                stop
        endif
enddo
gammcf=exp(-x+a*log(x)-gln)*h
return
end subroutine gcf
function gammln(xx)
use types
real(kind=dp) gammln,xx
integer j
real(kind=dp) ser,stp,tmp,x,y,cof(6)
save cof,stp
data cof,stp/76.18009172947146d0,-86.50532032941677d0,&
        24.01409824083091d0,-1.231739572450155d0,.1208650973866179d-2,&
        -.5395239384953d-5,2.5066282746310005d0/
x=xx
y=x
tmp=x+5.5d0
tmp=(x+0.5d0)*log(tmp)-tmp
ser=1.000000000190015d0
do j=1,6
        y=y+1.d0
        ser=ser+cof(j)/y
enddo
gammln=tmp+log(stp*ser/x)
return
end function gammln
