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
use IOunits
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
        write(stde,*) 'number of monomers =',PN
        write(stde,*) 'error: The number of monomers should be larger than 1.'
        write(stde,*) 'stop'
        stop
endif
if(df .lt. 0.9999_dp .or. df .gt. 3.00001) then
        write(stde,*) 'fractal dimension df =',df
        write(stde,*) 'error: Fractal dimension df should be 1<=df<=3.'
        write(stde,*) 'stop'
        stop
endif
if(iqcon .ne. 1 .and. iqcon .ne. 2) then
        write(stde,*) 'error: iqcon should be either 1 or 2.'
        write(stde,*) 'stop'
        stop
endif
if(iqcor .ne. 1 .and. iqcor .ne. 2 .and. iqcor .ne. 3) then
        write(stde,*) 'error: iqcor should be either 1, 2, or 3.'
        write(stde,*) 'stop'
        stop
endif
if(iqapp .ne. 1 .and. iqapp .ne. 2 .and. iqapp .ne. 3) then
        write(stde,*) 'error: iqapp should be either 1 or 3.'
        write(stde,*) 'stop'
        stop
endif
if(iqapp .eq. 2 .and. (debug .eqv. .false.)) then
        write(stde,*) 'error: iqapp = 2 is only allowed for a debug mode.'
        write(stde,*) 'stop'
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
real(kind=dp)           :: factor,aicgm,upper_incomplete_gamma
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
        sigma  = factor*upper_incomplete_gamma(aicgm,xmin)
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
! Upper incomplete gamma function
!
!                  /infinity
!   Gamma(a,x) =   |  t^{a-1} exp(-t) dt,       a > 0, x >= 0
!                  /x
!
! For x < a+1, use Gamma(a,x) = Gamma(a) - g(a,x) with the power series
! of the lower incomplete gamma function (Abramowitz & Stegun 6.5.29)
!
!   g(a,x) = x^a exp(-x) sum_{k=0}^{infinity} x^k / (a (a+1) ... (a+k))
!
! For x >= a+1, use the continued fraction of Legendre (A&S 6.5.31),
! in its even form
!
!   Gamma(a,x) = x^a exp(-x) / D,
!   D = b_0 + a_1/(b_1 + a_2/(b_2 + ...)),  b_k = x+2k+1-a,  a_k = -k(k-a)
!
! D is evaluated with the forward recurrence for the convergents
! A_k/B_k, A_k = b_k A_{k-1} + a_k A_{k-2} (and the same for B_k),
! rescaled in each step to avoid overflow.
!
!--------------------------------------------------------------------------------
function upper_incomplete_gamma(a,x)
use types; use IOunits
implicit none
integer,parameter       :: kmax = 1000
real(kind=dp)           :: a,x,upper_incomplete_gamma
real(kind=dp)           :: eps,term,series,ak,bk,A0,A1,B0,B1,A2,B2,D,Dold
integer                 :: k

eps = 4.0_dp*epsilon(x)
if(a .le. 0.0_dp .or. x .lt. 0.0_dp) then
        write(stde,*) 'upper_incomplete_gamma: need a>0 and x>=0. STOP'
        stop
endif
if(x .eq. 0.0_dp) then
        upper_incomplete_gamma = gamma(a)
        return
endif

if(x .lt. a+1.0_dp) then
        ! Power series for the lower incomplete gamma function
        term = 1.0_dp/a
        series = term
        do k=1,kmax
                term = term*x/(a+k)
                series = series+term
                if(term .lt. eps*series) exit
        enddo
        if(k .gt. kmax) then
                write(stde,*) 'upper_incomplete_gamma: series did not converge. STOP'
                stop
        endif
        upper_incomplete_gamma = gamma(a) - exp(a*log(x)-x)*series
else
        ! Continued fraction
        A0 = 1.0_dp                  ! A_{-1}
        B0 = 0.0_dp                  ! B_{-1}
        A1 = x+1.0_dp-a              ! A_0 = b_0
        B1 = 1.0_dp                  ! B_0
        D  = A1
        do k=1,kmax
                ak = -k*(k-a)
                bk = x+2*k+1.0_dp-a
                A2 = bk*A1 + ak*A0
                B2 = bk*B1 + ak*B0
                ! Rescale so that B_k = 1
                A0 = A1/B2
                B0 = B1/B2
                A1 = A2/B2
                B1 = 1.0_dp
                Dold = D
                D    = A1
                if(abs(D-Dold) .le. eps*abs(D)) exit
        enddo
        if(k .gt. kmax) then
                write(stde,*) 'upper_incomplete_gamma: continued fraction did not converge. STOP'
                stop
        endif
        upper_incomplete_gamma = exp(a*log(x)-x)/D
endif
return
end function upper_incomplete_gamma
