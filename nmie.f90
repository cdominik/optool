!==============================================================================
! nmie.f90  -  Fortran 90 port of scattnlay nMie for multilayered spheres
!
! Implements the recursive algorithm of:
!   W. Yang, Appl. Opt. 42, 1710 (2003)
!
! Equations referenced as (N) refer to:
!   O. Pena & U. Pal, Comput. Phys. Commun. 180, 2348 (2009)
!
! Original C++ authors:
!   Ovidio Pena <ovidio@bytesfall.com>  (2009-2022)
!   Konstantin Ladutenko <kostyfisik@gmail.com>  (2013-2022)
!
! License: GNU GPL v3+
!==============================================================================

module nmie_mod
  implicit none
  private

  integer,  parameter :: dp = kind(1.0d0)
  real(dp), parameter :: PI = 3.14159265358979323846d0

  public :: nMie

contains

  !----------------------------------------------------------------------------
  ! Numerically stable complex cotangent.
  ! When Im(z) < 0, uses cot(z) = conj(cot(conj(z))).
  ! See H. Du, Appl. Opt. 43, 1951 (2004), Eqs. 10-12.
  !----------------------------------------------------------------------------
  function complex_cot(z) result(res)
    complex(dp), intent(in) :: z
    complex(dp)             :: res
    real(dp) :: re, im, sgn, ev, tv, a, b, c, d, denom

    re = real(z,  dp)
    im = aimag(z)
    sgn = sign(1.0d0, im)
    if (im == 0.0d0) sgn = 1.0d0

    ev = exp(-2.0d0 * sgn * im)   ! exp(-2*|Im(z)|)
    tv = tan(re)

    a = tv - ev*tv               !  tan - exp*tan
    b = 1.0d0 + ev
    c = -1.0d0 + ev
    d = tv + ev*tv               !  tan + exp*tan
    denom = c*c + d*d

    res = cmplx((a*c + b*d)/denom,  sgn*(b*c - a*d)/denom,  dp)
  end function complex_cot

  !----------------------------------------------------------------------------
  ! Conservative starting index for the downward D1 recurrence.
  ! Ensures D1[nstar]=0 is a negligible seed for the backward recurrence
  ! down to D1[0..nmax].
  !----------------------------------------------------------------------------
  function calc_nstar(nmax, z) result(nstar)
    integer,     intent(in) :: nmax
    complex(dp), intent(in) :: z
    integer                 :: nstar
    real(dp)                :: absz, incr

    absz = abs(z)
    incr = max(16.0d0, 15.0d0 * max(1.0d0, absz)**(1.0d0/3.0d0))
    nstar = max(nmax, nint(absz)) + ceiling(incr)
  end function calc_nstar

  !----------------------------------------------------------------------------
  ! Downward recurrence for D1, the logarithmic derivative of the
  ! Riccati-Bessel function of the first kind.  Eqs. (16a), (16b).
  !
  ! D1[nstar] = 0
  ! D1[n-1]   = n/z - 1 / (D1[n] + n/z)     n = nstar..1
  ! D1[0]     = cot(z)   (stable upward value, see Du 2004)
  !----------------------------------------------------------------------------
  subroutine evalDownwardD1(z, nmax, D1)
    complex(dp), intent(in)  :: z
    integer,     intent(in)  :: nmax
    complex(dp), intent(out) :: D1(0:nmax)

    complex(dp), allocatable :: tmp(:)
    complex(dp) :: z_inv
    integer     :: nstar, n

    nstar = calc_nstar(nmax, z)
    allocate(tmp(0:nstar))

    z_inv     = 1.0d0 / z
    tmp(nstar) = (0.0d0, 0.0d0)
    do n = nstar, 1, -1
      tmp(n-1) = dble(n)*z_inv - 1.0d0/(tmp(n) + dble(n)*z_inv)
    end do
    tmp(0) = complex_cot(z)    ! replace with stable value

    D1(0:nmax) = tmp(0:nmax)
    deallocate(tmp)
  end subroutine evalDownwardD1

  !----------------------------------------------------------------------------
  ! Upward recurrence for the product PsiZeta and logarithmic derivative D3.
  ! Eqs. (18a)-(18d).
  !
  ! PsiZeta[0] = 0.5 * (1 - exp(2iz))
  ! D3[0]      = i
  ! PsiZeta[n] = PsiZeta[n-1] * (n/z - D1[n-1]) * (n/z - D3[n-1])
  ! D3[n]      = D1[n] + i / PsiZeta[n]
  !----------------------------------------------------------------------------
  subroutine evalUpwardD3(z, nmax, D1, D3, PsiZeta)
    complex(dp), intent(in)  :: z
    integer,     intent(in)  :: nmax
    complex(dp), intent(in)  :: D1(0:nmax)
    complex(dp), intent(out) :: D3(0:nmax)
    complex(dp), intent(out) :: PsiZeta(0:nmax)

    complex(dp) :: z_inv, ci
    real(dp)    :: zr, zi
    integer     :: n

    ci  = (0.0d0, 1.0d0)
    zr  = real(z, dp)
    zi  = aimag(z)
    z_inv = 1.0d0 / z

    ! 0.5*(1 - exp(2iz)),  exp(2iz) = exp(-2*zi)*(cos(2*zr)+i*sin(2*zr))
    PsiZeta(0) = 0.5d0 * cmplx(1.0d0 - cos(2.0d0*zr)*exp(-2.0d0*zi), &
                                       -sin(2.0d0*zr)*exp(-2.0d0*zi), dp)
    D3(0) = ci

    do n = 1, nmax
      PsiZeta(n) = PsiZeta(n-1) * (dble(n)*z_inv - D1(n-1)) &
                                 * (dble(n)*z_inv - D3(n-1))
      D3(n)      = D1(n) + ci / PsiZeta(n)
    end do
  end subroutine evalUpwardD3

  !----------------------------------------------------------------------------
  ! Upward recurrence for Psi, the Riccati-Bessel function of the first kind.
  ! Eqs. (20a)-(20b).
  !
  ! Psi[0] = sin(z)
  ! Psi[n] = Psi[n-1] * (n/z - D1[n-1])
  !----------------------------------------------------------------------------
  subroutine evalUpwardPsi(z, nmax, D1, Psi)
    complex(dp), intent(in)  :: z
    integer,     intent(in)  :: nmax
    complex(dp), intent(in)  :: D1(0:nmax)
    complex(dp), intent(out) :: Psi(0:nmax)

    complex(dp) :: z_inv
    integer     :: n

    z_inv  = 1.0d0 / z
    Psi(0) = sin(z)
    do n = 1, nmax
      Psi(n) = Psi(n-1) * (dble(n)*z_inv - D1(n-1))
    end do
  end subroutine evalUpwardPsi

  !----------------------------------------------------------------------------
  ! Maximum number of multipole terms needed.
  ! Uses the Wiscombe criterion (eq. 17) and the LeRu near-field cutoff,
  ! then checks that each layer's optical thickness is covered.
  !----------------------------------------------------------------------------
  function calcNmax(L, x, m) result(nmax)
    integer,     intent(in) :: L
    real(dp),    intent(in) :: x(L)
    complex(dp), intent(in) :: m(L)
    integer                 :: nmax

    real(dp) :: xL, ri, riM1
    integer  :: i, nwis, nLeRu

    xL = x(L)

    ! Wiscombe criterion
    if (xL <= 8.0d0) then
      nwis = nint(xL + 4.0d0*xL**(1.0d0/3.0d0) + 1.0d0)
    else if (xL <= 4200.0d0) then
      nwis = nint(xL + 4.05d0*xL**(1.0d0/3.0d0) + 2.0d0)
    else
      nwis = nint(xL + 4.0d0*xL**(1.0d0/3.0d0) + 2.0d0)
    end if

    ! LeRu near-field cutoff
    nLeRu = nint(xL + 11.0d0*xL**(1.0d0/3.0d0) + 1.0d0) + 1
    nmax  = max(nwis, nLeRu)

    ! Ensure optical thickness of each layer is covered
    do i = 1, L
      ri = abs(x(i) * m(i))
      nmax = max(nmax, nint(ri))
      if (i > 1) then
        riM1 = abs(x(i-1) * m(i))
        nmax = max(nmax, nint(riM1))
      end if
    end do

    nmax = nmax + 15
  end function calcNmax

  !----------------------------------------------------------------------------
  ! Compute Mie scattering coefficients an(1:nmax), bn(1:nmax) for a
  ! multilayered sphere using Yang's (2003) recursive Q/Ha/Hb algorithm.
  !
  ! Inputs:
  !   L       - number of concentric layers
  !   x(L)    - cumulative size parameters (outer radius of each shell)
  !   m(L)    - complex refractive indices relative to surrounding medium
  !   nmax    - number of multipole terms to compute
  !
  ! Outputs:
  !   an(nmax), bn(nmax)  - electric and magnetic Mie coefficients, n=1..nmax
  !----------------------------------------------------------------------------
  subroutine calcScattCoeffs(L, x, m, nmax, an, bn)
    integer,     intent(in)  :: L, nmax
    real(dp),    intent(in)  :: x(L)
    complex(dp), intent(in)  :: m(L)
    complex(dp), intent(out) :: an(nmax), bn(nmax)

    complex(dp), allocatable :: D1_z1(:), D3_z1(:), PsiZeta_z1(:)
    complex(dp), allocatable :: D1_z2(:), D3_z2(:), PsiZeta_z2(:)
    complex(dp), allocatable :: Psi(:), Zeta(:), PsiZeta_out(:)
    complex(dp), allocatable :: Q(:)
    complex(dp), allocatable :: Ha(:,:), Hb(:,:)   ! (0:nmax-1, 0:L-1)

    complex(dp) :: z1, z2, z_out, ml, mlm1
    complex(dp) :: G1_ha, G2_ha, G1_hb, G2_hb, Temp
    complex(dp) :: term1, term2, term3, term4, Num_Q, Denom_Q
    real(dp)    :: z1r, z1i, z2r, z2i, ef, ratio_sq
    integer     :: il, n                ! il = layer index (avoids clash with L)

    allocate(D1_z1(0:nmax), D3_z1(0:nmax), PsiZeta_z1(0:nmax))
    allocate(D1_z2(0:nmax), D3_z2(0:nmax), PsiZeta_z2(0:nmax))
    allocate(Psi(0:nmax), Zeta(0:nmax), PsiZeta_out(0:nmax))
    allocate(Q(0:nmax))
    allocate(Ha(0:nmax-1, 0:L-1), Hb(0:nmax-1, 0:L-1))

    !==========================================================================
    ! Layer 1 (index 0): initialise Ha, Hb from D1 of innermost layer
    !==========================================================================
    z1 = m(1) * x(1)
    call evalDownwardD1(z1, nmax, D1_z1)
    do n = 0, nmax-1
      Ha(n, 0) = D1_z1(n+1)
      Hb(n, 0) = D1_z1(n+1)
    end do

    !==========================================================================
    ! Layers 2..L  (il = 1..L-1, 0-indexed)
    ! Propagate Ha, Hb outward using Yang's Q recurrence
    !==========================================================================
    do il = 1, L-1
      ml   = m(il+1)
      mlm1 = m(il)
      z1   = ml * x(il+1)    ! m_l  * x_l
      z2   = ml * x(il)      ! m_l  * x_{l-1}

      call evalDownwardD1(z1, nmax, D1_z1)
      call evalUpwardD3  (z1, nmax, D1_z1, D3_z1, PsiZeta_z1)

      call evalDownwardD1(z2, nmax, D1_z2)
      call evalUpwardD3  (z2, nmax, D1_z2, D3_z2, PsiZeta_z2)

      !------------------------------------------------------------------------
      ! Q[0] = PsiZeta(z2)[0] / PsiZeta(z1)[0], in numerically stable form.
      ! Factors out exp(-2*(Im(z1)-Im(z2))) to keep exponents bounded.
      ! cos(-2*x) = cos(2*x), sin(-2*x) = -sin(2*x)
      !------------------------------------------------------------------------
      z1r = real(z1, dp);   z1i = aimag(z1)
      z2r = real(z2, dp);   z2i = aimag(z2)
      ef  = exp(-2.0d0*(z1i - z2i))

      Num_Q   = cmplx( ef*(cos(2.0d0*z2r) - exp(-2.0d0*z2i)), &
                      -ef* sin(2.0d0*z2r),                      dp)
      Denom_Q = cmplx( cos(2.0d0*z1r) - exp(-2.0d0*z1i),      &
                      -sin(2.0d0*z1r),                          dp)
      Q(0) = Num_Q / Denom_Q

      !------------------------------------------------------------------------
      ! Q[n] recurrence, n = 1..nmax  (Yang 2003, eq. for Q ratio)
      ! Q[n] = Q[n-1] * (x_{l-1}/x_l)^2
      !         * (z1*D1_z1[n] + n) * (n - z1*D3_z1[n-1])
      !         / ((z2*D1_z2[n] + n) * (n - z2*D3_z2[n-1]))
      !------------------------------------------------------------------------
      ratio_sq = (x(il)/x(il+1))**2
      do n = 1, nmax
        term1 = z1*D1_z1(n)   + dble(n)
        term2 = dble(n) - z1*D3_z1(n-1)
        term3 = z2*D1_z2(n)   + dble(n)
        term4 = dble(n) - z2*D3_z2(n-1)
        Q(n)  = Q(n-1) * ratio_sq * (term1*term2) / (term3*term4)
      end do

      !------------------------------------------------------------------------
      ! Propagate Ha and Hb to this layer  (Yang 2003, Ha/Hb recurrences)
      ! Ha[il][n-1] = (G2*D1_z1[n] - Q[n]*G1*D3_z1[n]) / (G2 - Q[n]*G1)
      !------------------------------------------------------------------------
      do n = 1, nmax
        ! --- Ha ---
        G1_ha = ml*Ha(n-1, il-1)   - mlm1*D1_z2(n)
        G2_ha = ml*Ha(n-1, il-1)   - mlm1*D3_z2(n)
        Temp  = Q(n) * G1_ha
        Ha(n-1, il) = (G2_ha*D1_z1(n) - Temp*D3_z1(n)) / (G2_ha - Temp)

        ! --- Hb ---
        G1_hb = mlm1*Hb(n-1, il-1) - ml*D1_z2(n)
        G2_hb = mlm1*Hb(n-1, il-1) - ml*D3_z2(n)
        Temp  = Q(n) * G1_hb
        Hb(n-1, il) = (G2_hb*D1_z1(n) - Temp*D3_z1(n)) / (G2_hb - Temp)
      end do
    end do

    !==========================================================================
    ! Psi and Zeta for the outermost medium (real argument x(L))
    !==========================================================================
    z_out = cmplx(x(L), 0.0d0, dp)
    call evalDownwardD1(z_out, nmax, D1_z1)
    call evalUpwardPsi (z_out, nmax, D1_z1, Psi)
    call evalUpwardD3  (z_out, nmax, D1_z1, D3_z1, PsiZeta_out)
    do n = 0, nmax
      Zeta(n) = PsiZeta_out(n) / Psi(n)   ! Zeta[n] = PsiZeta[n] / Psi[n]
    end do

    !==========================================================================
    ! Final an, bn from outermost Ha, Hb and Psi, Zeta  (eqs. 5, 6)
    ! an = (Ha/m_L + n/x_L)*Psi[n] - Psi[n-1]
    !    / (Ha/m_L + n/x_L)*Zeta[n] - Zeta[n-1]
    !==========================================================================
    do n = 1, nmax
      term1 = Ha(n-1, L-1)/m(L) + dble(n)/x(L)
      an(n) = (term1*Psi(n) - Psi(n-1)) / (term1*Zeta(n) - Zeta(n-1))

      term1 = m(L)*Hb(n-1, L-1) + dble(n)/x(L)
      bn(n) = (term1*Psi(n) - Psi(n-1)) / (term1*Zeta(n) - Zeta(n-1))
    end do

    deallocate(D1_z1, D3_z1, PsiZeta_z1, D1_z2, D3_z2, PsiZeta_z2)
    deallocate(Psi, Zeta, PsiZeta_out, Q, Ha, Hb)
  end subroutine calcScattCoeffs

  !----------------------------------------------------------------------------
  ! Main nMie subroutine.
  !
  ! Inputs:
  !   L            - number of concentric layers
  !   x(L)         - cumulative size parameters: x(i) = 2*pi*r_i / lambda
  !                  x must be strictly increasing; x(L) is the outermost shell
  !   m(L)         - complex refractive indices of each layer, relative to
  !                  the surrounding medium  (m = n_layer / n_medium)
  !   nTheta       - number of scattering angles
  !   Theta(nTheta)- scattering angles in radians  [0, pi]
  !   nmax_in      - (optional) override number of multipole terms
  !                  omit or pass <= 0 to use automatic Wiscombe criterion
  !
  ! Outputs:
  !   Qext    - extinction efficiency
  !   Qsca    - scattering efficiency
  !   Qabs    - absorption efficiency  (= Qext - Qsca)
  !   Qbk     - backscattering efficiency
  !   Qpr     - radiation pressure efficiency
  !   g       - asymmetry factor  (= (Qext-Qpr)/Qsca)
  !   Albedo  - single-scattering albedo  (= Qsca/Qext)
  !   S1(nTheta) - complex amplitude scattering function, perpendicular pol.
  !   S2(nTheta) - complex amplitude scattering function, parallel pol.
  !
  ! From S1, S2 the full Mueller matrix follows:
  !   S11 = 0.5*(|S2|^2 + |S1|^2)
  !   S12 = 0.5*(|S2|^2 - |S1|^2)
  !   S33 = Re(S2*conj(S1))
  !   S34 = Im(S1*conj(S2))
  !----------------------------------------------------------------------------
  subroutine nMie(L, x, m, nTheta, Theta, &
                  Qext, Qsca, Qabs, Qbk, Qpr, g, Albedo, &
                  S1, S2, nmax_in)

    integer,     intent(in)           :: L, nTheta
    real(dp),    intent(in)           :: x(L), Theta(nTheta)
    complex(dp), intent(in)           :: m(L)
    real(dp),    intent(out)          :: Qext, Qsca, Qabs, Qbk, Qpr, g, Albedo
    complex(dp), intent(out)          :: S1(nTheta), S2(nTheta)
    integer,     intent(in), optional :: nmax_in

    complex(dp), allocatable :: an(:), bn(:)
    real(dp),    allocatable :: mu(:)          ! cos(Theta(k))
    real(dp),    allocatable :: pi_prev(:), pi_curr(:)  ! per-angle Pi recurrence

    integer     :: nmax, n, k
    real(dp)    :: xL, xL2, norm
    real(dp)    :: n2p1, factor
    real(dp)    :: Qext_sum, Qsca_sum, Qpr_sum
    complex(dp) :: Qbk_tmp
    real(dp)    :: an_sq, bn_sq, re_ab, sign_n
    complex(dp) :: an_n, bn_n, an_prev, bn_prev
    real(dp)    :: pi_next, tau_n_k

    !--- Number of multipole terms -------------------------------------------
    if (present(nmax_in)) then
      if (nmax_in > 0) then
        nmax = nmax_in
      else
        nmax = calcNmax(L, x, m)
      end if
    else
      nmax = calcNmax(L, x, m)
    end if

    allocate(an(nmax), bn(nmax))
    allocate(mu(nTheta), pi_prev(nTheta), pi_curr(nTheta))

    !--- Precompute cos(Theta) once -------------------------------------------
    do k = 1, nTheta
      mu(k) = cos(Theta(k))
    end do

    !--- Scattering coefficients ----------------------------------------------
    call calcScattCoeffs(L, x, m, nmax, an, bn)

    !--- Initialise accumulators ----------------------------------------------
    Qext_sum = 0.0d0
    Qsca_sum = 0.0d0
    Qpr_sum  = 0.0d0
    Qbk_tmp  = (0.0d0, 0.0d0)
    S1       = (0.0d0, 0.0d0)
    S2       = (0.0d0, 0.0d0)

    ! Pi angular-function recurrence state: Pi_{n-1}=0, Pi_n=1 at start
    pi_prev  = 0.0d0
    pi_curr  = 1.0d0
    an_prev  = (0.0d0, 0.0d0)
    bn_prev  = (0.0d0, 0.0d0)

    !==========================================================================
    ! Main Mie series loop  n = 1 .. nmax
    ! Accumulates Qext, Qsca, Qpr, Qbk, S1, S2 simultaneously.
    !==========================================================================
    do n = 1, nmax
      an_n = an(n)
      bn_n = bn(n)

      n2p1   = dble(2*n + 1)
      factor = n2p1 / dble(n*(n+1))   ! (2n+1)/(n*(n+1))

      !--- Extinction: eq. (27) -----------------------------------------------
      Qext_sum = Qext_sum + n2p1 * real(an_n + bn_n, dp)

      !--- Scattering: eq. (28) -----------------------------------------------
      an_sq    = real(an_n,dp)**2 + aimag(an_n)**2
      bn_sq    = real(bn_n,dp)**2 + aimag(bn_n)**2
      Qsca_sum = Qsca_sum + n2p1 * (an_sq + bn_sq)

      !--- Radiation pressure: eq. (29) ---------------------------------------
      if (n > 1) then
        Qpr_sum = Qpr_sum + dble((n-1)*(n+1))/dble(n) * &
                  real(an_prev*conjg(an_n) + bn_prev*conjg(bn_n), dp)
      end if
      re_ab   = real(an_n,dp)*real(bn_n,dp) + aimag(an_n)*aimag(bn_n)
      Qpr_sum = Qpr_sum + factor * re_ab

      !--- Backscattering: eq. (33), sign = (-1)^n ----------------------------
      if (mod(n,2) == 0) then
        sign_n =  1.0d0
      else
        sign_n = -1.0d0
      end if
      Qbk_tmp = Qbk_tmp + n2p1 * sign_n * (an_n - bn_n)

      !--- Angular functions Pi_n, Tau_n and accumulate S1, S2 ---------------
      ! Pi recurrence (eqs. 26a-26c), updated inline per angle.
      ! For n=1: Pi_1=1 (already in pi_curr), Tau_1=mu. No state update.
      ! For n>=2: update pi_prev, pi_curr before use.
      do k = 1, nTheta
        if (n == 1) then
          !  pi_curr(k) = 1 (initialised), tau = mu
          tau_n_k = mu(k)
          ! pi_prev, pi_curr not updated for n=1
        else
          pi_next   = (dble(2*n-1)*mu(k)*pi_curr(k) - dble(n)*pi_prev(k)) &
                      / dble(n-1)
          tau_n_k   = dble(n)*mu(k)*pi_next - dble(n+1)*pi_curr(k)
          pi_prev(k) = pi_curr(k)
          pi_curr(k) = pi_next
        end if

        S1(k) = S1(k) + factor * (an_n*pi_curr(k) + bn_n*tau_n_k)
        S2(k) = S2(k) + factor * (an_n*tau_n_k     + bn_n*pi_curr(k))
      end do

      an_prev = an_n
      bn_prev = bn_n
    end do

    !==========================================================================
    ! Finalise efficiency factors
    !==========================================================================
    xL   = x(L)
    xL2  = xL * xL
    norm = 2.0d0 / xL2

    Qext = norm * Qext_sum
    Qsca = norm * Qsca_sum
    Qabs = Qext - Qsca
    Qpr  = Qext - 4.0d0/xL2 * Qpr_sum
    Qbk  = (real(Qbk_tmp,dp)**2 + aimag(Qbk_tmp)**2) / xL2

    if (Qsca > 1.0d-12) then
      g = (Qext - Qpr) / Qsca
    else
      g = 0.0d0
    end if
    if (Qext > 1.0d-12) then
      Albedo = Qsca / Qext
    else
      Albedo = 0.0d0
    end if

    deallocate(an, bn, mu, pi_prev, pi_curr)
  end subroutine nMie

end module nmie_mod
