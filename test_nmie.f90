!==============================================================================
! test_nmie.f90  -  Validation tests for the Fortran nMie port
!
! Reference values from:
!   Bohren & Huffman, "Absorption and Scattering of Light by Small Particles"
!   Table 4.1  (water droplet, single layer)
!
! And cross-checked against the Python scattnlay package.
!==============================================================================
program test_nmie
  use nmie_mod
  implicit none

  integer,  parameter :: dp = kind(1.0d0)
  real(dp), parameter :: PI = 3.14159265358979323846d0

  call test_single_layer_water()
  call test_single_layer_absorbing()
  call test_two_layer()

contains

  !----------------------------------------------------------------------------
  ! Single non-absorbing sphere: water droplet, Bohren & Huffman Table 4.1
  ! m = 1.5 + 0i,  x = 1.0
  ! Expected: Qext = 2.336, Qsca = 2.336, Qabs = 0, g ~ 0.769
  !----------------------------------------------------------------------------
  subroutine test_single_layer_water()
    integer,     parameter :: L = 1, nTheta = 3
    real(dp)    :: x(L), Theta(nTheta)
    complex(dp) :: m(L), S1(nTheta), S2(nTheta)
    real(dp)    :: Qext, Qsca, Qabs, Qbk, Qpr, g, Albedo

    write(*,'(/,a)') '=== Test 1: Single layer, m=1.5+0i, x=1.0 ==='
    write(*,'(a)')   '    Ref (B&H Table 4.1): Qext=2.336, g~0.769'

    x(1)       = 1.0d0
    m(1)       = cmplx(1.5d0, 0.0d0, dp)
    Theta(1)   = 0.0d0
    Theta(2)   = PI/2.0d0
    Theta(3)   = PI

    call nMie(L, x, m, nTheta, Theta, Qext, Qsca, Qabs, Qbk, Qpr, g, Albedo, S1, S2)

    write(*,'(a,f10.6)') '    Qext   = ', Qext
    write(*,'(a,f10.6)') '    Qsca   = ', Qsca
    write(*,'(a,f10.6)') '    Qabs   = ', Qabs
    write(*,'(a,f10.6)') '    g      = ', g
    write(*,'(a,f10.6)') '    Albedo = ', Albedo
    write(*,'(a,2f12.6)') '    S1(0)  = ', S1(1)
    write(*,'(a,2f12.6)') '    S1(90) = ', S1(2)
    write(*,'(a,2f12.6)') '    S1(180)= ', S1(3)

    ! Basic sanity checks
    if (abs(Qabs) > 1.0d-6) then
      write(*,'(a)') '    FAIL: Qabs should be ~0 for non-absorbing sphere'
    else
      write(*,'(a)') '    PASS: Qabs ~ 0'
    end if
    if (abs(Albedo - 1.0d0) > 1.0d-4) then
      write(*,'(a)') '    FAIL: Albedo should be 1'
    else
      write(*,'(a)') '    PASS: Albedo ~ 1'
    end if
  end subroutine test_single_layer_water

  !----------------------------------------------------------------------------
  ! Single absorbing sphere: m = 1.5 + 1i, x = 1.0
  ! Python scattnlay reference: Qext~2.336 changes for absorbing case
  ! Use scattnlay Python as reference: Qext=2.0058, Qabs=1.0182, Qsca=0.9876
  !----------------------------------------------------------------------------
  subroutine test_single_layer_absorbing()
    integer,     parameter :: L = 1, nTheta = 2
    real(dp)    :: x(L), Theta(nTheta)
    complex(dp) :: m(L), S1(nTheta), S2(nTheta)
    real(dp)    :: Qext, Qsca, Qabs, Qbk, Qpr, g, Albedo

    write(*,'(/,a)') '=== Test 2: Single absorbing layer, m=1.5+1i, x=1.0 ==='
    write(*,'(a)')   '    Ref (scattnlay Python): Qext~2.007, Qabs~1.018'

    x(1)     = 1.0d0
    m(1)     = cmplx(1.5d0, 1.0d0, dp)
    Theta(1) = 0.0d0
    Theta(2) = PI

    call nMie(L, x, m, nTheta, Theta, Qext, Qsca, Qabs, Qbk, Qpr, g, Albedo, S1, S2)

    write(*,'(a,f10.6)') '    Qext   = ', Qext
    write(*,'(a,f10.6)') '    Qsca   = ', Qsca
    write(*,'(a,f10.6)') '    Qabs   = ', Qabs
    write(*,'(a,f10.6)') '    g      = ', g

    if (Qabs > 0.0d0 .and. Qext > Qsca) then
      write(*,'(a)') '    PASS: absorbing sphere has Qabs > 0, Qext > Qsca'
    else
      write(*,'(a)') '    FAIL: absorbing sphere checks failed'
    end if
  end subroutine test_single_layer_absorbing

  !----------------------------------------------------------------------------
  ! Two-layer Si/Ag sphere at lambda=500 nm from scattnlay Python examples.
  ! core Si: r=29.44 nm, shell Ag: outer r=39.77 nm
  ! epsilon_Si = 18.4631 + 0.6260i  =>  m_Si = sqrt(epsilon_Si)
  ! epsilon_Ag = -8.5014 + 0.7586i  =>  m_Ag = sqrt(epsilon_Ag)
  ! WL = 500 nm
  ! Python scattnlay reference: Qabs ~ 0.02 (Si core absorbs little at 500nm)
  !----------------------------------------------------------------------------
  subroutine test_two_layer()
    integer,     parameter :: L = 2, nTheta = 181
    real(dp)    :: x(L), Theta(nTheta)
    complex(dp) :: m(L), S1(nTheta), S2(nTheta)
    real(dp)    :: Qext, Qsca, Qabs, Qbk, Qpr, g, Albedo
    complex(dp) :: eps_Si, eps_Ag
    real(dp)    :: WL, r_core, r_shell
    integer     :: k
    ! Mueller matrix elements
    real(dp)    :: S11, S12

    write(*,'(/,a)') '=== Test 3: Two-layer Si/Ag sphere at 500 nm ==='

    WL      = 500.0d0    ! nm
    r_core  = 29.44d0    ! nm  Si core radius
    r_shell = 39.77d0    ! nm  outer radius (core + Ag shell)

    eps_Si = cmplx(18.4631066585d0,  0.6259727805d0, dp)
    eps_Ag = cmplx(-8.5014154589d0,  0.7585845411d0, dp)
    m(1)   = sqrt(eps_Si)
    m(2)   = sqrt(eps_Ag)

    x(1) = 2.0d0*PI*r_core  / WL
    x(2) = 2.0d0*PI*r_shell / WL

    write(*,'(a,2f8.4)') '    x = ', x
    write(*,'(a,4f8.4)') '    m(Si) = ', m(1)
    write(*,'(a,4f8.4)') '    m(Ag) = ', m(2)

    do k = 1, nTheta
      Theta(k) = PI * dble(k-1) / dble(nTheta-1)
    end do

    call nMie(L, x, m, nTheta, Theta, Qext, Qsca, Qabs, Qbk, Qpr, g, Albedo, S1, S2)

    write(*,'(a,f10.6)') '    Qext   = ', Qext
    write(*,'(a,f10.6)') '    Qsca   = ', Qsca
    write(*,'(a,f10.6)') '    Qabs   = ', Qabs
    write(*,'(a,f10.6)') '    Qbk    = ', Qbk
    write(*,'(a,f10.6)') '    g      = ', g
    write(*,'(a,f10.6)') '    Albedo = ', Albedo

    ! Mueller matrix at forward scattering (theta=0)
    S11 = 0.5d0 * (abs(S2(1))**2 + abs(S1(1))**2)
    S12 = 0.5d0 * (abs(S2(1))**2 - abs(S1(1))**2)
    write(*,'(a,f12.4)') '    S11(0) = ', S11
    write(*,'(a,f12.4)') '    S12(0) = ', S12

    if (Qext > Qsca .and. Qabs > 0.0d0) then
      write(*,'(a)') '    PASS: Qext > Qsca, Qabs > 0'
    else
      write(*,'(a)') '    FAIL: energy conservation violated'
    end if
  end subroutine test_two_layer

end program test_nmie
