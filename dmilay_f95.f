c-----------------------------------------------------------------------
c-----------------------------------------------------------------------


      SUBROUTINE DMiLay( RCORE, RSHELL, WVNO, RINDSH, RINDCO, MU,
     &                   NUMANG, QEXT, QSCA, QBS, GQSC,
     &                   M1, M2, S21, D21, MAXANG, Err)

c **********************************************************************
c    DOUBLE PRECISION version of MieLay, which computes electromagnetic
c    scattering by a stratified sphere, i.e. a particle consisting of a
c    spherical core surrounded by a spherical shell.  The surrounding
c    medium is assumed to have refractive index unity.  The formulas,
c    manipulated to avoid the ill-conditioning that plagued earlier
c    formulations, were published in:

c        Toon, O. and T. Ackerman, Applied Optics 20, 3657 (1981)

c    Further documentation, including definitons of input and output
c    arguments, is inside the single precision version of this program
c    (SUBROUTINE MieLay, available by anonymous ftp from
c    climate.gsfc.nasa.gov in directory pub/wiscombe).

c    It is recommended to use this DOUBLE PRECISION version for IEEE
c    arithmetic (32-bit floating point) computers, just to be safe.
c    If computer time is critical, back-to-back tests with the single
c    precision version should be done for ranges of radii and refractive
c    index relevant to your particular problem, before adopting the
c    single precision version.  This version is also recommended for
c    cases of large size parameter (bigger than 10 or so) and/or large
c    imaginary refractive index (bigger than 1 or so) and also whenever
c    overflows or strange behavior are encountered in running the
c    single precision version.  Sometimes the bigger exponent range in
c    DOUBLE PRECISION is as important as the added precision.

c    This version is designed to be interchangeable with the single
c    precision version:  all the floating-point input arguments are
c    still single precision.  Only the name of the routine has been
c    changed to prevent confusion (and it is strongly urged not to
c    change it to MieLay for the same reason).

c **********************************************************************

c     .. Parameters ..

      INTEGER   MXANG, LL, Err
      DOUBLE PRECISION ZERO, ONE, TWO
      PARAMETER ( MXANG = 200, LL = 200000, ZERO = 0.0D0, ONE = 1.0D0,
     &            TWO = 2.0D0 )
c     ..
c     .. Scalar Arguments ..

      INTEGER   MAXANG, NUMANG
      INTEGER, PARAMETER        :: dp      = SELECTED_REAL_KIND(P=15)
      REAL (KIND=dp)      GQSC, QBS, QEXT, QSCA, RCORE, RSHELL, WVNO
      COMPLEX (KIND=dp)  RINDCO, RINDSH
c     ..
c     .. Array Arguments ..

      REAL (KIND=dp) MU( NUMANG ), D21( MAXANG, 2 ), M1( MAXANG, 2 ),
     &          M2( MAXANG, 2 ), S21( MAXANG, 2 )
c     ..
c     .. Local Scalars ..

      LOGICAL   INPERR, PASS1
      INTEGER   J, K, M, N, NMX1, NMX2, NN

      DOUBLE PRECISION  AA, AIM, AM1IM, AM1RE, ARE, BB, BIM, BM1IM,
     &                  BM1RE, BRE, CC, COSX1, COSX4, DD, DENOM,
     &                  DGQSC, DQEXT, DQSCA, E2Y1,
     &                  EY1, EY1MY4, EY1PY4, EY4, FOURPI, PINUM,
     &                  RMM, RX, SINX1, SINX4, TOLER, X1, X4,
     &                  XCORE, XSHELL, Y1, Y4

      DOUBLE COMPLEX  AC, ACOE, ACOEM1, BC, BCOE, BCOEM1, CI, CZERO,
     &                DH1, DH2, DH4, DUMMY, DUMSQ, K1, K2, K3,
     &                P24H21, P24H24, RRFX, SBACK, WM1
c     ..
c     .. Local Arrays ..

      DOUBLE PRECISION  PI( MXANG, 3 ), SI2THT( MXANG ), T( 5 ),
     &                  TA( 4 ), TAU( MXANG, 3 )

      DOUBLE COMPLEX  S1( MXANG, 2 ), S2( MXANG, 2 ),
     &                U( 8 ), WFN( 2 ), Z( 4 )
	double complex,allocatable :: W(:,:),acap(:)
c     ..
c     .. External Functions ..

      LOGICAL   WRTBAD, WRTDIM
      EXTERNAL  WRTBAD, WRTDIM
c     ..
c     .. External Subroutines ..

      EXTERNAL  ERRMSG
c     ..
c     .. Intrinsic Functions ..

      INTRINSIC ABS, DIMAG, ASIN, DCMPLX, COS, EXP, MOD, DBLE, SIN
c     ..
c     .. Save statement ..

      SAVE  PINUM, PASS1
c     ..
c     .. Data statements ..

      DATA      PASS1 / .True. / , TOLER / 1.D-6 / ,
     &          CZERO / ( 0.D0, 0.D0 ) / , CI / ( 0.D0, 1.D0 ) /
c     ..

	allocate(w(3,LL))
	allocate(acap(LL))

      IF( PASS1 ) THEN

         PINUM  = TWO*ASIN( ONE )
         PASS1  = .False.

      END IF

      XSHELL = RSHELL*WVNO
      XCORE  = RCORE*WVNO
      T( 1 ) = XSHELL*ABS( RINDSH )
      NMX1   = 1.1D0*T( 1 )
      NMX2   = T( 1 )

      IF( NMX1.LE.150 ) THEN

         NMX1   = 150
         NMX2   = 135

      END IF

c                        ** Check input arguments for gross errors
      INPERR = .False.

      IF( WVNO.LE.0.0 ) INPERR = WRTBAD( 'WVNO' )

      IF( RSHELL.LE.0.0 ) INPERR = WRTBAD( 'Rshell' )

      IF( RCORE.LE.0.0 .OR. RCORE.GT.RSHELL )
     &    INPERR = WRTBAD( 'Rcore' )

      IF( REAL(RINDSH).LE.0.0 .OR. AIMAG(RINDSH).GT.0.0 )
     &    INPERR = WRTBAD( 'RindSh' )

      IF( REAL(RINDCO).LE.0.0 .OR. AIMAG(RINDCO).GT.0.0 )
     &    INPERR = WRTBAD( 'RindCo' )

      IF( NUMANG.LT.0 ) INPERR = WRTBAD( 'NumAng' )

      IF( NUMANG.GT.MXANG ) INPERR = WRTDIM( 'MxAng', NUMANG )

      IF( NUMANG.GT.MAXANG ) INPERR = WRTDIM( 'MaxAng', NUMANG )

      IF( NMX1 + 1 .GT. LL ) INPERR = WRTDIM( 'LL', NMX1 + 1 )

      DO 10 J = 1, NUMANG
         IF( MU(J).LT.- TOLER .OR. MU(J).GT. 1.0+TOLER )
     &        INPERR = WRTBAD( 'MU' )
   10 CONTINUE

      IF( INPERR ) then
		Err=1
		return
		CALL ERRMSG(
     &    'MIELAY--Input argument errors.  Aborting...', .True. )
	endif

      K1     = RINDCO*WVNO
      K2     = RINDSH*WVNO
      K3     = DCMPLX( WVNO )
      Z( 1 ) = RINDSH*XSHELL
      Z( 2 ) = XSHELL
      Z( 3 ) = RINDCO*XCORE
      Z( 4 ) = RINDSH*XCORE
      X1     =  DBLE( Z(1) )
      Y1     = DIMAG( Z(1) )
      X4     =  DBLE( Z(4) )
      Y4     = DIMAG( Z(4) )
      RX     = ONE / XSHELL

c                                ** Down-recurrence for A function
      ACAP( NMX1 + 1 ) = CZERO
      DO 20 M = 1, 3
         W( M, NMX1 + 1 ) = CZERO
   20 CONTINUE

      RRFX  = ONE / ( RINDSH*XSHELL)
      DO 40 NN = NMX1, 1, - 1

         ACAP( NN ) = ( ( NN + 1)*RRFX ) -
     &                ONE / ( ( (NN + 1)*RRFX) + ACAP( NN + 1) )

         DO 30 M = 1, 3

            W( M, NN ) = ( ( NN + 1) / Z( M + 1) ) -
     &                   ONE / ( ( (NN + 1)/Z(M + 1)) + W( M, NN + 1) )

   30    CONTINUE

   40 CONTINUE


      DO 50 J = 1, NUMANG

         SI2THT( J ) = ONE - MU( J )**2
         PI( J, 1 ) = ZERO
         PI( J, 2 ) = ONE
         TAU( J, 1 ) = ZERO
         TAU( J, 2 ) = MU( J )

   50 CONTINUE

c                          ** Initialization of homogeneous sphere

      T( 1 ) = COS( XSHELL )
      T( 2 ) = SIN( XSHELL )
      WM1      = DCMPLX( T(1), - T(2) )
      WFN( 1 ) = DCMPLX( T(2), T(1) )
      TA( 1 ) = T( 2 )
      TA( 2 ) = T( 1 )
      WFN( 2 ) = RX*WFN( 1 ) - WM1
      TA( 3 ) =  DBLE( WFN(2) )
      TA( 4 ) = DIMAG( WFN(2) )

c                      ** Initialization procedure for stratified sphere
      N      = 1
      SINX1  = SIN( X1 )
      SINX4  = SIN( X4 )
      COSX1  = COS( X1 )
      COSX4  = COS( X4 )
      EY1    = EXP( Y1 )
      E2Y1   = EY1**2
      EY4    = EXP( Y4 )
      EY1MY4 = EXP( Y1 - Y4 )
      EY1PY4 = EY1*EY4
      AA     = SINX4*( EY1PY4 + EY1MY4 )
      BB     = COSX4*( EY1PY4 - EY1MY4 )
      CC     = SINX1*( E2Y1 + ONE )
      DD     = COSX1*( E2Y1 - ONE )
      DENOM  = ONE + E2Y1*( 4.0D0*SINX1**2 - TWO + E2Y1 )
      DUMMY  = DCMPLX( ( AA*CC + BB*DD) / DENOM,
     &                 ( BB*CC - AA*DD) / DENOM )
      DUMMY  = DUMMY*( ACAP(N) + N / Z(1) ) / ( W(3, N) + N / Z(4) )
      DUMSQ  = DUMMY**2

      P24H24 = 0.5D0 + DCMPLX( SINX4**2 - 0.5D0, COSX4*SINX4 )*EY4**2
      P24H21 = 0.5D0*DCMPLX( SINX1*SINX4 - COSX1*COSX4,
     &                       SINX1*COSX4 + COSX1*SINX4 )*EY1PY4
     &       + 0.5D0*DCMPLX( SINX1*SINX4 + COSX1*COSX4,
     &                     - SINX1*COSX4 + COSX1*SINX4 )*EY1MY4
      DH1    = Z( 1 ) / ( ONE + CI*Z( 1) ) - ONE / Z( 1 )
      DH2    = Z( 2 ) / ( ONE + CI*Z( 2) ) - ONE / Z( 2 )
      DH4    = Z( 4 ) / ( ONE + CI*Z( 4) ) - ONE / Z( 4 )
      P24H24 = P24H24 / ( ( DH4 + N/Z(4))*( W(3, N) + N/Z(4)) )
      P24H21 = P24H21 / ( ( DH1 + N/Z(1))*( W(3, N) + N/Z(4)) )

      U( 1 ) = K3*ACAP( N ) - K2*W( 1, N )
      U( 2 ) = K3*ACAP( N ) - K2*DH2
      U( 3 ) = K2*ACAP( N ) - K3*W( 1, N )
      U( 4 ) = K2*ACAP( N ) - K3*DH2
      U( 5 ) = K1*W( 3, N ) - K2*W( 2, N )
      U( 6 ) = K2*W( 3, N ) - K1*W( 2, N )
      U( 7 ) = - CI*( DUMMY*P24H21 - P24H24 )
      U( 8 ) = TA( 3 ) / WFN( 2 )

      ACOE  = U( 8 )*( U(1)*U(5)*U(7) + K1*U(1) - DUMSQ*K3*U(5) ) /
     &               ( U(2)*U(5)*U(7) + K1*U(2) - DUMSQ*K3*U(5) )

      BCOE  = U( 8 )*( U(3)*U(6)*U(7) + K2*U(3) - DUMSQ*K2*U(6) ) /
     &               ( U(4)*U(6)*U(7) + K2*U(4) - DUMSQ*K2*U(6) )

      ACOEM1 = ACOE
      BCOEM1 = BCOE
      ARE    =  DBLE( ACOE )
      AIM    = DIMAG( ACOE )
      BRE    =  DBLE( BCOE )
      BIM    = DIMAG( BCOE )

      DQEXT  = 3.D0*( ARE + BRE )
      DQSCA  = 3.D0*( ARE**2 + AIM**2 + BRE**2 + BIM**2 )
      DGQSC  = ZERO
      SBACK  = 3.D0*( ACOE - BCOE )
      RMM    = ONE

      AC  = 1.5D0*ACOE
      BC  = 1.5D0*BCOE
      DO 60 J = 1, NUMANG

         S1( J, 1 ) = AC*PI( J, 2 ) + BC*TAU( J, 2 )
         S1( J, 2 ) = AC*PI( J, 2 ) - BC*TAU( J, 2 )
         S2( J, 1 ) = BC*PI( J, 2 ) + AC*TAU( J, 2 )
         S2( J, 2 ) = BC*PI( J, 2 ) - AC*TAU( J, 2 )

   60 CONTINUE

c ***************** Start of Mie summing loop ******************

      N  = 2
   70 CONTINUE
c                              ** Recurrences for functions little-pi,
c                                 little-tau of Mie theory
      T( 1 ) = 2*N - 1
      T( 2 ) = N - 1
      DO 80 J = 1, NUMANG

         PI( J, 3 ) = ( T( 1)*PI( J, 2)*MU( J) - N*PI( J, 1) ) / T( 2 )

         TAU( J, 3 ) = MU( J )*( PI( J, 3) - PI( J, 1) ) -
     &                 T( 1 )*SI2THT( J )*PI( J, 2 ) + TAU( J, 1 )

   80 CONTINUE

c                                 ** Here set up homogeneous sphere
      WM1    = WFN( 1 )
      WFN( 1 ) = WFN( 2 )
      WFN( 2 ) = T( 1 )*RX*WFN( 1 ) - WM1
      TA( 1 ) =  DBLE( WFN( 1) )
      TA( 2 ) = DIMAG( WFN( 1) )
      TA( 3 ) =  DBLE( WFN( 2) )
      TA( 4 ) = DIMAG( WFN( 2) )

c                                 ** Here set up stratified sphere

      DH1    = - N / Z( 1 ) + ONE / ( N / Z( 1) - DH1 )
      DH2    = - N / Z( 2 ) + ONE / ( N / Z( 2) - DH2 )
      DH4    = - N / Z( 4 ) + ONE / ( N / Z( 4) - DH4 )
      P24H24 = P24H24 / ( ( DH4 + N/Z(4))*( W(3, N) + N/Z(4)) )
      P24H21 = P24H21 / ( ( DH1 + N/Z(1))*( W(3, N) + N/Z(4)) )
      DUMMY  = DUMMY*( ACAP(N) + N / Z(1) ) / ( W(3, N) + N / Z(4) )
      DUMSQ  = DUMMY**2

      U( 1 ) = K3*ACAP( N ) - K2*W( 1, N )
      U( 2 ) = K3*ACAP( N ) - K2*DH2
      U( 3 ) = K2*ACAP( N ) - K3*W( 1, N )
      U( 4 ) = K2*ACAP( N ) - K3*DH2
      U( 5 ) = K1*W( 3, N ) - K2*W( 2, N )
      U( 6 ) = K2*W( 3, N ) - K1*W( 2, N )
      U( 7 ) = - CI*( DUMMY*P24H21 - P24H24 )
      U( 8 ) = TA( 3 ) / WFN( 2 )

      ACOE  = U( 8 )*( U(1)*U(5)*U(7) + K1*U(1) - DUMSQ*K3*U(5) ) /
     &               ( U(2)*U(5)*U(7) + K1*U(2) - DUMSQ*K3*U(5) )

      BCOE  = U( 8 )*( U(3)*U(6)*U(7) + K2*U(3) - DUMSQ*K2*U(6) ) /
     &               ( U(4)*U(6)*U(7) + K2*U(4) - DUMSQ*K2*U(6) )
      ARE  =  DBLE( ACOE )
      AIM  = DIMAG( ACOE )
      BRE  =  DBLE( BCOE )
      BIM  = DIMAG( BCOE )

c                           ** Increment sums for efficiency factors

      AM1RE  =  DBLE( ACOEM1 )
      AM1IM  = DIMAG( ACOEM1 )
      BM1RE  =  DBLE( BCOEM1 )
      BM1IM  = DIMAG( BCOEM1 )
      T( 4 ) = (2*N - ONE) / ( N*(N - ONE) )
      T( 2 ) = (N - ONE)*(N + ONE) / N
      DGQSC  = DGQSC + T( 2 )*( AM1RE*ARE + AM1IM*AIM +
     &                          BM1RE*BRE + BM1IM*BIM ) +
     &                 T( 4 )*( AM1RE*BM1RE + AM1IM*BM1IM )

      T( 3 )  = 2*N + 1
      DQEXT   = DQEXT + T( 3 )*( ARE + BRE )
      T( 4 )  = ARE**2 + AIM**2 + BRE**2 + BIM**2
      DQSCA   = DQSCA + T( 3 )*T( 4 )
      RMM     = - RMM
      SBACK  = SBACK + T( 3 ) * RMM *( ACOE - BCOE )

      T( 2 ) = N*( N + 1 )
      T( 1 ) = T( 3 ) / T( 2 )

      AC  = T( 1 )*ACOE
      BC  = T( 1 )*BCOE
      DO 90 J = 1, NUMANG
         S1( J, 1 ) = S1( J, 1 ) + AC*PI( J, 3 ) + BC*TAU( J, 3 )
         S2( J, 1 ) = S2( J, 1 ) + BC*PI( J, 3 ) + AC*TAU( J, 3 )
   90 CONTINUE

c                               ** Scattering matrix elements for
c                                  supplements of 0-90 degree scattering
c                                  angles submitted by user
      IF( MOD(N, 2).EQ.0 ) THEN

         DO 100 J = 1, NUMANG
            S1( J, 2 ) = S1( J, 2 ) - AC*PI( J, 3 ) + BC*TAU( J, 3 )
            S2( J, 2 ) = S2( J, 2 ) - BC*PI( J, 3 ) + AC*TAU( J, 3 )
  100    CONTINUE

      ELSE

         DO 110 J = 1, NUMANG
            S1( J, 2 ) = S1( J, 2 ) + AC*PI( J, 3 ) - BC*TAU( J, 3 )
            S2( J, 2 ) = S2( J, 2 ) + BC*PI( J, 3 ) - AC*TAU( J, 3 )
  110    CONTINUE

      END IF

c                                      ** Test for convergence of sums
      IF( T(4).GE.1.0D-14 ) THEN

         N  = N + 1

         IF( N.GT.NMX2 ) then
		Err=1
		return
	   	CALL ERRMSG(
     &       'MIELAY--Dimensions for W,ACAP not enough. Suggest'//
     &       ' get detailed output, modify routine', .True. )
	   endif

         DO 120 J = 1, NUMANG

            PI( J, 1 ) = PI( J, 2 )
            PI( J, 2 ) = PI( J, 3 )
            TAU( J, 1 ) = TAU( J, 2 )
            TAU( J, 2 ) = TAU( J, 3 )

  120    CONTINUE

         ACOEM1 = ACOE
         BCOEM1 = BCOE

         GO TO 70

      END IF

c ***************** End of summing loop ******************

c                            ** Transform complex scattering amplitudes
c                               into elements of real scattering matrix

      DO 140 J = 1, NUMANG

         DO 130 K = 1, 2

            M1( J, K ) = DBLE( S1(J, K) )**2 + DIMAG( S1(J, K) )**2
            M2( J, K ) = DBLE( S2(J, K) )**2 + DIMAG( S2(J, K) )**2
            S21( J, K ) = DBLE(  S1(J, K) )*DBLE(  S2(J, K) ) +
     &                    DIMAG( S1(J, K) )*DIMAG( S2(J, K) )
            D21( J, K ) = DIMAG( S1(J, K) )*DBLE( S2(J, K) ) -
     &                    DIMAG( S2(J, K) )*DBLE( S1(J, K) )

  130    CONTINUE

  140 CONTINUE


      T( 1 ) = TWO*RX**2
      QEXT   = T( 1 )*DQEXT
      QSCA   = T( 1 )*DQSCA
      GQSC   = TWO*T( 1 )*DGQSC
      SBACK  = 0.5*SBACK
      QBS    = ( DBLE(SBACK)**2 + DIMAG(SBACK)**2 ) / (PINUM*XSHELL**2)

	return
      END




      SUBROUTINE ErrMsg( MESSAG, FATAL )

c        Print out a warning or error message;  abort if error
c        after making symbolic dump (machine-specific)

c       Provenance:  the 3 error-handlers ErrMsg, WrtBad, WrtDim are
c                    borrowed from MIEV, the Wiscombe Mie program

c     .. Scalar Arguments ..

      CHARACTER MESSAG*( * )
      LOGICAL   FATAL
c     ..
c     .. Local Scalars ..

      LOGICAL   MSGLIM
      INTEGER   MAXMSG, NUMMSG
c     ..
c     .. External Subroutines ..

c cccc EXTERNAL  SYMDUMP
c     ..
c     .. Save statement ..

      SAVE      MAXMSG, NUMMSG, MSGLIM
c     ..
c     .. Data statements ..

      DATA      NUMMSG / 0 / , MAXMSG / 100 / , MSGLIM / .FALSE. /
c     ..

      IF( FATAL ) THEN

         WRITE( *, '(//,2A,//)' ) ' ****** ERROR *****  ', MESSAG

c                                 ** Example symbolic dump call for Cray
c cccc    CALL SYMDUMP( '-B -c3' )

         write(*,*) 'I should actually stop, but... whatever'!STOP

      END IF


      NUMMSG = NUMMSG + 1

      IF( MSGLIM ) RETURN

      IF( NUMMSG.LE.MAXMSG ) THEN

         WRITE( *, '(/,2A,/)' ) ' ****** WARNING *****  ', MESSAG

      ELSE

         WRITE( *, '(//,A,//)' )
     &      ' ****** TOO MANY WARNING MESSAGES --  ' //
     &      'They will no longer be printed *******'

         MSGLIM = .True.

      END IF

      END

      LOGICAL FUNCTION WrtBad( VARNAM )

c          Write names of erroneous variables and return 'TRUE'

c      INPUT :   VarNam = Name of erroneous variable to be written
c                         ( CHARACTER, any length )

c     .. Scalar Arguments ..

      CHARACTER VARNAM*( * )
c     ..
c     .. Local Scalars ..

      INTEGER   MAXMSG, NUMMSG
c     ..
c     .. External Subroutines ..

      EXTERNAL  ERRMSG
c     ..
c     .. Save statement ..

      SAVE      NUMMSG, MAXMSG
c     ..
c     .. Data statements ..

      DATA      NUMMSG / 0 / , MAXMSG / 50 /
c     ..

      WRTBAD = .TRUE.
      NUMMSG = NUMMSG + 1
      WRITE( *, '(3A)' ) ' ****  Input variable  ', VARNAM,
     &   '  in error  ****'

      IF( NUMMSG.EQ.MAXMSG ) CALL ERRMSG(
     &    'Too many input errors.  Aborting...', .TRUE. )

      END

      LOGICAL FUNCTION WrtDim( DIMNAM, MINVAL )

c          Write name of too-small symbolic dimension and
c          the value it should be increased to;  return 'TRUE'

c      INPUT :  DimNam = Name of symbolic dimension which is too small
c                        ( CHARACTER, any length )
c               Minval = Value to which that dimension should be
c                        increased (at least)

c     .. Scalar Arguments ..

      CHARACTER DIMNAM*( * )
      INTEGER   MINVAL
c     ..

      WRITE( *, '(3A,I7)' ) ' ****  Symbolic dimension  ',
     &   DIMNAM, '  should be increased to at least ', MINVAL

      WRTDIM = .TRUE.

      END
