      SUBROUTINE ludcmp(a,n,np,indx,d)
      INTEGER n,np,indx(n),NMAX
      REAL d,a(np,np),TINY
      PARAMETER (NMAX=500,TINY=1.0e-20)
      INTEGER i,imax,j,k
      REAL aamax,dum,sum,vv(NMAX)
      d=1.
      imax=1
      do 12 i=1,n
        aamax=0.
        do 11 j=1,n
          if (abs(a(i,j)).gt.aamax) aamax=abs(a(i,j))
11      continue
c        if (aamax.eq.0.) pause 'singular matrix in ludcmp' !RT
        vv(i)=1./aamax
12    continue
      do 19 j=1,n
        do 14 i=1,j-1
          sum=a(i,j)
          do 13 k=1,i-1
            sum=sum-a(i,k)*a(k,j)
13        continue
          a(i,j)=sum
14      continue
        aamax=0.
        do 16 i=j,n
          sum=a(i,j)
          do 15 k=1,j-1
            sum=sum-a(i,k)*a(k,j)
15        continue
          a(i,j)=sum
          dum=vv(i)*abs(sum)
          if (dum.ge.aamax) then
            imax=i
            aamax=dum
          endif
16      continue
        if (j.ne.imax)then
          do 17 k=1,n
            dum=a(imax,k)
            a(imax,k)=a(j,k)
            a(j,k)=dum
17        continue
          d=-d
          vv(imax)=vv(j)
        endif
        indx(j)=imax
        if(a(j,j).eq.0.)a(j,j)=TINY
        if(j.ne.n)then
          dum=1./a(j,j)
          do 18 i=j+1,n
            a(i,j)=a(i,j)*dum
18        continue
        endif
19    continue
      return
      END
      SUBROUTINE lubksb(a,n,np,indx,b)
      INTEGER n,np,indx(n)
      REAL a(np,np),b(n)
      INTEGER i,ii,j,ll
      REAL sum
      ii=0
      do 12 i=1,n
        ll=indx(i)
        sum=b(ll)
        b(ll)=b(i)
        if (ii.ne.0)then
          do 11 j=ii,i-1
            sum=sum-a(i,j)*b(j)
11        continue
        else if (sum.ne.0.) then
          ii=i
        endif
        b(i)=sum
12    continue
      do 14 i=n,1,-1
        sum=b(i)
        do 13 j=i+1,n
          sum=sum-a(i,j)*b(j)
13      continue
        b(i)=sum/a(i,i)
14    continue
      return
      END
      SUBROUTINE gauleg(x1,x2,x,w,n)
      INTEGER n
      REAL x1,x2,x(n),w(n)
      DOUBLE PRECISION EPS
      PARAMETER (EPS=3.d-14)
      INTEGER i,j,m
      DOUBLE PRECISION p1,p2,p3,pp,xl,xm,z,z1
      m=(n+1)/2
      xm=0.5d0*(x2+x1)
      xl=0.5d0*(x2-x1)
      do 12 i=1,m
        z=cos(3.141592654d0*(i-.25d0)/(n+.5d0))
1       continue
          p1=1.d0
          p2=0.d0
          do 11 j=1,n
            p3=p2
            p2=p1
            p1=((2.d0*j-1.d0)*z*p2-(j-1.d0)*p3)/j
11        continue
          pp=n*(z*p1-p2)/(z*z-1.d0)
          z1=z
          z=z1-p1/pp
        if(abs(z-z1).gt.EPS)goto 1
c        x(i)=xm-xl*z
        x(i)=real(xm-xl*z) !RT
c        x(n+1-i)=xm+xl*z
        x(n+1-i)=real(xm+xl*z) !RT
        !w(i)=2.d0*xl/((1.d0-z*z)*pp*pp)
        w(i)=real(2.d0*xl/((1.d0-z*z)*pp*pp)) !RT
        w(n+1-i)=w(i)
12    continue
      return
      END
C        PROGRAM MLPMN
C
C       ==========================================================
C       Purpose: This program computes the associated Legendre 
C                functions Pmn(x) and their derivatives Pmn'(x) 
C                using subroutine LPMN
C       Input :  x --- Argument of Pmn(x)
C                m --- Order of Pmn(x),  m = 0,1,2,...,n
C                n --- Degree of Pmn(x), n = 0,1,2,...,N
C       Output:  PM(m,n) --- Pmn(x)
C                PD(m,n) --- Pmn'(x)
C       Example: x = 0.50
C          Pmn(x):
C          m\n        1            2            3            4
C         --------------------------------------------------------
C           0      .500000     -.125000     -.437500     -.289063
C           1     -.866025    -1.299038     -.324760     1.353165
C           2      .000000     2.250000     5.625000     4.218750
C           3      .000000      .000000    -9.742786   -34.099750
C           4      .000000      .000000      .000000    59.062500
C
C          Pmn'(x):
C          m\n        1            2            3            4
C         --------------------------------------------------------
C           0     1.000000     1.500000      .375000    -1.562500
C           1      .577350    -1.732051    -6.278684    -5.773503
C           2      .000000    -3.000000     3.750000    33.750000
C           3      .000000      .000000    19.485572      .000000
C           4      .000000      .000000      .000000  -157.500000
C       ==========================================================
C
C        IMPLICIT DOUBLE PRECISION (P,X)
C        DIMENSION PM(0:100,0:100),PD(0:100,0:100)
C        WRITE(*,*)'  Please enter m, n and x'
C        READ(*,*) M,N,X
C        WRITE(*,*)
C        WRITE(*,*)'  m     n      x          Pmn(x)         Pmn''(x)'
C        WRITE(*,*)' ---------------------------------------------------'
C        CALL LPMN(100,M,N,X,PM,PD)
C        DO 15 J=0,N
C           WRITE(*,10)M,J,X,PM(M,J),PD(M,J)
C15      CONTINUE
C10      FORMAT(1X,I3,3X,I3,3X,F5.1,2E17.8)
C        END
C

        SUBROUTINE LPMN(MM,M,N,X,PM,PD)
C
C       =====================================================
C       Purpose: Compute the associated Legendre functions 
C                Pmn(x) and their derivatives Pmn'(x)
C       Input :  x  --- Argument of Pmn(x)
C                m  --- Order of Pmn(x),  m = 0,1,2,...,n
C                n  --- Degree of Pmn(x), n = 0,1,2,...,N
C                mm --- Physical dimension of PM and PD
C       Output:  PM(m,n) --- Pmn(x)
C                PD(m,n) --- Pmn'(x)
C       =====================================================
C
        IMPLICIT DOUBLE PRECISION (P,X)
        DIMENSION PM(0:MM,0:N),PD(0:MM,0:N)
        DO 10 I=0,N
        DO 10 J=0,M
           PM(J,I)=0.0D0
10         PD(J,I)=0.0D0
        PM(0,0)=1.0D0
        IF (DABS(X).EQ.1.0D0) THEN
           DO 15 I=1,N
              PM(0,I)=X**I
15            PD(0,I)=0.5D0*I*(I+1.0D0)*X**(I+1)
           DO 20 J=1,N
           DO 20 I=1,M
              IF (I.EQ.1) THEN
                 PD(I,J)=1.0D+300
              ELSE IF (I.EQ.2) THEN
                 PD(I,J)=-0.25D0*(J+2)*(J+1)*J*(J-1)*X**(J+1)
              ENDIF
20         CONTINUE
           RETURN
        ENDIF
        LS=1
        IF (DABS(X).GT.1.0D0) LS=-1
        XQ=DSQRT(LS*(1.0D0-X*X))
        XS=LS*(1.0D0-X*X)
        DO 30 I=1,M
30         PM(I,I)=-LS*(2.0D0*I-1.0D0)*XQ*PM(I-1,I-1)
        DO 35 I=0,M
35         PM(I,I+1)=(2.0D0*I+1.0D0)*X*PM(I,I)
        DO 40 I=0,M
        DO 40 J=I+2,N
           PM(I,J)=((2.0D0*J-1.0D0)*X*PM(I,J-1)-
     &             (I+J-1.0D0)*PM(I,J-2))/(J-I)
40      CONTINUE
        PD(0,0)=0.0D0
        DO 45 J=1,N
45         PD(0,J)=LS*J*(PM(0,J-1)-X*PM(0,J))/XS
        DO 50 I=1,M
        DO 50 J=I,N
           PD(I,J)=LS*I*X*PM(I,J)/XS+(J+I)
     &             *(J-I+1.0D0)/XQ*PM(I-1,J)
50      CONTINUE
        RETURN
        END
C       PROGRAM MLPN
C
C       ========================================================
C       Purpose: This program computes the Legendre polynomials 
C                Pn(x) and their derivatives Pn'(x) using
C                subroutine LPN
C       Input :  x --- Argument of Pn(x)
C                n --- Degree of Pn(x) ( n = 0,1,...)
C       Output:  PN(n) --- Pn(x)
C                PD(n) --- Pn'(x)
C       Example:    x = 0.5
C                  n          Pn(x)            Pn'(x)
C                ---------------------------------------
C                  0       1.00000000        .00000000
C                  1        .50000000       1.00000000
C                  2       -.12500000       1.50000000
C                  3       -.43750000        .37500000
C                  4       -.28906250      -1.56250000
C                  5        .08984375      -2.22656250
C       ========================================================
C
C        DOUBLE PRECISION PN,PD,X
C        DIMENSION PN(0:100),PD(0:100)
C        WRITE(*,*)'  Please enter Nmax and x '
C        READ(*,*)N,X
C        WRITE(*,30)X
C        WRITE(*,*)
C        CALL LPN(N,X,PN,PD)
C        WRITE(*,*)'  n         Pn(x)           Pn''(x)'
C        WRITE(*,*)'---------------------------------------'
C        DO 10 K=0,N
C10         WRITE(*,20)K,PN(K),PD(K)
C20      FORMAT(1X,I3,2E17.8)
C30      FORMAT(3X,'x =',F5.1)
C        END
C

        SUBROUTINE LPN(N,X,PN,PD)
C
C       ===============================================
C       Purpose: Compute Legendre polynomials Pn(x)
C                and their derivatives Pn'(x)
C       Input :  x --- Argument of Pn(x)
C                n --- Degree of Pn(x) ( n = 0,1,...)
C       Output:  PN(n) --- Pn(x)
C                PD(n) --- Pn'(x)
C       ===============================================
C
        IMPLICIT DOUBLE PRECISION (P,X)
        DIMENSION PN(0:N),PD(0:N)
        PN(0)=1.0D0
        PN(1)=X
        PD(0)=0.0D0
        PD(1)=1.0D0
        P0=1.0D0
        P1=X
        DO 10 K=2,N
           PF=(2.0D0*K-1.0D0)/K*X*P1-(K-1.0D0)/K*P0
           PN(K)=PF
           IF (DABS(X).EQ.1.0D0) THEN
              PD(K)=0.5D0*X**(K+1)*K*(K+1.0D0)
           ELSE
              PD(K)=K*(P1-X*PF)/(1.0D0-X*X)
           ENDIF
           P0=P1
10         P1=PF
        RETURN
        END
