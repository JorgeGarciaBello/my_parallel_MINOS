      SUBROUTINE EIGRS (A,N,JOBN,D,Z,IZ,WK)
      IMPLICIT REAL*8 (A-H,O-Z)
C
C-----------------------------------------------------------------------
C
C   PURPOSE             - EIGENVALUES AND (OPTIONALLY) EIGENVECTORS OF
C                           A REAL SYMMETRIC MATRIX
C
C   USAGE               - CALL EIGRS (A,N,JOBN,D,Z,IZ,WK)
C
C   ARGUMENTS    A      - INPUT REAL SYMMETRIC MATRIX OF ORDER N,
C                           WHOSE EIGENVALUES AND EIGENVECTORS
C                           ARE TO BE COMPUTED. INPUT A IS
C                           DESTROYED IF IJOB IS EQUAL TO 0 OR 1.
C                N      - INPUT ORDER OF THE MATRIX A.
C                JOBN   - INPUT OPTION PARAMETER.  IF JOBN.GE.10
C                         A IS ASSUMED TO BE IN FULL STORAGE MODE
C                         (IN THIS CASE, A MUST BE DIMENSIONED EXACTLY
C                         N BY N IN THE CALLING PROGRAM).
C                         IF JOBN.LT.10 THEN A IS ASSUMED TO BE IN
C                         SYMMETRIC STORAGE MODE.  DEFINE
C                         IJOB=MOD(JOBN,10).  THEN WHEN
C                           IJOB = 0, COMPUTE EIGENVALUES ONLY
C                           IJOB = 1, COMPUTE EIGENVALUES AND EIGEN-
C                             VECTORS.
C                           IJOB = 2, COMPUTE EIGENVALUES, EIGENVECTORS
C                             AND PERFORMANCE INDEX.
C                           IJOB = 3, COMPUTE PERFORMANCE INDEX ONLY.
C                           IF THE PERFORMANCE INDEX IS COMPUTED, IT IS
C                           RETURNED IN WK(1). THE ROUTINES HAVE
C                           PERFORMED (WELL, SATISFACTORILY, POORLY) IF
C                           WK(1) IS (LESS THAN 1, BETWEEN 1 AND 100,
C                           GREATER THAN 100).
C                D      - OUTPUT VECTOR OF LENGTH N,
C                           CONTAINING THE EIGENVALUES OF A.
C                Z      - OUTPUT N BY N MATRIX CONTAINING
C                           THE EIGENVECTORS OF A.
C                           THE EIGENVECTOR IN COLUMN J OF Z CORRES-
C                           PONDS TO THE EIGENVALUE D(J).
C                           IF IJOB = 0, Z IS NOT USED.
C                IZ     - INPUT ROW DIMENSION OF MATRIX Z EXACTLY AS
C                           SPECIFIED IN THE DIMENSION STATEMENT IN THE
C                           CALLING PROGRAM.
C                WK     - WORK AREA, THE LENGTH OF WK DEPENDS
C                           ON THE VALUE OF IJOB, WHEN
C                           IJOB = 0, THE LENGTH OF WK IS AT LEAST N.
C                           IJOB = 1, THE LENGTH OF WK IS AT LEAST N.
C                           IJOB = 2, THE LENGTH OF WK IS AT LEAST
C                             N(N+1)/2+N.
C                           IJOB = 3, THE LENGTH OF WK IS AT LEAST 1.
C
C-----------------------------------------------------------------------
C
      DIMENSION A(1),D(1),WK(1),Z(IZ,1)
      DATA RDELP/0.222045D-15/
      DATA ZERO,ONE/0.0D0,1.0D0/,TEN/10.0D0/,THOUS/1000.0D0/
C
      IF (JOBN.LT.10) GO TO 15
C                                  CONVERT TO SYMMETRIC STORAGE MODE
      K = 1
      JI = N-1
      IJ = 1
      DO 10 J=1,N
         DO 5 I=1,J
            A(K) = A(IJ)
            IJ = IJ+1
            K = K+1
    5    CONTINUE
         IJ = IJ + JI
         JI = JI - 1
   10 CONTINUE
   15 IJOB = MOD(JOBN,10)
      IF (IJOB.GE.0.AND.IJOB.LE.3) GO TO 20
      STOP 'EIGRS: IJOB NOT IN RANGE'
   20 IF (IJOB.EQ.0) GO TO 35
   25 IF (IZ.GE.N) GO TO 30
      STOP 'EIGRS: IZ .LT. N'
   30 IF (IJOB.EQ.3) GO TO 75
   35 NA = (N*(N+1))/2
      IF (IJOB.NE.2) GO TO 45
      DO 40 I=1,NA
         WK(I) = A(I)
   40 CONTINUE
C                                  SAVE INPUT A IF IJOB = 2
   45 ND = 1
      IF (IJOB.EQ.2) ND = NA+1
C                                  REDUCE A TO SYMMETRIC TRIDIAGONAL
C                                    FORM
      CALL EHOUSS (A,N,D,WK(ND),WK(ND))
      IIZ = 1
      IF (IJOB.EQ.0) GO TO 60
      IIZ = IZ
C                                  SET Z TO THE IDENTITY MATRIX
      DO 55 I=1,N
         DO 50 J=1,N
            Z(I,J) = ZERO
   50    CONTINUE
         Z(I,I) = ONE
   55 CONTINUE
C                                  COMPUTE EIGENVALUES AND EIGENVECTORS
   60 CALL EQRT2S (D,WK(ND),N,Z,IIZ)
      IF (IJOB.EQ.0) GO TO 9000
C                                  BACK TRANSFORM EIGENVECTORS
      CALL EHOBKS (A,N,1,N,Z,IZ)
   65 IF (IJOB.LE.1) GO TO 9000
C                                  MOVE INPUT MATRIX BACK TO A
      DO 70 I=1,NA
         A(I) = WK(I)
   70 CONTINUE
      WK(1) = THOUS
C                                  COMPUTE 1 - NORM OF A
   75 ANORM = ZERO
      IBEG = 1
      DO 85 I=1,N
         ASUM = ZERO
         IL = IBEG
         KK = 1
         DO 80 L=1,N
            ASUM = ASUM+DABS(A(IL))
            IF (L.GE.I) KK = L
            IL = IL+KK
   80    CONTINUE
         ANORM = DMAX1(ANORM,ASUM)
         IBEG = IBEG+I
   85 CONTINUE
      IF (ANORM.EQ.ZERO) ANORM = ONE
C                                  COMPUTE PERFORMANCE INDEX
      PI = ZERO
      DO 100 I=1,N
         IBEG = 1
         S = ZERO
         SUMZ = ZERO
         DO 95 L=1,N
            LK = IBEG
            KK = 1
            SUMZ = SUMZ+DABS(Z(L,I))
            SUMR = -D(I)*Z(L,I)
            DO 90 K=1,N
               SUMR = SUMR+A(LK)*Z(K,I)
               IF (K.GE.L) KK = K
               LK = LK+KK
   90       CONTINUE
            S = S+DABS(SUMR)
            IBEG = IBEG+L
   95    CONTINUE
         IF (SUMZ.EQ.ZERO) GO TO 100
         PI = DMAX1(PI,S/SUMZ)
  100 CONTINUE
      AN = N
      PI = PI/(ANORM*TEN*AN*RDELP)
      WK(1) = PI
      IF (JOBN.LT.10) GO TO 9000
C                                  CONVERT BACK TO FULL STORAGE MODE
      NP1 = N+1
      IJ = (N-1)*NP1 + 2
      K = (N*(NP1))/2
      DO 110 JR=1,N
         J = NP1-JR
         DO 105 IR=1,J
            IJ = IJ-1
            A(IJ) = A(K)
            K = K-1
  105    CONTINUE
         IJ = IJ-JR
  110 CONTINUE
      JI = 0
      K = N-1
      DO 120 I=1,N
         IJ = I-N
         DO 115 J=1,I
            IJ = IJ+N
            JI = JI+1
            A(IJ) = A(JI)
  115    CONTINUE
         JI = JI + K
         K = K-1
  120 CONTINUE
 9000 CONTINUE
      RETURN
      END
      SUBROUTINE EHOUSS (A,N,D,E,E2)
      IMPLICIT REAL*8 (A-H,O-Z)
C
C     CALLED ONLY BY EIGRS
C
      DIMENSION A(1),D(N),E(N),E2(N)
      DATA ZERO/0.0D0/
C
      NP1 = N+1
      NN = (N*NP1)/2-1
      NBEG = NN+1-N
      DO 70 II = 1,N
         I = NP1-II
         L = I-1
         H = ZERO
         SCALE = ZERO
         IF (L .LT. 1) GO TO 10
         NK = NN
         DO 5 K = 1,L
            SCALE = SCALE+DABS(A(NK))
            NK = NK-1
    5    CONTINUE
         IF (SCALE .NE. ZERO) GO TO 15
   10    E(I) = ZERO
         E2(I) = ZERO
         GO TO 65
   15    NK = NN
         DO 20 K = 1,L
            A(NK) = A(NK)/SCALE
            H = H+A(NK)*A(NK)
            NK = NK-1
   20    CONTINUE
         E2(I) = SCALE*SCALE*H
         F = A(NN)
         G = -DSIGN(DSQRT(H),F)
         E(I) = SCALE*G
         H = H-F*G
         A(NN) = F-G
         IF (L .EQ. 1) GO TO 55
         F = ZERO
         JK1 = 1
         DO 40 J = 1,L
            G = ZERO
            IK = NBEG+1
            JK = JK1
            DO 25 K = 1,J
               G = G+A(JK)*A(IK)
               JK = JK+1
               IK = IK+1
   25       CONTINUE
            JP1 = J+1
            IF (L .LT. JP1) GO TO 35
            JK = JK+J-1
            DO 30 K = JP1,L
               G = G+A(JK)*A(IK)
               JK = JK+K
               IK = IK+1
   30       CONTINUE
   35       E(J) = G/H
            F = F+E(J)*A(NBEG+J)
            JK1 = JK1+J
   40    CONTINUE
         HH = F/(H+H)
         JK = 1
         DO 50 J = 1,L
            F = A(NBEG+J)
            G = E(J)-HH*F
            E(J) = G
            DO 45 K = 1,J
               A(JK) = A(JK)-F*E(K)-G*A(NBEG+K)
               JK = JK+1
   45       CONTINUE
   50    CONTINUE
   55    DO 60 K = 1,L
            A(NBEG+K) = SCALE*A(NBEG+K)
   60    CONTINUE
   65    D(I) = A(NBEG+I)
         A(NBEG+I) = H*SCALE*SCALE
         NBEG = NBEG-I+1
         NN = NN-I
   70 CONTINUE
      RETURN
      END
      SUBROUTINE EHOBKS (A,N,M1,M2,Z,IZ)
      IMPLICIT REAL*8 (A-H,O-Z)
C
C-----------------------------------------------------------------------
C   PURPOSE     DERIVES EIGENVECTORS M1 TO M2 OF
C               THE ORIGINAL MATRIX FROM EIGENVECTORS
C               M1 TO M2 OF THE SYMMETRIC TRIDIAGONAL MATRIX.
C               USED ONLY BY EIGRS.
C-----------------------------------------------------------------------
C
C
      DIMENSION A(1),Z(IZ,1)
C
      IF (N .EQ. 1) GO TO 30
      DO 25 I=2,N
         L = I-1
         IA = (I*L)/2
         H = A(IA+I)
         IF (H.EQ.0.D0) GO TO 25
         DO 20 J = M1,M2
            S = 0.0D0
            DO 10 K = 1,L
               S = S+A(IA+K)*Z(K,J)
   10       CONTINUE
            S = S/H
            DO 15 K=1,L
               Z(K,J) = Z(K,J)-S*A(IA+K)
   15       CONTINUE
   20    CONTINUE
   25 CONTINUE
   30 RETURN
      END
      SUBROUTINE EQRT2S (D,E,N,Z,IZ)
      IMPLICIT REAL*8 (A-H,O-Z)
C
C-----------------------------------------------------------------------
C
C   PURPOSE             - EIGENVALUES AND (OPTIONALLY) EIGENVECTORS OF
C                           A SYMMETRIC TRIDIAGONAL MATRIX USING THE
C                           QL METHOD.
C
C   USAGE               - CALL EQRT2S (D,E,N,Z,IZ)
C
C   ARGUMENTS    D      - ON INPUT, THE VECTOR D OF LENGTH N CONTAINS
C                           THE DIAGONAL ELEMENTS OF THE SYMMETRIC
C                           TRIDIAGONAL MATRIX T.
C                           ON OUTPUT, D CONTAINS THE EIGENVALUES OF
C                           T IN ASCENDING ORDER.
C                E      - ON INPUT, THE VECTOR E OF LENGTH N CONTAINS
C                           THE SUB-DIAGONAL ELEMENTS OF T IN POSITION
C                           2,...,N. ON OUTPUT, E IS DESTROYED.
C                N      - ORDER OF TRIDIAGONAL MATRIX T.(INPUT)
C                Z      - ON INPUT, Z CONTAINS THE IDENTITY MATRIX OF
C                           ORDER N.
C                           ON OUTPUT, Z CONTAINS THE EIGENVECTORS
C                           OF T. THE EIGENVECTOR IN COLUMN J OF Z
C                           CORRESPONDS TO THE EIGENVALUE D(J).
C                IZ     - INPUT ROW DIMENSION OF MATRIX Z EXACTLY AS
C                           SPECIFIED IN THE DIMENSION STATEMENT IN THE
C                           CALLING PROGRAM. IF IZ IS LESS THAN N, THE
C                           EIGENVECTORS ARE NOT COMPUTED. IN THIS CASE
C                           Z IS NOT USED.
C
C-----------------------------------------------------------------------
C
      DIMENSION D(1),E(1),Z(IZ,1)
      DATA RDELP/0.222045D-15/
      DATA ZERO,ONE/0.0D0,1.0D0/
C                                  MOVE THE LAST N-1 ELEMENTS
C                                  OF E INTO THE FIRST N-1 LOCATIONS
      IF (N .EQ. 1) GO TO 9005
      DO 5  I=2,N
         E(I-1) = E(I)
    5 CONTINUE
      E(N) = ZERO
      B = ZERO
      F = ZERO
      DO  60  L=1,N
         J = 0
         H = RDELP*(DABS(D(L))+DABS(E(L)))
         IF (B.LT.H) B = H
C                                  LOOK FOR SMALL SUB-DIAGONAL ELEMENT
         DO 10  M=L,N
            K=M
            IF (DABS(E(K)) .LE. B) GO TO 15
   10    CONTINUE
   15    M = K
         IF (M.EQ.L) GO TO 55
   20    IF (J .EQ. 30) STOP 'EQRT2S: NON-CONVERGENCE'
         J = J+1
         L1 = L+1
         G = D(L)
         P = (D(L1)-G)/(E(L)+E(L))
         R = DABS(P)
         IF (RDELP*DABS(P) .LT. 1.0D0) R = DSQRT(P*P+ONE)
         D(L) = E(L)/(P+DSIGN(R,P))
         H = G-D(L)
         DO 25 I = L1,N
            D(I) = D(I)-H
   25    CONTINUE
         F = F+H
C                                  QL TRANSFORMATION
         P = D(M)
         C = ONE
         S = ZERO
         MM1 = M-1
         MM1PL = MM1+L
         IF (L.GT.MM1) GO TO 50
         DO 45 II=L,MM1
            I = MM1PL-II
            G = C*E(I)
            H = C*P
            IF (DABS(P).LT.DABS(E(I))) GO TO 30
            C = E(I)/P
            R = DSQRT(C*C+ONE)
            E(I+1) = S*P*R
            S = C/R
            C = ONE/R
            GO TO 35
   30       C = P/E(I)
            R = DSQRT(C*C+ONE)
            E(I+1) = S*E(I)*R
            S = ONE/R
            C = C*S
   35       P = C*D(I)-S*G
            D(I+1) = H+S*(C*G+S*D(I))
            IF (IZ .LT. N) GO TO 45
C                                  FORM VECTOR
            DO 40 K=1,N
               H = Z(K,I+1)
               Z(K,I+1) = S*Z(K,I)+C*H
               Z(K,I) = C*Z(K,I)-S*H
   40       CONTINUE
   45    CONTINUE
   50    E(L) = S*P
         D(L) = C*P
         IF (DABS(E(L)) .GT.B) GO TO 20
   55    D(L) = D(L) + F
   60 CONTINUE
C                                  ORDER EIGENVALUES AND EIGENVECTORS
      DO  80  I=1,N
         K = I
         P = D(I)
         IP1 = I+1
         IF (IP1.GT.N) GO TO 70
         DO 65  J=IP1,N
            IF (D(J) .GE. P) GO TO 65
            K = J
            P = D(J)
   65    CONTINUE
   70    IF (K.EQ.I) GO TO 80
         D(K) = D(I)
         D(I) = P
         IF (IZ .LT. N) GO TO 80
         DO 75 J = 1,N
            P = Z(J,I)
            Z(J,I) = Z(J,K)
            Z(J,K) = P
   75    CONTINUE
   80 CONTINUE
 9005 RETURN
      END
      REAL FUNCTION GGUBFS (DSEED)
C
C-----------------------------------------------------------------------
C
C   PURPOSE             - BASIC UNIFORM (0,1) RANDOM NUMBER GENERATOR -
C                           FUNCTION FORM OF GGUBS
C
C   USAGE               - FUNCTION GGUBFS (DSEED)
C
C   ARGUMENTS    GGUBFS - RESULTANT DEVIATE.
C                DSEED  - INPUT/OUTPUT DOUBLE PRECISION VARIABLE
C                           ASSIGNED AN INTEGER VALUE IN THE
C                           EXCLUSIVE RANGE (1.D0, 2147483647.D0).
C                           DSEED IS REPLACED BY A NEW VALUE TO BE
C                           USED IN A SUBSEQUENT CALL.
C
C   PRECISION/HARDWARE  - SINGLE/ALL
C
C-----------------------------------------------------------------------
C
C                                  SPECIFICATIONS FOR ARGUMENTS
      DOUBLE PRECISION   DSEED
C                                  SPECIFICATIONS FOR LOCAL VARIABLES
      DOUBLE PRECISION   D2P31M,D2P31
C                                  D2P31M=(2**31) - 1
C                                  D2P31 =(2**31)(OR AN ADJUSTED VALUE)
      DATA               D2P31M/2147483647.D0/
      DATA               D2P31 /2147483648.D0/
C                                  FIRST EXECUTABLE STATEMENT
      DSEED = DMOD(16807.D0*DSEED,D2P31M)
      GGUBFS = DSEED / D2P31
      RETURN
      END
         SUBROUTINE MATINV (A,NDM,N)
         IMPLICIT REAL*8 (A-H,O-Z)
*
*
************************************************************************
*                                                                      *
*        SUBROUTINE FOR INVERSION OF A GENERAL-REAL MATRIX.            *
*        MODIFIED BY A.S. UMAR                                         *
*                                                                      *
*        A - ON OUTPUT CONTAINS THE INVERSE MATRIX                     *
*        NDM - THE MAXIMUM DIMENSION OF A IN THE CALLING ROUTINE       *
*        N - THE ACTUAL DIMENSION USED IN CALCULATIONS (N<=NDM)        *
*                                                                      *
************************************************************************
*
*
         DIMENSION A(NDM,NDM),PIVOT(1000),INDEX(1000)
*
*        INITIALIZE PIVOT ELEMENT ARRAY
*
         DO 20 I=1,N
         PIVOT(I)=0.0D0
   20    INDEX(I)=0
*
*        PERFORM SUCCESSIVE PIVOT OPERATIONS (GRAND LOOP)
*
         DO 550 I=1,N
*
*        SEARCH FOR PIVOT ELEMENT
*
         AMAX=0.0D0
         DO 105 J=1,N
         IF (PIVOT(J).NE.0.0) GO TO 105
         DO 100 K=1,N
         IF (PIVOT(K).NE.0.0) GO TO 100
         TEMP=DABS(A(J,K))
         IF (TEMP.LT.AMAX) GO TO 100
         IROW=J
         ICOLUM=K
         AMAX=TEMP
  100    CONTINUE
  105    CONTINUE
         INDEX(I)=4096*IROW+ICOLUM
         J=IROW
         AMAX=A(J,ICOLUM)
         PIVOT(ICOLUM)=AMAX
*
*        INTERCHANGE ROWS TO PUT PIVOT ELEMENT ON DIAGONAL
*
         IF (IROW.EQ.ICOLUM) GO TO 260
         DO 200 K=1,N
         SWAP=A(J,K)
         A(J,K)=A(ICOLUM,K)
         A(ICOLUM,K)=SWAP
  200    CONTINUE
*
*        DIVIDE PIVOT ROW BY PIVOT ELEMENT
*
  260    K=ICOLUM
         A(ICOLUM,K)=1.0D0
         DO 350 K=1,N
         A(ICOLUM,K)=A(ICOLUM,K)/AMAX
  350    CONTINUE
*
*        REDUCE NON-PIVOT ROWS
*
         DO 550 J=1,N
         IF (J.EQ.ICOLUM) GO TO 550
         T=A( J,ICOLUM)
         A( J,ICOLUM)=0.0D0
         DO 450 K=1,N
         A( J,K)=A( J,K)-A(ICOLUM,K)*T
  450    CONTINUE
  550    CONTINUE
*
*     INTERCHANGE COLUMNS AFTER ALL PIVOT OPERATIONS HAVE BEEN PERFORMED
*
         DO 710 I=1,N
         I1=N+1-I
         K=INDEX(I1)/4096
         ICOLUM=INDEX(I1)-4096*K
         IF (K.EQ.ICOLUM) GO TO 710
         DO 705 J=1,N
         SWAP=A(J,K)
         A(J,K)=A(J,ICOLUM)
         A(J,ICOLUM)=SWAP
  705    CONTINUE
  710    CONTINUE
         RETURN
         END
      SUBROUTINE TB04AD(N,X,F,D,A)
      IMPLICIT REAL*8 (A-H,O-Z)
*
***********************************************************************
*     CALCULATES THE CUBIC SPLINE DERIVATIVES.                        *
*                                                                     *
*           N    -   DIMENSION OF THE ARRAYS X,F,D                    *
*           X    -   ARRAY CONTAINING KNOT VALUES IN ASCENDING ORDER  *
*           F    -   ARRAY CONTAINING THE CORRESPONDING FUNCTION VALS.*
*           D    -   OUTPUT: CONTAINS THE SPLINE FIRST DERIVATIVES    *
*           A    -   WORK ARRAY OF LENGHT 2000                        *
*                                                                     *
***********************************************************************
*
      DIMENSION A(2000),X(N),F(N),D(N)
      DO 5 I = 2,N
      IF(X(I) - X(I - 1))1,1,5
1     WRITE(6,3)I
3     FORMAT(' RETURN FROM TB04AD BECAUSE X(',I3,') OUT OF ORDER')
      A(1) = 1D0
      RETURN
5     CONTINUE
      DO 30 I = 1,N
      J = 2
      IF(I - 1)6,10,6
6     J = N - 1
      IF(I.EQ.N)GO TO 10
      H1 = 1D0/(X(I) - X(I - 1))
      H2 = 1D0/(X(I + 1) - X(I))
      A(3*I - 2) = H1
      A(3*I - 1) = 2D0*(H1 + H2)
      A(3*I) = H2
      D(I) = 3*(F(I + 1)*H2*H2 + F(I)*(H1*H1 - H2*H2) - F(I - 1)*H1*H1)
      GO TO 30
10    H1 = 1D0/(X(J) - X(J - 1))
      H2 = 1D0/(X(J + 1) - X(J))
      A(3*I - 2) = H1*H1
      A(3*I - 1) = H1*H1 - H2*H2
      A(3*I) =  - H2*H2
      D(I)=2*(F(J)*(H2*H2*H2+H1*H1*H1)-F(J+1)*H2*H2*H2-F(J-1)*H1*H1*H1)
30    CONTINUE
      P = A(4)/A(1)
      A(5) = A(5) - P*A(2)
      A(6) = A(6) - P*A(3)
      D(2) = D(2) - P*D(1)
      DO 50 I = 3,N
      K = 3*I - 4
      P = A(K + 2)/A(K)
      A(K + 3) = A(K + 3) - P*A(K + 1)
      D(I) = D(I) - P*D(I - 1)
      IF(I.NE.N - 1)GO TO 50
      P = A(K + 5)/A(K)
      A(K + 5) = A(K + 6) - P*A(K + 1)
      A(K + 6) = A(K + 7)
      D(N) = D(N) - P*D(N - 2)
50    CONTINUE
      D(N) = D(N)/A(3*N - 1)
      DO 60 I = 3,N
      J = N + 2 - I
60    D(J) = (D(J) - A(3*J)*D(J + 1))/A(3*J - 1)
      D(1) = (D(1) - D(2)*A(2) - D(3)*A(3))/A(1)
      A(1) = 0D0
      RETURN
      END
      SUBROUTINE TG02AD(N,U,S,D,X,V)
      IMPLICIT REAL*8 (A-H,O-Z)
*
***********************************************************************
*                                                                     *
*       CUBIC SPLINE INTERPOLATION                                    *
*                                                                     *
*       N   -   DIMENSION OF THE ARRAYS U,S,D                         *
*       U   -   ARRAY CONTAINING KNOT VALUES IN ASCENDING ORDER       *
*       S   -   ARRAY CONTAINING SPLINE (FUNCTION) VALUES             *
*       D   -   ARRAY CONTAINING THE FIRST DERIVATIVES FROM TB04AD    *
*       X   -   SCALAR, THE VALUE AT WHICH INTERPOLATION WILL OCCUR   *
*       V   -   OUTPUT: ARRAY OF LENGTH 4. V(1)= INTERPOLATED FUNCTION*
*               VALUE. V(2),V(3),V(4)=FIRST,SECOND,THIRD DERIVATIVEs  *
*               AT THIS POINT                                         *
*                                                                     *
***********************************************************************
*
      DIMENSION U(N),S(N),D(N),V(4)
      DATA IFLG/0/,IEPS/ - 50/
      DATA IX/-1/
      K = 0
      IF(X.LT.U(1)) GO TO 990
      IF(X.GT.U(N)) GO TO 991
      IF(IX.LT.0.OR.IFLG.EQ.0) GO TO 12
      IF(X.GT.U(J + 1)) GO TO 1
      IF(X.GE.U(J)) GO TO 18
      GO TO 2
    1 J = J + 1
   11 IF(X.GT.U(J + 1)) GO TO 1
      GO TO 7
   12 J = ABS(X - U(1))/(U(N) - U(1))*(N - 1) + 1
      J = MIN0(J,N - 1)
      IFLG = 1
      IF(X.GE.U(J)) GO TO 11
    2 J = J - 1
      IF(X.LT.U(J)) GO TO 2
    7 K = J
      H = U(J + 1) - U(J)
      HR = 1D0/H
      HRR = (HR + HR)*HR
      S0 = S(J)
      S1 = S(J + 1)
      D0 = D(J)
      D1 = D(J + 1)
      A = S1 - S0
      B = A - H*D1
      A = A - H*D0
      C = A + B
      C3 = C*3D0
   18 THETA = (X - U(J))*HR
      PHI = 1D0 - THETA
      T = THETA*PHI
      GAMA = THETA*B - PHI*A
      V(1) = THETA*S1 + PHI*S0 + T*GAMA
      V(2) = THETA*D1 + PHI*D0 + T*C3*HR
      V(3) = (C*(PHI - THETA) - GAMA)*HRR
      V(4) =  - C3*HRR*HR
      RETURN
  990 IF(X.LE.U(1) - 2D0**IEPS*MAX(ABS(U(1)),ABS(U(N)))) GO TO 99
      J = 1
      GO TO 7
  991 IF(X.GE.U(N) + 2D0**IEPS*MAX(ABS(U(1)),ABS(U(N)))) GO TO 995
      J = N - 1
      GO TO 7
  995 K = N
   99 IFLG = 0
      DO 6 I = 1,4
    6 V(I) = 0D0
      RETURN
      END
      SUBROUTINE ZSPOW (FCN,NSIG,N,ITMAX,PAR,X,FNORM,WK)
C-----------------------------------------------------------------------
C
C   USAGE               - CALL ZSPOW (FCN,NSIG,N,ITMAX,PAR,X,FNORM,
C                           WK)
C
C   ARGUMENTS    FCN    - THE NAME OF A USER-SUPPLIED SUBROUTINE WHICH
C                           EVALUATES THE SYSTEM OF EQUATIONS TO BE
C                           SOLVED. FCN MUST BE DECLARED EXTERNAL IN
C                           THE CALLING PROGRAM AND MUST HAVE THE
C                           FOLLOWING FORM,
C                             SUBROUTINE FCN(X,F,N,PAR)
C                             IMPLICIT REAL*8 (A-H,O-Z)
C                             REAL*8 X(N),F(N),PAR(1)
C                             F(1)=
C                              .
C                             F(N)=
C                             RETURN
C                             END
C                           GIVEN X(1)...X(N), FCN MUST EVALUATE THE
C                           FUNCTIONS F(1)...F(N) WHICH ARE TO BE MADE
C                           ZERO. X SHOULD NOT BE ALTERED BY FCN. THE
C                           PARAMETERS IN VECTOR PAR (SEE ARGUMENT
C                           PAR BELOW) MAY ALSO BE USED IN THE
C                           CALCULATION OF F(1)...F(N).
C                NSIG   - THE NUMBER OF DIGITS OF ACCURACY DESIRED
C                           IN THE COMPUTED ROOT. (INPUT)
C                N      - THE NUMBER OF EQUATIONS TO BE SOLVED AND
C                           THE NUMBER OF UNKNOWNS. (INPUT)
C                ITMAX  - THE MAXIMUM ALLOWABLE NUMBER OF ITERATIONS.
C                           (INPUT) THE MAXIMUM NUMBER OF CALLS TO FCN
C                           IS ITMAX*(N+1). SUGGESTED VALUE = 200.
C                PAR    - PAR CONTAINS A PARAMETER SET WHICH IS
C                           PASSED TO THE USER-SUPPLIED FUNCTION FCN.
C                           PAR MAY BE USED TO PASS ANY AUXILIARY
C                           PARAMETERS NECESSARY FOR COMPUTATION OF
C                           THE FUNCTION FCN. (INPUT)
C                X      - A VECTOR OF LENGTH N. (INPUT/OUTPUT) ON INPUT,
C                           X IS THE INITIAL APPROXIMATION TO THE ROOT.
C                           ON OUTPUT, X IS THE BEST APPROXIMATION TO
C                           THE ROOT FOUND BY ZSPOW.
C                FNORM  - ON OUTPUT, FNORM IS EQUAL TO
C                           F(1)**2+...F(N)**2 AT THE POINT X.
C                WK     - WORK VECTOR OF LENGTH N*(3*N+15)/2
C-----------------------------------------------------------------------
C                                  SPECIFICATIONS FOR ARGUMENTS
      INTEGER            NSIG,N,ITMAX,IER
      DOUBLE PRECISION   PAR(1),X(N),FNORM,WK(1)
C                                  SPECIFICATIONS FOR LOCAL VARIABLES
      INTEGER            INDEX2,INDEX,INFO,I,J,LR,MAXFEV,ML,MODE,MU,
     1                   NFEV,NPRINT
      DOUBLE PRECISION   EPSFCN,FACTOR,ONE,XTOL,ZERO
      EXTERNAL           FCN
      DATA               FACTOR,ONE,ZERO /1.0D2,1.0D0,0.0D0/
C                                  FIRST EXECUTABLE STATEMENT
      INFO = 0
C                                  CALL ZSPWA
      MAXFEV = ITMAX*(N + 1)
      XTOL = 0.1D0**NSIG
      ML = N - 1
      MU = N - 1
      EPSFCN = ZERO
      MODE = 2
      DO 5 J = 1, N
         WK(J) = ONE
    5 CONTINUE
      NPRINT = 0
      LR = (N*(N + 1))/2
      INDEX = 7*N + LR
      CALL ZSPWA(FCN,N,X,WK(6*N+1),XTOL,MAXFEV,ML,MU,EPSFCN,WK(1),
     * MODE,FACTOR,NPRINT,INFO,NFEV,WK(INDEX+1),N,WK(7*N+1),LR,
     * WK(N+1),WK(2*N+1),WK(3*N+1),WK(4*N+1),WK(5*N+1),PAR)
      IF (INFO .EQ. 5) INFO = 4
      FNORM = 0.0D0
      DO 10 I=1,N
         INDEX2 = 6*N+I
         FNORM = FNORM+WK(INDEX2)*WK(INDEX2)
   10 CONTINUE
      IF (INFO .EQ. 2) WRITE(6,*) ' ZSPOW: NUMBER OF CALLS > ITMAX*(N+1)
     X TRY A NEW INITIAL GUESS'
      IF (INFO .EQ. 3) WRITE(6,*) ' ZSPOW: NSIG IS TOO LARGE'
      IF (INFO .EQ. 4) WRITE(6,*) ' ZSPOW: ITERATION HAS NOT MADE GOOD P
     XROGRESS (TRY A NEW GUESS)'
      RETURN
      END
      SUBROUTINE ZSPWA (FCN,N,X,FVEC,XTOL,MAXFEV,ML,MU,EPSFCN,DIAG,MODE,
     *                   FACTOR,NPRINT,INFO,NFEV,FJAC,LDFJAC,R,LR,QTF,
     *                   WA1,WA2,WA3,WA4,PAR)
C                                  SPECIFICATIONS FOR ARGUMENTS
      INTEGER            N,MAXFEV,ML,MU,MODE,NPRINT,INFO,NFEV,LDFJAC,LR
      DOUBLE PRECISION   X(N),FVEC(N),XTOL,EPSFCN,DIAG(N),FACTOR,
     *                   FJAC(LDFJAC,N),R(LR),QTF(N),WA1(N),WA2(N),
     *                   WA3(N),WA4(N),PAR(1)
C                                  SPECIFICATIONS FOR LOCAL VARIABLES
      INTEGER            IFLAG,ITER,IWA(1),I,JM1,J,L,MSUM,NCFAIL,NCSUC,
     *                   NSLOW1,NSLOW2
      DOUBLE PRECISION   ACTRED,DELTA,EPSMCH,FNORM1,FNORM,ONE,P0001,
     *                   P001,P1,P5,PNORM,PRERED,RATIO,SPMPAR,SUM,TEMP,
     *                   XNORM,ZERO
      DOUBLE PRECISION   DNRM2
      LOGICAL            JEVAL,SING
      EXTERNAL           FCN
      DATA               SPMPAR /0.222045D-15/
      DATA               ONE,P1,P5,P001,P0001,ZERO /1.0D0,1.0D-1,5.0D-1,
     *                   1.0D-3,1.0D-4,0.0D0/
C                                  EPSMCH IS THE MACHINE PRECISION.
C                                  FIRST EXECUTABLE STATEMENT
      EPSMCH = SPMPAR
      INFO = 0
      IFLAG = 0
      NFEV = 0
C                                  CHECK THE INPUT PARAMETERS FOR
C                                  ERRORS.
      IF (N.LE.0 .OR. XTOL.LT.ZERO .OR. MAXFEV.LE.0 .OR. ML.LT.0 .OR.
     *MU.LT.0 .OR. FACTOR.LE.ZERO .OR. LDFJAC.LT.N .OR.
     *LR.LT.(N*(N+1))/2) GO TO 150
      IF (MODE.NE.2) GO TO 10
      DO 5 J=1,N
         IF (DIAG(J).LE.ZERO) GO TO 150
    5 CONTINUE
   10 CONTINUE
C                                  EVALUATE THE FUNCTION AT THE STARTING
C                                  POINT AND CALCULATE ITS NORM.
      IFLAG = 1
      CALL FCN(X,FVEC,N,PAR)
      NFEV = 1
      IF (IFLAG.LT.0) GO TO 150
      FNORM = DNRM2(N,FVEC,1)
C                                  DETERMINE THE NUMBER OF CALLS TO FCN
C                                  NEEDED TO COMPUTE THE JACOBIAN
C                                  MATRIX.
C
      MSUM = MIN0(ML+MU+1,N)
C
C                                  INITIALIZE ITERATION COUNTER AND
C                                  MONITORS.
      ITER = 1
      NCSUC = 0
      NCFAIL = 0
      NSLOW1 = 0
      NSLOW2 = 0
C                                  BEGINNING OF THE OUTER LOOP.
   15 CONTINUE
      JEVAL = .TRUE.
C                                  CALCULATE THE JACOBIAN MATRIX.
      IFLAG = 2
      CALL ZSPWB(FCN,N,X,FVEC,FJAC,LDFJAC,IFLAG,ML,MU,EPSFCN,WA1,WA2,
     *PAR)
      NFEV = NFEV+MSUM
      IF (IFLAG.LT.0) GO TO 150
C                                  COMPUTE THE QR FACTORIZATION OF THE
C                                  JACOBIAN.
      CALL ZSPWG(N,N,FJAC,LDFJAC,.FALSE.,IWA,1,WA1,WA2,WA3)
C                                  ON THE FIRST ITERATION AND IF MODE IS
C                                  1, SCALE ACCORDING TO THE NORMS OF
C                                  THE COLUMNS OF THE INITIAL JACOBIAN.
      IF (ITER.NE.1) GO TO 35
      IF (MODE.EQ.2) GO TO 25
      DO 20 J=1,N
         DIAG(J) = WA2(J)
         IF (WA2(J).EQ.ZERO) DIAG(J) = ONE
   20 CONTINUE
   25 CONTINUE
C                                  ON THE FIRST ITERATION, CALCULATE THE
C                                  NORM OF THE SCALED X AND INITIALIZE
C                                  THE STEP BOUND DELTA.
      DO 30 J=1,N
         WA3(J) = DIAG(J)*X(J)
   30 CONTINUE
      XNORM = DNRM2(N,WA3,1)
      DELTA = FACTOR*XNORM
      IF (DELTA.EQ.ZERO) DELTA = FACTOR
   35 CONTINUE
C                                  FORM (Q TRANSPOSE)*FVEC AND STORE IN
C                                  QTF.
      DO 40 I=1,N
         QTF(I) = FVEC(I)
   40 CONTINUE
      DO 60 J=1,N
         IF (FJAC(J,J).EQ.ZERO) GO TO 55
         SUM = ZERO
         DO 45 I=J,N
            SUM = SUM+FJAC(I,J)*QTF(I)
   45    CONTINUE
         TEMP = -SUM/FJAC(J,J)
         DO 50 I=J,N
            QTF(I) = QTF(I)+FJAC(I,J)*TEMP
   50    CONTINUE
   55    CONTINUE
   60 CONTINUE
C                                  COPY THE TRIANGULAR FACTOR OF THE QR
C                                  FACTORIZATION INTO R.
      SING = .FALSE.
      DO 75 J=1,N
         L = J
         JM1 = J-1
         IF (JM1.LT.1) GO TO 70
         DO 65 I=1,JM1
            R(L) = FJAC(I,J)
            L = L+N-I
   65    CONTINUE
   70    CONTINUE
         R(L) = WA1(J)
         IF (WA1(J).EQ.ZERO) SING = .TRUE.
   75 CONTINUE
C                                  ACCUMULATE THE ORTHOGONAL FACTOR IN
C                                  FJAC.
      CALL ZSPWF(N,N,FJAC,LDFJAC,WA1)
C                                  RESCALE IF NECESSARY.
      IF (MODE.EQ.2) GO TO 85
      DO 80 J=1,N
         DIAG(J) = DMAX1(DIAG(J),WA2(J))
   80 CONTINUE
   85 CONTINUE
C                                  BEGINNING OF THE INNER LOOP.
   90 CONTINUE
C                                  IF REQUESTED, CALL FCN TO ENABLE
C                                  PRINTING OF ITERATES.
      IF (NPRINT.LE.0) GO TO 95
      IFLAG = 0
      IF (IFLAG.LT.0) GO TO 150
   95 CONTINUE
C                                  DETERMINE THE DIRECTION P.
      CALL ZSPWC(N,R,LR,DIAG,QTF,DELTA,WA1,WA2,WA3)
C                                  STORE THE DIRECTION P AND X + P.
C                                  CALCULATE THE NORM OF P.
      DO 100 J=1,N
         WA1(J) = -WA1(J)
         WA2(J) = X(J)+WA1(J)
         WA3(J) = DIAG(J)*WA1(J)
  100 CONTINUE
      PNORM = DNRM2(N,WA3,1)
C                                  ON THE FIRST ITERATION, ADJUST THE
C                                  INITIAL STEP BOUND.
      IF (ITER.EQ.1) DELTA = DMIN1(DELTA,PNORM)
C                                  EVALUATE THE FUNCTION AT X + P AND
C                                  CALCULATE ITS NORM.
      IFLAG = 1
      CALL FCN(WA2,WA4,N,PAR)
      NFEV = NFEV+1
      IF (IFLAG.LT.0) GO TO 150
      FNORM1 = DNRM2(N,WA4,1)
C                                  COMPUTE THE SCALED ACTUAL REDUCTION.
      ACTRED = -ONE
      IF (FNORM1.LT.FNORM) ACTRED = ONE-(FNORM1/FNORM)**2
C                                  COMPUTE THE SCALED PREDICTED
C                                  REDUCTION.
      L = 1
      DO 110 I=1,N
         SUM = ZERO
         DO 105 J=I,N
            SUM = SUM+R(L)*WA1(J)
            L = L+1
  105    CONTINUE
         WA3(I) = QTF(I)+SUM
  110 CONTINUE
      TEMP = DNRM2(N,WA3,1)
      PRERED = ONE
      IF (TEMP.LT.FNORM) PRERED = ONE-(TEMP/FNORM)**2
C                                  COMPUTE THE RATIO OF THE ACTUAL TO
C                                  THE PREDICTED REDUCTION.
      RATIO = ZERO
      IF (PRERED.GT.ZERO) RATIO = ACTRED/PRERED
C                                  UPDATE THE STEP BOUND.
      IF (RATIO.GE.P1) GO TO 115
      NCSUC = 0
      NCFAIL = NCFAIL+1
      DELTA = P5*DELTA
      GO TO 120
  115 CONTINUE
      NCFAIL = 0
      NCSUC = NCSUC+1
      IF (RATIO.GE.P5 .OR. NCSUC.GT.1) DELTA = DMAX1(DELTA,PNORM/P5)
      IF (DABS(RATIO-ONE).LE.P1) DELTA = PNORM/P5
  120 CONTINUE
C                                  TEST FOR SUCCESSFUL ITERATION.
      IF (RATIO.LT.P0001) GO TO 130
C                                  SUCCESSFUL ITERATION. UPDATE X, FVEC,
C                                  AND THEIR NORMS.
      DO 125 J=1,N
         X(J) = WA2(J)
         WA2(J) = DIAG(J)*X(J)
         FVEC(J) = WA4(J)
  125 CONTINUE
      XNORM = DNRM2(N,WA2,1)
      FNORM = FNORM1
      ITER = ITER+1
  130 CONTINUE
C                                  DETERMINE THE PROGRESS OF THE
C                                  ITERATION.
      NSLOW1 = NSLOW1+1
      IF (ACTRED.GE.P001) NSLOW1 = 0
      IF (JEVAL) NSLOW2 = NSLOW2+1
      IF (ACTRED.GE.P1) NSLOW2 = 0
C                                  TEST FOR CONVERGENCE.
      IF (DELTA.LE.XTOL*XNORM .OR. FNORM.EQ.ZERO) INFO = 1
      IF (INFO.NE.0) GO TO 150
C                                  TESTS FOR TERMINATION AND STRINGENT
C                                  TOLERANCES.
      IF (NFEV.GE.MAXFEV) INFO = 2
      IF (P1*DMAX1(P1*DELTA,PNORM).LE.EPSMCH*XNORM) INFO = 3
      IF (NSLOW2.EQ.5) INFO = 4
      IF (NSLOW1.EQ.10) INFO = 5
      IF (INFO.NE.0) GO TO 150
C                                  CRITERION FOR RECALCULATING JACOBIAN
C                                  APPROXIMATION BY FORWARD DIFFERENCES.
      IF (NCFAIL.EQ.2) GO TO 145
C                                  CALCULATE THE RANK ONE MODIFICATION
C                                  TO THE JACOBIAN AND UPDATE QTF IF
C                                  NECESSARY.
      DO 140 J=1,N
         SUM = ZERO
         DO 135 I=1,N
            SUM = SUM+FJAC(I,J)*WA4(I)
  135    CONTINUE
         WA2(J) = (SUM-WA3(J))/PNORM
         WA1(J) = DIAG(J)*((DIAG(J)*WA1(J))/PNORM)
         IF (RATIO.GE.P0001) QTF(J) = SUM
  140 CONTINUE
C                                  COMPUTE THE QR FACTORIZATION OF THE
C                                  UPDATED JACOBIAN.
      CALL ZSPWE(N,N,R,LR,WA1,WA2,WA3,SING)
      CALL ZSPWD(N,N,FJAC,LDFJAC,WA2,WA3)
      CALL ZSPWD(1,N,QTF,1,WA2,WA3)
C                                  END OF THE INNER LOOP.
      JEVAL = .FALSE.
      GO TO 90
  145 CONTINUE
C                                  END OF THE OUTER LOOP.
      GO TO 15
  150 CONTINUE
C                                  TERMINATION, EITHER NORMAL OR USER
C                                  IMPOSED.
      IF (IFLAG.LT.0) INFO = IFLAG
      IFLAG = 0
      RETURN
      END
      SUBROUTINE ZSPWB (FCN,N,X,FVEC,FJAC,LDFJAC,IFLAG,ML,MU,EPSFCN,WA1,
     *                   WA2,PAR)
C                                  SPECIFICATIONS FOR ARGUMENTS
      INTEGER            N,LDFJAC,IFLAG,ML,MU
      DOUBLE PRECISION   X(N),FVEC(N),FJAC(LDFJAC,N),EPSFCN,WA1(N),
     *                   WA2(N),PAR(1)
C                                  SPECIFICATIONS FOR LOCAL VARIABLES
      INTEGER            I,J,K,MSUM
      DOUBLE PRECISION   EPSMCH,EPS,H,SPMPAR,TEMP,ZERO
      DATA               SPMPAR /0.222045D-15/
      DATA               ZERO /0.0D0/
C                                  EPSMCH IS THE MACHINE PRECISION.
C                                  FIRST EXECUTABLE STATEMENT
      EPSMCH = SPMPAR
      EPS = DSQRT(DMAX1(EPSFCN,EPSMCH))
      MSUM = ML+MU+1
      IF (MSUM.LT.N) GO TO 20
C                                  COMPUTATION OF DENSE APPROXIMATE
C                                  JACOBIAN.
      DO 10 J=1,N
         TEMP = X(J)
         H = EPS*DABS(TEMP)
         IF (H.EQ.ZERO) H = EPS
         X(J) = TEMP+H
         CALL FCN(X,WA1,N,PAR)
         IF (IFLAG.LT.0) GO TO 15
         X(J) = TEMP
         DO 5 I=1,N
            FJAC(I,J) = (WA1(I)-FVEC(I))/H
    5    CONTINUE
   10 CONTINUE
   15 CONTINUE
      GO TO 50
   20 CONTINUE
C                                  COMPUTATION OF BANDED APPROXIMATE
C                                  JACOBIAN.
      DO 40 K=1,MSUM
         DO 25 J=K,N,MSUM
            WA2(J) = X(J)
            H = EPS*DABS(WA2(J))
            IF (H.EQ.ZERO) H = EPS
            X(J) = WA2(J)+H
   25    CONTINUE
         CALL FCN(X,WA1,N,PAR)
         IF (IFLAG.LT.0) GO TO 45
         DO 35 J=K,N,MSUM
            X(J) = WA2(J)
            H = EPS*DABS(WA2(J))
            IF (H.EQ.ZERO) H = EPS
            DO 30 I=1,N
               FJAC(I,J) = ZERO
               IF (I.GE.J-MU .AND. I.LE.J+ML) FJAC(I,J) =
     *         (WA1(I)-FVEC(I))/H
   30       CONTINUE
   35    CONTINUE
   40 CONTINUE
   45 CONTINUE
   50 CONTINUE
      RETURN
      END
      SUBROUTINE ZSPWC (N,R,LR,DIAG,QTB,DELTA,X,WA1,WA2)
C                                  SPECIFICATIONS FOR ARGUMENTS
      INTEGER            N,LR
      DOUBLE PRECISION   R(LR),DIAG(N),QTB(N),DELTA,X(N),WA1(N),WA2(N)
C                                  SPECIFICATIONS FOR LOCAL VARIABLES
      INTEGER            I,JJ,JP1,J,K,L
      DOUBLE PRECISION   ALPHA,BNORM,EPSMCH,GNORM,ONE,QNORM,SGNORM,
     *                   SPMPAR,SUM,TEMP,ZERO
      DOUBLE PRECISION   DNRM2
      DATA               SPMPAR /0.222045D-15/
      DATA               ONE,ZERO /1.0D0,0.0D0/
C                                  EPSMCH IS THE MACHINE PRECISION.
C                                  FIRST EXECUTABLE STATEMENT
      EPSMCH = SPMPAR
C                                  FIRST, CALCULATE THE GAUSS-NEWTON
C                                  DIRECTION.
      JJ = (N*(N+1))/2+1
      DO 25 K=1,N
         J = N-K+1
         JP1 = J+1
         JJ = JJ-K
         L = JJ+1
         SUM = ZERO
         IF (N.LT.JP1) GO TO 10
         DO 5 I=JP1,N
            SUM = SUM+R(L)*X(I)
            L = L+1
    5    CONTINUE
   10    CONTINUE
         TEMP = R(JJ)
         IF (TEMP.NE.ZERO) GO TO 20
         L = J
         DO 15 I=1,J
            TEMP = DMAX1(TEMP,DABS(R(L)))
            L = L+N-I
   15    CONTINUE
         TEMP = EPSMCH*TEMP
         IF (TEMP.EQ.ZERO) TEMP = EPSMCH
   20    CONTINUE
         X(J) = (QTB(J)-SUM)/TEMP
   25 CONTINUE
C                                  TEST WHETHER THE GAUSS-NEWTON
C                                  DIRECTION IS ACCEPTABLE.
      DO 30 J=1,N
         WA1(J) = ZERO
         WA2(J) = DIAG(J)*X(J)
   30 CONTINUE
      QNORM = DNRM2(N,WA2,1)
      IF (QNORM.LE.DELTA) GO TO 70
C                                  THE GAUSS-NEWTON DIRECTION IS NOT
C                                  ACCEPTABLE. NEXT, CALCULATE THE
C                                  SCALED GRADIENT DIRECTION.
      L = 1
      DO 40 J=1,N
         TEMP = QTB(J)
         DO 35 I=J,N
            WA1(I) = WA1(I)+R(L)*TEMP
            L = L+1
   35    CONTINUE
         WA1(J) = WA1(J)/DIAG(J)
   40 CONTINUE
C                                  CALCULATE THE NORM OF THE SCALED
C                                  GRADIENT AND TEST FOR THE SPECIAL
C                                  CASE IN WHICH THE SCALED GRADIENT IS
C                                  ZERO.
      GNORM = DNRM2(N,WA1,1)
      SGNORM = ZERO
      ALPHA = DELTA/QNORM
      IF (GNORM.EQ.ZERO) GO TO 60
C                                  CALCULATE THE POINT ALONG THE SCALED
C                                  GRADIENT AT WHICH THE QUADRATIC IS
C                                  MINIMIZED.
      DO 45 J=1,N
         WA1(J) = (WA1(J)/GNORM)/DIAG(J)
   45 CONTINUE
      L = 1
      DO 55 J=1,N
         SUM = ZERO
         DO 50 I=J,N
            SUM = SUM+R(L)*WA1(I)
            L = L+1
   50    CONTINUE
         WA2(J) = SUM
   55 CONTINUE
      TEMP = DNRM2(N,WA2,1)
      SGNORM = (GNORM/TEMP)/TEMP
C                                  TEST WHETHER THE SCALED GRADIENT
C                                  DIRECTION IS ACCEPTABLE.
      ALPHA = ZERO
      IF (SGNORM.GE.DELTA) GO TO 60
C                                  THE SCALED GRADIENT DIRECTION IS NOT
C                                  ACCEPTABLE. FINALLY, CALCULATE THE
C                                  POINT ALONG THE DOGLEG AT WHICH THE
C                                  QUADRATIC IS MINIMIZED.
      BNORM = DNRM2(N,QTB,1)
      TEMP = (BNORM/GNORM)*(BNORM/QNORM)*(SGNORM/DELTA)
      TEMP = TEMP-(DELTA/QNORM)*(SGNORM/DELTA)**2+DSQRT((TEMP-(DELTA
     */QNORM))**2+(ONE-(DELTA/QNORM)**2)*(ONE-(SGNORM/DELTA)**2))
      ALPHA = ((DELTA/QNORM)*(ONE-(SGNORM/DELTA)**2))/TEMP
   60 CONTINUE
C                                  FORM APPROPRIATE CONVEX COMBINATION
C                                  OF THE GAUSS-NEWTON DIRECTION AND THE
C                                  SCALED GRADIENT DIRECTION.
      TEMP = (ONE-ALPHA)*DMIN1(SGNORM,DELTA)
      DO 65 J=1,N
         X(J) = TEMP*WA1(J)+ALPHA*X(J)
   65 CONTINUE
   70 CONTINUE
      RETURN
      END
      SUBROUTINE ZSPWD (M,N,A,LDA,V,W)
C                                  SPECIFICATIONS FOR ARGUMENTS
      INTEGER            M,N,LDA
      DOUBLE PRECISION   A(LDA,N),V(N),W(N)
C                                  SPECIFICATIONS FOR LOCAL VARIABLES
      INTEGER            I,J,NM1,NMJ
      DOUBLE PRECISION   TEMP1,ONE,TEMP2,TEMP
      DATA               ONE /1.0D0/
C                                  APPLY THE FIRST SET OF GIVENS
C                                  ROTATIONS TO A.
C                                  FIRST EXECUTABLE STATEMENT
      NM1 = N-1
      IF (NM1.LT.1) GO TO 25
      DO 10 NMJ=1,NM1
         J = N-NMJ
         IF (DABS(V(J)).GT.ONE) TEMP1 = ONE/V(J)
         IF (DABS(V(J)).GT.ONE) TEMP2 = DSQRT(ONE-TEMP1**2)
         IF (DABS(V(J)).LE.ONE) TEMP2 = V(J)
         IF (DABS(V(J)).LE.ONE) TEMP1 = DSQRT(ONE-TEMP2**2)
         DO 5 I=1,M
            TEMP = TEMP1*A(I,J)-TEMP2*A(I,N)
            A(I,N) = TEMP2*A(I,J)+TEMP1*A(I,N)
            A(I,J) = TEMP
    5    CONTINUE
   10 CONTINUE
C                                  APPLY THE SECOND SET OF GIVENS
C                                  ROTATIONS TO A.
      DO 20 J=1,NM1
         IF (DABS(W(J)).GT.ONE) TEMP1 = ONE/W(J)
         IF (DABS(W(J)).GT.ONE) TEMP2 = DSQRT(ONE-TEMP1**2)
         IF (DABS(W(J)).LE.ONE) TEMP2 = W(J)
         IF (DABS(W(J)).LE.ONE) TEMP1 = DSQRT(ONE-TEMP2**2)
         DO 15 I=1,M
            TEMP = TEMP1*A(I,J)+TEMP2*A(I,N)
            A(I,N) = -TEMP2*A(I,J)+TEMP1*A(I,N)
            A(I,J) = TEMP
   15    CONTINUE
   20 CONTINUE
   25 CONTINUE
      RETURN
      END
      SUBROUTINE ZSPWE (M,N,S,LS,U,V,W,SING)
C                                  SPECIFICATIONS FOR ARGUMENTS
      INTEGER            M,N,LS
      DOUBLE PRECISION   S(LS),U(M),V(N),W(M)
      LOGICAL            SING
C                                  SPECIFICATIONS FOR LOCAL VARIABLES
      INTEGER            I,JJ,J,L,NM1,NMJ
      DOUBLE PRECISION   TEMP1,TEMP2,GIANT,ONE,P25,P5,TEMP3,SPMPAR,
     *                   TEMP4,TAU,TEMP,ZERO
      DATA               GIANT /0.7237005577D+76/
      DATA               ONE,P5,P25,ZERO /1.0D0,5.0D-1,2.5D-1,0.0D0/
C                                  INITIALIZE THE DIAGONAL ELEMENT
C                                  POINTER.
C                                  FIRST EXECUTABLE STATEMENT
      JJ = (N*(2*M-N+1))/2-(M-N)
C                                  MOVE THE NONTRIVIAL PART OF THE LAST
C                                  COLUMN OF S INTO W.
      L = JJ
      DO 5 I=N,M
         W(I) = S(L)
         L = L+1
    5 CONTINUE
C                                  ROTATE THE VECTOR V INTO A MULTIPLE
C                                  OF THE N-TH UNIT VECTOR IN SUCH A WAY
C                                  THAT A SPIKE IS INTRODUCED INTO W.
      NM1 = N-1
      IF (NM1.LT.1) GO TO 35
      DO 30 NMJ=1,NM1
         J = N-NMJ
         JJ = JJ-(M-J+1)
         W(J) = ZERO
         IF (V(J).EQ.ZERO) GO TO 25
C                                  DETERMINE A GIVENS ROTATION WHICH
C                                  ELIMINATES THE J-TH ELEMENT OF V.
         IF (DABS(V(N)).GE.DABS(V(J))) GO TO 10
         TEMP2 = V(N)/V(J)
         TEMP3 = P5/DSQRT(P25+P25*TEMP2**2)
         TEMP1 = TEMP3*TEMP2
         TAU = ONE
         IF (DABS(TEMP1)*GIANT.GT.ONE) TAU = ONE/TEMP1
         GO TO 15
   10    CONTINUE
         TEMP4 = V(J)/V(N)
         TEMP1 = P5/DSQRT(P25+P25*TEMP4**2)
         TEMP3 = TEMP1*TEMP4
         TAU = TEMP3
   15    CONTINUE
C                                  APPLY THE TRANSFORMATION TO V AND
C                                  STORE THE INFORMATION NECESSARY TO
C                                  RECOVER THE GIVENS ROTATION.
         V(N) = TEMP3*V(J)+TEMP1*V(N)
         V(J) = TAU
C                                  APPLY THE TRANSFORMATION TO S AND
C                                  EXTEND THE SPIKE IN W.
         L = JJ
         DO 20 I=J,M
            TEMP = TEMP1*S(L)-TEMP3*W(I)
            W(I) = TEMP3*S(L)+TEMP1*W(I)
            S(L) = TEMP
            L = L+1
   20    CONTINUE
   25    CONTINUE
   30 CONTINUE
   35 CONTINUE
C                                  ADD THE SPIKE FROM THE RANK 1 UPDATE
C                                  TO W.
      DO 40 I=1,M
         W(I) = W(I)+V(N)*U(I)
   40 CONTINUE
C                                  ELIMINATE THE SPIKE.
      SING = .FALSE.
      IF (NM1.LT.1) GO TO 70
      DO 65 J=1,NM1
         IF (W(J).EQ.ZERO) GO TO 60
C                                  DETERMINE A GIVENS ROTATION WHICH
C                                  ELIMINATES THE J-TH ELEMENT OF THE
C                                  SPIKE.
         IF (DABS(S(JJ)).GE.DABS(W(J))) GO TO 45
         TEMP2 = S(JJ)/W(J)
         TEMP3 = P5/DSQRT(P25+P25*TEMP2**2)
         TEMP1 = TEMP3*TEMP2
         TAU = ONE
         IF (DABS(TEMP1)*GIANT.GT.ONE) TAU = ONE/TEMP1
         GO TO 50
   45    CONTINUE
         TEMP4 = W(J)/S(JJ)
         TEMP1 = P5/DSQRT(P25+P25*TEMP4**2)
         TEMP3 = TEMP1*TEMP4
         TAU = TEMP3
   50    CONTINUE
C                                  APPLY THE TRANSFORMATION TO S AND
C                                  REDUCE THE SPIKE IN W.
         L = JJ
         DO 55 I=J,M
            TEMP = TEMP1*S(L)+TEMP3*W(I)
            W(I) = -TEMP3*S(L)+TEMP1*W(I)
            S(L) = TEMP
            L = L+1
   55    CONTINUE
C                                  STORE THE INFORMATION NECESSARY TO
C                                  RECOVER THE GIVENS ROTATION.
         W(J) = TAU
   60    CONTINUE
C                                  TEST FOR ZERO DIAGONAL ELEMENTS IN
C                                  THE OUTPUT S.
         IF (S(JJ).EQ.ZERO) SING = .TRUE.
         JJ = JJ+(M-J+1)
   65 CONTINUE
   70 CONTINUE
C                                  MOVE W BACK INTO THE LAST COLUMN OF
C                                  THE OUTPUT S.
      L = JJ
      DO 75 I=N,M
         S(L) = W(I)
         L = L+1
   75 CONTINUE
      IF (S(JJ).EQ.ZERO) SING = .TRUE.
      RETURN
      END
      SUBROUTINE ZSPWF (M,N,Q,LDQ,WA)
C                                  SPECIFICATIONS FOR ARGUMENTS
      INTEGER            M,N,LDQ
      DOUBLE PRECISION   Q(LDQ,M),WA(M)
C                                  SPECIFICATIONS FOR LOCAL VARIABLES
      INTEGER            I,JM1,J,K,L,MINMN,NP1
      DOUBLE PRECISION   ONE,SUM,TEMP,ZERO
      DATA               ONE,ZERO /1.0D0,0.0D0/
C                                  ZERO OUT UPPER TRIANGLE OF Q IN THE
C                                  FIRST MIN(M,N) COLUMNS.
C                                  FIRST EXECUTABLE STATEMENT
      MINMN = MIN0(M,N)
      IF (MINMN.LT.2) GO TO 15
      DO 10 J=2,MINMN
         JM1 = J-1
         DO 5 I=1,JM1
            Q(I,J) = ZERO
    5    CONTINUE
   10 CONTINUE
   15 CONTINUE
C                                  INITIALIZE REMAINING COLUMNS TO THOSE
C                                  OF THE IDENTITY MATRIX.
      NP1 = N+1
      IF (M.LT.NP1) GO TO 30
      DO 25 J=NP1,M
         DO 20 I=1,M
            Q(I,J) = ZERO
   20    CONTINUE
         Q(J,J) = ONE
   25 CONTINUE
   30 CONTINUE
C                                  ACCUMULATE Q FROM ITS FACTORED FORM.
      DO 60 L=1,MINMN
         K = MINMN-L+1
         DO 35 I=K,M
            WA(I) = Q(I,K)
            Q(I,K) = ZERO
   35    CONTINUE
         Q(K,K) = ONE
         IF (WA(K).EQ.ZERO) GO TO 55
         DO 50 J=K,M
            SUM = ZERO
            DO 40 I=K,M
               SUM = SUM+Q(I,J)*WA(I)
   40       CONTINUE
            TEMP = SUM/WA(K)
            DO 45 I=K,M
               Q(I,J) = Q(I,J)-TEMP*WA(I)
   45       CONTINUE
   50    CONTINUE
   55    CONTINUE
   60 CONTINUE
      RETURN
      END
      SUBROUTINE ZSPWG (M,N,A,LDA,PIVOT,IPVT,LIPVT,RDIAG,ACNORM,WA)
C                                  SPECIFICATIONS FOR ARGUMENTS
      INTEGER            M,N,LDA,LIPVT,IPVT(LIPVT)
      DOUBLE PRECISION   A(LDA,N),RDIAG(N),ACNORM(N),WA(N)
      LOGICAL            PIVOT
C                                  SPECIFICATIONS FOR LOCAL VARIABLES
      INTEGER            I,JP1,J,KMAX,K,MINMN
      DOUBLE PRECISION   AJNORM,EPSMCH,ONE,P05,SPMPAR,SUM,TEMP,ZERO
      DOUBLE PRECISION   DNRM2
      DATA               SPMPAR /0.222045D-15/
      DATA               ONE,P05,ZERO /1.0D0,5.0D-2,0.0D0/
C                                  EPSMCH IS THE MACHINE PRECISION.
C                                  FIRST EXECUTABLE STATEMENT
      EPSMCH = SPMPAR
C                                  COMPUTE THE INITIAL COLUMN NORMS AND
C                                  INITIALIZE SEVERAL ARRAYS.
      DO 5 J=1,N
         ACNORM(J) = DNRM2(M,A(1,J),1)
         RDIAG(J) = ACNORM(J)
         WA(J) = RDIAG(J)
         IF (PIVOT) IPVT(J) = J
    5 CONTINUE
C                                  REDUCE A TO R WITH HOUSEHOLDER
C                                  TRANSFORMATIONS.
      MINMN = MIN0(M,N)
      DO 55 J=1,MINMN
         IF (.NOT.PIVOT) GO TO 20
C                                  BRING THE COLUMN OF LARGEST NORM INTO
C                                  THE PIVOT POSITION.
         KMAX = J
         DO 10 K=J,N
            IF (RDIAG(K).GT.RDIAG(KMAX)) KMAX = K
   10    CONTINUE
         IF (KMAX.EQ.J) GO TO 20
         DO 15 I=1,M
            TEMP = A(I,J)
            A(I,J) = A(I,KMAX)
            A(I,KMAX) = TEMP
   15    CONTINUE
         RDIAG(KMAX) = RDIAG(J)
         WA(KMAX) = WA(J)
         K = IPVT(J)
         IPVT(J) = IPVT(KMAX)
         IPVT(KMAX) = K
   20    CONTINUE
C                                  COMPUTE THE HOUSEHOLDER
C                                  TRANSFORMATION TO REDUCE THE J-TH
C                                  COLUMN OF A TO A MULTIPLE OF THE J-TH
C                                  UNIT VECTOR.
         AJNORM = DNRM2(M-J+1,A(J,J),1)
         IF (AJNORM.EQ.ZERO) GO TO 50
         IF (A(J,J).LT.ZERO) AJNORM = -AJNORM
         DO 25 I=J,M
            A(I,J) = A(I,J)/AJNORM
   25    CONTINUE
         A(J,J) = A(J,J)+ONE
C                                  APPLY THE TRANSFORMATION TO THE
C                                  REMAINING COLUMNS AND UPDATE THE
C                                  NORMS.
         JP1 = J+1
         IF (N.LT.JP1) GO TO 50
         DO 45 K=JP1,N
            SUM = ZERO
            DO 30 I=J,M
               SUM = SUM+A(I,J)*A(I,K)
   30       CONTINUE
            TEMP = SUM/A(J,J)
            DO 35 I=J,M
               A(I,K) = A(I,K)-TEMP*A(I,J)
   35       CONTINUE
            IF (.NOT.PIVOT .OR. RDIAG(K).EQ.ZERO) GO TO 40
            TEMP = A(J,K)/RDIAG(K)
            RDIAG(K) = RDIAG(K)*DSQRT(DMAX1(ZERO,ONE-TEMP**2))
            IF (P05*(RDIAG(K)/WA(K))**2.GT.EPSMCH) GO TO 40
            RDIAG(K) = DNRM2(M-J,A(JP1,K),1)
            WA(K) = RDIAG(K)
   40       CONTINUE
   45    CONTINUE
   50    CONTINUE
         RDIAG(J) = -AJNORM
   55 CONTINUE
      RETURN
      END
      DOUBLE PRECISION FUNCTION DNRM2 (N,DX,INCX)
C
C-----------------------------------------------------------------------
C
C   PURPOSE             - COMPUTE THE EUCLIDEAN LENGTH OR L2 NORM
C                           OF A DOUBLE PRECISION VECTOR
C
C   USAGE               - FUNCTION DNRM2 (N,DX,INCX)
C
C   ARGUMENTS    DNRM2  - DOUBLE PRECISION SQUARE ROOT OF THE SUM FROM
C                           I=1 TO N OF X(I)**2. (OUTPUT)
C                           X(I) REFERS TO A SPECIFIC ELEMENT OF DX.
C                           SEE INCX ARGUMENT DESCRIPTION.
C                N      - LENGTH OF VECTOR X. (INPUT)
C                DX     - DOUBLE PRECISION VECTOR OF LENGTH N*INCX.
C                           (INPUT)
C                INCX   - DISPLACEMENT BETWEEN ELEMENTS OF DX. (INPUT)
C                           X(I) IS DEFINED TO BE DX(1+(I-1)*INCX).
C                           INCX MUST BE GREATER THAN ZERO.
C
C-----------------------------------------------------------------------
C                                  SPECIFICATIONS FOR ARGUMENTS
      INTEGER            N,INCX
      DOUBLE PRECISION   DX(1)
C                                  SPECIFICATIONS FOR LOCAL VARIABLES
      INTEGER            I,J,NEXT,NN
      DOUBLE PRECISION   CUTLO,CUTHI,SUM,XMAX,ZERO,ONE,HITEST
      DATA               ZERO, ONE /0.0D0, 1.0D0/
      DATA               CUTLO, CUTHI / 8.232D-11,  1.304D19 /
C                                  FIRST EXECUTABLE STATEMENT
      IF (N.GT.0) GO TO 5
      DNRM2 = ZERO
      GO TO 70
C
    5 ASSIGN 15 TO NEXT
      SUM = ZERO
      NN = N*INCX
C                                  BEGIN MAIN LOOP
      I = 1
   10 GO TO NEXT, (15,20,35,40)
   15 IF (DABS(DX(I)).GT.CUTLO) GO TO 55
      ASSIGN 20 TO NEXT
      XMAX = ZERO
C                                  PHASE 1. SUM IS ZERO
   20 IF (DX(I).EQ.ZERO) GO TO 65
      IF (DABS(DX(I)).GT.CUTLO) GO TO 55
C                                  PREPARE FOR PHASE 2.
      ASSIGN 35 TO NEXT
      GO TO 30
C                                  PREPARE FOR PHASE 4.
   25 I = J
      ASSIGN 40 TO NEXT
      SUM = (SUM/DX(I))/DX(I)
   30 XMAX = DABS(DX(I))
      GO TO 45
C                                  PHASE 2. SUM IS SMALL. SCALE TO
C                                    AVOID DESTRUCTIVE UNDERFLOW.
   35 IF (DABS(DX(I)).GT.CUTLO) GO TO 50
C                                  COMMON CODE FOR PHASES 2 AND 4. IN
C                                    PHASE 4 SUM IS LARGE. SCALE TO
C                                    AVOID OVERFLOW.
   40 IF (DABS(DX(I)).LE.XMAX) GO TO 45
      SUM = ONE+SUM*(XMAX/DX(I))**2
      XMAX = DABS(DX(I))
      GO TO 65
C
   45 SUM = SUM+(DX(I)/XMAX)**2
      GO TO 65
C                                  PREPARE FOR PHASE 3.
   50 SUM = (SUM*XMAX)*XMAX
C                                  FOR REAL OR D.P. SET HITEST =
C                                    CUTHI/N FOR COMPLEX SET HITEST =
C                                    CUTHI/(2*N)
   55 HITEST = CUTHI/FLOAT(N)
C                                  PHASE 3. SUM IS MID-RANGE. NO
C                                    SCALING.
      DO 60 J=I,NN,INCX
         IF (DABS(DX(J)).GE.HITEST) GO TO 25
   60 SUM = SUM+DX(J)**2
      DNRM2 = DSQRT(SUM)
      GO TO 70
C
   65 CONTINUE
      I = I+INCX
      IF (I.LE.NN) GO TO 10
C                                  END OF MAIN LOOP. COMPUTE SQUARE
C                                    ROOT AND ADJUST FOR SCALING.
      DNRM2 = XMAX*DSQRT(SUM)
   70 CONTINUE
      RETURN
      END
      SUBROUTINE XFORM(YY,TI,NT,IZERO,FR,FI,FRE,NW,SUMZ)
C
C----------------------------------------------------------------------
C     SUBROUTINE FOR FAST FOURIER TRANSFORM (FFT) (GENERALIZED VERS.)
C
C     GIVEN NT POINTS FOR FOURIER TRANSFORM THE ROUTINE CHOOSES THE
C     LAST 2**N POINTS FOR FOURIER TRANSFORM.
C
C     INPUT:
C            YY   -   REAL ARRAY OF LENGTH NT. FUNCTION VALUES FOR FFT
C            TI   -   REAL ARRAY OF LENGTH NT. CORRESPONDING KNOTS
C            NT   -   NUMBER OF POINTS
C            IZERO-   IF IZERO.EQ.0 WE SUBTRACT FROM YY THE AVERAGE OF
C                     YY. (ONLY FOR THE SELECTED POINTS)
C
C     OUTPUT:
C            FR   -   REAL ARRAY OF LENGTH NT. FIRST NW ELEMENTS
C                     CONTAIN THE REAL PART OF THE FFT
C            FI   -   REAL ARRAY OF LENGTH NT. FIRST NW ELEMENTS
C                     CONTAIN THE IMAGINARY PART OF THE FFT
C            FRE  -   REAL ARRAY OF LENGTH NT. FIRST NW ELEMENTS
C                     CONTAIN THE FREQUENCIES CALCULATED BY
C                     FRE(I)=(I-1)*(2*PI/(TI(MAX)-TI(MIN))) I=1..NW
C                     THIS MAY NEED TO BE FURTHER MULTIPLIED BY A
C                     CONSTANT DEPENDING ON THE UNITS OF TI.
C            NW   -   NUMBER OF POINTS SELECTED FROM NT POINTS
C                     SATISFYING NW=2**N
C            SUMZ -   AVERAGE VALUE OF THE SELECTED PART OF YY.
C
C----------------------------------------------------------------------
C
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION TI(NT),FRE(NT),YY(NT),FR(NT),FI(NT)
      DATA PI/3.14159265358979D0/
C
C
      DFRE=2.0D0*PI
      N=NT
      FN=N
      FLNN=DLOG(FN)
      FLN2=DLOG(2.0D0)
      FN=FLNN/FLN2
      NP=FN
      NW=2**NP
      NEND=NT
      NBEG=NEND-NW+1
C
C*****INITILIZE THE REAL AND IMAGINARY PARTS OF THE FUNCTION TO
C*****BE FAST FOURIER TRANSFORMED
C
      II=0
      SUMZ=0.0D0
      DO 40 I=NBEG,NEND
      II=II+1
      FR(II)=YY(I)
      FI(II)=0.0D0
      SUMZ=SUMZ+FR(II)
 40   CONTINUE
      SUMZ=SUMZ/NW
      IF(IZERO.NE.0) GO TO 44
      DO 43 I=1,NW
      FR(I)=FR(I)-SUMZ
43    CONTINUE
44    CONTINUE
C
C*****CALCULATE THE FOURIER DIFFERENTIAL.
C
      TFOU=TI(NEND)-TI(NBEG)
      DFRE=DFRE/TFOU
C
C*****CALCULATE THE FOURIER TRANSFORM
C
      INV=2
      CALL FFT(NW,INV,FR,FI)
      DO 80 I=1,NW
      FRE(I)=(I-1)*DFRE
 80   CONTINUE
      RETURN
      END
      SUBROUTINE FFT(IT,INV,TR,TI)
C
C----------------------------------------------------------------------
C
C      THIS ROUTINE CALCULATES THE FOURIER TRANSFORM OF EQUALLY SPACE
C      F(N)  N=0,1,...,IT-1
C      THE DATA IS TAKEN TO BE PERIODIC IE.   F(N+IT) = F(N)
C***** ARGUMENTS SET BY THE CALLING PROGRAM ******
C      IT IS THE PROBLEM SIZE AND MUST BE A POWER OF 2
C      INV = 2 FOR DIRECT TRANSFORM IE.
C      G(M) = SUM OVER N=0,1,..,IT-1 OF F(N)*EXP(2PI*SQRT(-1)*N*M/IT)
C      FOR M=0,1,...,IT-1
C      INV = 1 FOR INVERSE TRANSFORM IE.
C      F(N) = (1./IT)*(SUM OVER M=0,1,..,IT-1 OF G(M)*EXP(-2PI*SQRT(-1)*
C      FOR N =0,1,...,IT-1
C      TR(I)    I=1,2,..,IT MUST CONTAIN REAL PART OF DATA
C      TI(I)    I=1,2,..,IT MUST CONTAIN THE IMAGINARY PART OF DATA
C***** ARGUMENTS SET BY ROUTINE ******
C      IF IT IS NOT A POWER OF 2 INV IS SET TO -1 FOR ERROR RETURN
C      TR(I)    I=1,2,..,IT IS SET TO REAL PART OF TRANSFORM
C      TI(I)    I=1,2,..,IT IS SET TO THE IMAGINARY PART OF TRANSFORM
C      THE METHOD USED IN THIS ROUTINE IS DISCRIBED IN
C
C----------------------------------------------------------------------
C
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION TR(4),TI(4),UR(15),UI(15)
      DATA KJUMP/1/
      GO TO (100,200),KJUMP
  100 UM=0.5D0
      DO 50 I=1,15
      UM=0.5D0*UM
      TH=6.283185307178D0*UM
      UR(I)=DCOS(TH)
   50 UI(I)=DSIN(TH)
  200 UM=1.0D0
      GO TO(1,2),INV
    1 UM=-1.0D0
    2 I0=2
      DO 3 I=2,16
      I0=I0+I0
      IF(I0-IT)3,4,5
    3 CONTINUE
C     ERROR IN IT - SET INV=-1 AND RETURN
    5 INV=-1
      RETURN
C     IT= 2**I - INITIALISE OUTER LOOP
    4 I0=I
      II=I0
      I1=IT/2
      I3=1
C     START MIDDLE LOOP
   10 K=0
      I2=I1+I1
C     CALCULATE TWIDDLE FACTOR E(K/I2)
   11 WR=1.
      WI=0.
      KK=K
      J0=I0
   24 IF(KK)21,22,21
   21 J0=J0-1
      KK1=KK
      KK=KK/2
      IF(KK1-2*KK)23,21,23
   23 WS=WR*UR(J0)-WI*UI(J0)
      WI=WR*UI(J0)+WI*UR(J0)
      WR=WS
      GO TO 24
   22 WI=WI*UM
C     START INNER LOOP
      J=0
C     DO 2*2 TRANSFORM
   31 L=J*I2+K
      L1=L+I1
      ZR=TR(L+1)+TR(L1+1)
      ZI=TI(L+1)+TI(L1+1)
      Z=WR*(TR(L+1)-TR(L1+1))-WI*(TI(L+1)-TI(L1+1))
      TI(L1+1)=WR*(TI(L+1)-TI(L1+1))+WI*(TR(L+1)-TR(L1+1))
      TR(L+1)=ZR
      TR(L1+1)=Z
      TI(L+1)=ZI
C     INDEX J LOOP
      J=J+1
      IF(J-I3)31,12,12
C     INDEX K LOOP
   12 K=K+1
      IF(K-I1)11,6,6
C     INDEX OUTER LOOP
    6 I3=I3+I3
      I0=I0-1
      I1=I1/2
      IF(I1)51,51,10
C     UNSCRAMBLE
   51 J=1
      UM=1.
      GO TO(61,52),INV
   61 UM=1./DFLOAT(IT)
   52 K=0
      J1=J
      DO 53 I=1,II
      J2=J1/2
      K=2*(K-J2)+J1
   53 J1=J2
   54 IF(K-J)66,56,55
   56 TR(J+1)=TR(J+1)*UM
      TI(J+1)=TI(J+1)*UM
      GO TO 66
   55 ZR=TR(J+1)
      ZI=TI(J+1)
      TR(J+1)=TR(K+1)*UM
      TI(J+1)=TI(K+1)*UM
      TR(K+1)=ZR*UM
      TI(K+1)=ZI*UM
   66 J=J+1
      IF(J-IT+1)52,57,57
   57 TR(1)=TR(1)*UM
      TI(1)=TI(1)*UM
      TR(IT)=TR(IT)*UM
      TI(IT)=TI(IT)*UM
      RETURN
      END
      DOUBLE PRECISION FUNCTION DGAMM(X)
C
C     DOUBLE PRECISION GAMMA FUNCTION
C
      DOUBLE PRECISION X,G,GBIG,T,XBIG,XMINV,XSMALL,Y
C
C08   DATA XSMALL/1.0D-08/
C09   DATA XSMALL/3.0D-09/
C12   DATA XSMALL/1.0D-12/
C15   DATA XSMALL/3.0D-15/
      DATA XSMALL/1.0D-17/
C19   DATA XSMALL/1.7D-18/
C
C     XBIG = LARGEST X SUCH THAT  GAMMA(X) .LT. MAXREAL
C                            AND  1.0/GAMMA(X+1.0) .GT. MINREAL
C             (ROUNDED DOWN TO AN INTEGER)
C     GBIG = GAMMA(XBIG)
C     XMINV = MAX(1.0/MAXREAL,MINREAL)  (ROUNDED UP)
C     FOR IBM 360/370 AND SIMILAR MACHINES
      DATA XBIG,GBIG,XMINV /57.0D0,7.1D+74,1.4D-76/
C     FOR DEC-10, HONEYWELL, UNIVAC 1100 (S.P.)
CR2   DATA XBIG,GBIG,XMINV /34.0D0,8.7D+36,5.9D-39/
C     FOR ICL 1900
CR3   DATA XBIG,GBIG,XMINV /58.0D0,4.0D+76,1.8D-77/
C     FOR CDC 7600/CYBER
CR4   DATA XBIG,GBIG,XMINV /164.0D0,2.0D+291,3.2D-294/
C     FOR UNIVAC 1100 (D.P.)
CR5   DATA XBIG,GBIG,XMINV /171.0D0,7.3D+306,1.2D-308/
C     FOR DEC VAX 11/780
CR6   DATA XBIG,GBIG,XMINV /34.0,8.6D+36,5.9D-39/
C
C     ERROR 1 AND 2 TEST
      T = DABS(X)
      IF (T.GT.XBIG) GO TO 160
C     SMALL RANGE TEST
      IF (T.LE.XSMALL) GO TO 140
C     MAIN RANGE REDUCTION
      M = X
      IF (X.LT.0.0D0) GO TO 80
      T = X - DBLE(FLOAT(M))
      M = M - 1
      G = 1.0D0
      IF (M) 20, 120, 40
   20 G = G/X
      GO TO 120
   40 DO 60 I=1,M
         G = G*(X-DBLE(FLOAT(I)))
   60 CONTINUE
      GO TO 120
   80 T = X - DBLE(FLOAT(M-1))
C     ERROR 4 TEST
      IF (T.EQ.1.0D0) GO TO 220
      M = 1 - M
      G = X
      DO 100 I=1,M
         G = G*(X+DBLE(FLOAT(I)))
  100 CONTINUE
      G = 1.0D0/G
  120 T = 2.0D0*T - 1.0D0
C
C
C     EXPANSION (0026) EVALUATED AS Y(T)  --PRECISION 08E.09
C08   Y=  +8.86226925D-1+T*(  +1.61692007D-2+T*(  +1.03703361D-1+
C08  AT*(  -1.34119055D-2+T*(  +9.04037536D-3+T*(  -2.42216251D-3+
C08  BT*(  +9.15547391D-4+T*(  -2.98340924D-4+T*(  +1.01593694D-4+
C08  CT*(  -3.13088821D-5+T*(  +1.03144033D-5+T*(  -5.48272091D-6+
C08  DT*(  +1.88278283D-6))))))))))))
C
C     EXPANSION (0026) EVALUATED AS Y(T)  --PRECISION 09E.10
C09   Y=  +8.862269255D-1+T*(  +1.616919866D-2+T*(  +1.037033609D-1+
C09  AT*(  -1.341184808D-2+T*(  +9.040375355D-3+T*(  -2.422622002D-3+
C09  BT*(  +9.155473906D-4+T*(  -2.967655076D-4+T*(  +1.015936944D-4+
C09  CT*(  -3.393457634D-5+T*(  +1.031440334D-5+T*(  -3.382165478D-6+
C09  DT*(  +1.882782826D-6+T*(  -6.463247484D-7)))))))))))))
C
C     EXPANSION (0026) EVALUATED AS Y(T)  --PRECISION 12E.13
C12   Y=  +8.862269254520D-1+T*(  +1.616919872669D-2+
C12  AT*(  +1.037033635205D-1+T*(  -1.341185067782D-2+
C12  BT*(  +9.040332894085D-3+T*(  -2.422593898516D-3+
C12  CT*(  +9.158021574033D-4+T*(  -2.968993359366D-4+
C12  DT*(  +1.008657892262D-4+T*(  -3.360744031186D-5+
C12  ET*(  +1.138199762073D-5+T*(  -3.810416284805D-6+
C12  FT*(  +1.106350622249D-6+T*(  -3.608242105549D-7+
C12  GT*(  +2.218377726362D-7+T*(  -7.613347676160D-8)))))))))))))))
C
C     EXPANSION (0026) EVALUATED AS Y(T)  --PRECISION 15E.16
C15   Y=  +8.862269254527580D-1+T*(  +1.616919872444243D-2+
C15  AT*(  +1.037033634220705D-1+T*(  -1.341185057058971D-2+
C15  BT*(  +9.040334940477911D-3+T*(  -2.422595384546340D-3+
C15  CT*(  +9.157859942174304D-4+T*(  -2.968901194293069D-4+
C15  DT*(  +1.009281733953869D-4+T*(  -3.363759801664768D-5+
C15  ET*(  +1.125234962812416D-5+T*(  -3.754930502328320D-6+
C15  FT*(  +1.253148247777280D-6+T*(  -4.179652784537600D-7+
C15  GT*(  +1.387603440435200D-7+T*(  -4.620920340480000D-8+
C15  HT*(  +1.613133578240000D-8+T*(  -5.419466096640000D-9+
C15  IT*(  +1.265236705280000D-9+T*( -4.030909644800000D-10+
C15  JT*( +3.622882508800000D-10+
C15  KT*( -1.243191705600000D-10)))))))))))))))))))))
C
C     EXPANSION (0026) EVALUATED AS Y(T)  --PRECISION 17E.18
      Y=  +8.86226925452758013D-1+T*(  +1.61691987244425092D-2+
     AT*(  +1.03703363422075456D-1+T*(  -1.34118505705967765D-2+
     BT*(  +9.04033494028101968D-3+T*(  -2.42259538436268176D-3+
     CT*(  +9.15785997288933120D-4+T*(  -2.96890121633200000D-4+
     DT*(  +1.00928148823365120D-4+T*(  -3.36375833240268800D-5+
     ET*(  +1.12524642975590400D-5+T*(  -3.75499034136576000D-6+
     FT*(  +1.25281466396672000D-6+T*(  -4.17808776355840000D-7+
     GT*(  +1.39383522590720000D-7+T*(  -4.64774927155200000D-8+
     HT*(  +1.53835215257600000D-8+T*(  -5.11961333760000000D-9+
     IT*(  +1.82243164160000000D-9+T*( -6.13513953280000000D-10+
     JT*( +1.27679856640000000D-10+T*( -4.01499750400000000D-11+
     KT*( +4.26560716800000000D-11+
     LT*( -1.46381209600000000D-11)))))))))))))))))))))))
C
C     EXPANSION (0026) EVALUATED AS Y(T)  --PRECISION 19E.20
C19   Y=  +8.8622692545275801366D-1+T*(  +1.6169198724442506740D-2+
C19  AT*(  +1.0370336342207529018D-1+T*(  -1.3411850570596516480D-2+
C19  BT*(  +9.0403349402888779200D-3+T*(  -2.4225953843706897600D-3+
C19  CT*(  +9.1578599714350784000D-4+T*(  -2.9689012151880000000D-4+
C19  DT*(  +1.0092815021080832000D-4+T*(  -3.3637584239226880000D-5+
C19  ET*(  +1.1252456510305280000D-5+T*(  -3.7549858152857600000D-6+
C19  FT*(  +1.2528422549504000000D-6+T*(  -4.1782339907584000000D-7+
C19  GT*(  +1.3931958370304000000D-7+T*(  -4.6445740523520000000D-8+
C19  HT*(  +1.5481318277120000000D-8+T*(  -5.1663077376000000000D-9+
C19  IT*(  +1.7255563264000000000D-9+T*( -5.6763875328000000000D-10+
C19  JT*( +1.8595971072000000000D-10+T*( -6.8985815040000000000D-11+
C19  KT*( +2.4998051840000000000D-11+T*( -4.1523609600000000000D-12+
C19  LT*( +6.7108864000000000000D-13+T*( -1.6777216000000000000D-12+
C19  MT*( +6.7108864000000000000D-13))))))))))))))))))))))))))
C
      DGAMM = Y*G
      GO TO 240
C
C     ERROR 3 TEST
  140 IF (T.LT.XMINV) GO TO 200
      DGAMM = 1.0D0/X
      WRITE(6,*) 'ARGUMENT SMALLER THAN THE SMALLEST MACHINE NUMBER'
      GO TO 240
C
C     ERROR EXITS
  160 IF (X.LT.0.0D0) GO TO 180
      DGAMM = GBIG
      WRITE(6,*) 'ARGUMENT LARGER THAN THE LARGEST MACHINE NUMBER'
      GO TO 240
C
  180 CONTINUE
      DGAMM = 0.0D0
      WRITE(6,*) 'ERROR 2 IN DGAMM FUNCTION'
      GO TO 240
C
  200 CONTINUE
      T = X
      IF (X.EQ.0.0D0) T = 1.0D0
      DGAMM = DSIGN(1.0D0/XMINV,T)
      WRITE(6,*) 'ERROR 3 IN DGAMM FUNCTION'
      GO TO 240
C
  220 CONTINUE
      DGAMM = GBIG
      WRITE(6,*) 'ARGUMENT LESS THEN ZERO'
C
  240 CONTINUE
      RETURN
      END
      DOUBLE PRECISION FUNCTION SHINT(YI,ZI,NPT1,NPT2,H1,H2,FX)
C
C-----------------------------------------------------------------------
C         SUBROUTINE TO DO A 2-D LINEAR INTERPOLATION
C         YI    -   Y VALUE OF THE INTERPOLATION POINT
C         ZI    -   Z VALUE OF THE INTERPOLATION POINT
C         NPT1  -   NUMBER OF POINTS IN Y DIRECTION
C         NPT2  -   NUMBER OF POINTS IN Z DIRECTION
C         H1    -   MESH SPACING IN Y DIRECTION
C         H2    -   MESH SPACING IN Z DIRECTION
C         FX    -   MATRIX OF (NPT1,NPT2) TO BE INTERPOLATED AT YI,ZI
C-----------------------------------------------------------------------
C
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION RX(2),FX(NPT1,NPT2)
      DATA   ZERO/0.0D0/
      SHINT = ZERO
      Y = YI
      Z = ZI
      RX(1) = Y
      RX(2) = Z
      X1 = RX(1)/H1
      X2 = RX(2)/H2
      I1 = X1 + 1
      I2 = X2 + 1
      X3 = (I1 - 1)*H1
      X4 = (I2 - 1)*H2
      P = (RX(1) - X3)/H1
      Q = (RX(2) - X4)/H2
      IF(I1.GE.NPT1.OR.I2.GE.NPT2) GO TO 300
      IF(I1.LT.1.OR.I2.LT.1) GO TO 300
      SHINT = (1.0D0 - P)*(1.0D0 - Q)*FX(I1,I2) + P*(1.0D0 - Q)*
     .FX(I1 + 1,I2) + Q*(1.D0 - P)*FX(I1,I2 + 1) + P*Q*FX(I1 + 1,I2 + 1)
 300  CONTINUE
      RETURN
      END
      SUBROUTINE GAUXW (NG2,NPTG,XG,WG)
C
C-------------------------------------------------------------------
C
C     THIS SUBROUTINE WILL RETURN ZEROS AND WEIGHTS OF THE
C     GAUSS INTEGRATION FORMULA BETWEEN (-1,+1). NG2 IS INPUTTED
C     SUCH THAT NG = 2**NG2 IS THE NUMBER OF INTEGRATION POINTS
C     (MAXIMUM NG2 = 8 , 256 POINTS).
C     NPTG=1 WILL PRINT OUT THE CHOSEN NUMBER OF ZEROS AND WEIGHTS.
C     ZEROS AND WEIGHTS ARE RETURNED IN XG,WG (XG(256),WG(256).
C
C-------------------------------------------------------------------
C
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION T(255),W(255),XG(256),WG(256)
C
C----WE HAVE ALL OF THE POSITIVE ZEROS AND WEIGHTS IN DATA STATEMENTS.
C
      DATA T(1),T(2),T(3),T(4),T(5),T(6),T(7),T(8),T(9),T(10),T(11),
     *     T(12),T(13),T(14),T(15)/.577350269189626D0,
     *  .861136311594053D0,.339981043584856D0,.960289856497536D0,
     *  .796666477413627D0,.525532409916329D0,.183434642495650D0,
     *  .989400934991650D0,.944575023073233D0,.865631202387832D0,
     *  .755404408355003D0,.617876244402644D0,.458016777657227D0,
     *  .281603550779259D0,.950125098376374D-1/
      DATA T(16),T(17),T(18),T(19),T(20),T(21),T(22),T(23),T(24),T(25),
     *     T(26),T(27),T(28),T(29),T(30),T(31)/.997263861849482D0,
     *  .985611511545268D0,.964762255587506D0,.934906075937740D0,
     *  .896321155766052D0,.849367613732570D0,.794483795967942D0,
     *  .732182118740290D0,.663044266930215D0,.587715757240762D0,
     *  .506899908932229D0,.421351276130635D0,.331868602282128D0,
     *  .239287362252137D0,.144471961582796D0,.483076656877383D-1/
      DATA T(32),T(33),T(34),T(35),T(36),T(37),T(38),T(39),T(40),T(41),
     *     T(42),T(43),T(44),T(45),T(46),T(47)/.999305041735772D0,
     *  .996340116771955D0,.991013371476744D0,.983336253884626D0,
     *  .973326827789911D0,.961008799652054D0,.946411374858403D0,
     *  .929569172131940D0,.910522137078503D0,.889315445995114D0,
     *  .865999398154093D0,.840629296252580D0,.813265315122798D0,
     *  .783972358943341D0,.752819907260532D0,.719881850171611D0/
      DATA T(48),T(49),T(50),T(51),T(52),T(53),T(54),T(55),T(56),T(57),
     *     T(58),T(59),T(60),T(61),T(62),T(63)/.685236313054233D0,
     *  .648965471254657D0,.611155355172393D0,.571895646202634D0,
     *  .531279464019895D0,.489403145707053D0,.446366017253464D0,
     *  .402270157963992D0,.357220158337668D0,.311322871990211D0,
     *  .264687162208767D0,.217423643740007D0,.169644420423993D0,
     *  .121462819296121D0,.729931217877990D-1,.243502926634244D-1/
      DATA T(64),T(65),T(66),T(67),T(68),T(69),T(70),T(71),T(72),T(73),
     *     T(74),T(75),T(76),T(77),T(78),T(79)/.999824887947132D0,
     *  .999077459977376D0,.997733248625514D0,.995792758534981D0,
     *  .993257112900213D0,.990127818491734D0,.986406742724586D0,
     *  .982096108435719D0,.977198491463907D0,.971716818747137D0,
     *  .965654366431965D0,.959014757853700D0,.951801961341264D0,
     *  .944020287830220D0,.935674388277916D0,.926769250878948D0/
      DATA T(80),T(81),T(82),T(83),T(84),T(85),T(86),T(87),T(88),T(89),
     *     T(90),T(91),T(92),T(93),T(94),T(95),T(96),T(97)/
     *  .917310198080961D0,.907302883401757D0,.896753288049158D0,
     *  .885667717345397D0,.874052796958032D0,.861915468939548D0,
     *  .849262987577969D0,.836102915060907D0,.822443116955644D0,
     *  .808291757507914D0,.793657294762193D0,.778548475506412D0,
     *  .762974330044095D0,.746944166797062D0,.730467566741909D0,
     *  .713554377683587D0,.696214708369514D0,.678458922447719D0/
      DATA T(98),T(99),T(100),T(101),T(102),T(103),T(104),T(105),T(106),
     *     T(107),T(108),T(109),T(110),T(111),T(112),T(113),T(114)/
     *  .660297632272646D0,.641741692562308D0,.622802193910585D0,
     *  .603490456158549D0,.583818021628763D0,.563796648226618D0,
     *  .543438302412810D0,.522755152051175D0,.501759559136144D0,
     *  .480464072404172D0,.458881419833552D0,.437024501037104D0,
     *  .414906379552275D0,.392540275033267D0,.369939555349859D0,
     *  .347117728597636D0,.324088435024413D0/
      DATA T(115),T(116),T(117),T(118),T(119),T(120),T(121),T(122),
     *     T(123),T(124),T(125),T(126),T(127)/.300865438877677D0,
     *  .277462620177904D0,.253893966422694D0,.230173564226660D0,
     *  .206315590902079D0,.182334305985337D0,.158244042714225D0,
     *  .134059199461188D0,.109794231127644D0,
     *  .854636405045155D-1,.610819696041396D-1,.366637909687335D-1,
     *  .122236989606158D-1/
      DATA T(128),T(129),T(130),T(131),T(132),T(133),T(134),T(135),
     *     T(136),T(137),T(138),T(139),T(140),T(141),T(142)/
     *  .999956050018992D0,.999768437409263D0,.999430937466261D0,
     *  .998943525843409D0,.998306266473006D0,.997519252756721D0,
     *  .996582602023382D0,.995496454481096D0,.994260972922410D0,
     *  .992876342608822D0,.991342771207583D0,.989660488745065D0,
     *  .987829747564861D0,.985850822286126D0,.983724009760315D0/
      DATA T(143),T(144),T(145),T(146),T(147),T(148),T(149),T(150),
     *     T(151),T(152),T(153),T(154),T(155),T(156),T(157)/
     *  .981449629025464D0,.979028021257622D0,.976459549719234D0,
     *  .973744599704370D0,.970883578480743D0,.967876915228489D0,
     *  .964725060975706D0,.961428488530732D0,.957987692411178D0,
     *  .954403188769716D0,.950675515316628D0,.946805231239127D0,
     *  .942792917117462D0,.938639174837814D0, .934344627502003D0/
      DATA T(158),T(159),T(160),T(161),T(162),T(163),T(164),T(165),
     *     T(166),T(167),T(168),T(169),T(170),T(171),T(172)/
     *  .929909919334006D0,.925335715583316D0,.920622702425146D0,
     *  .915771586857490D0,.910783096595065D0,.905657979960145D0,
     *  .900397005770304D0,.895000963223085D0,.889470661777611D0,
     *  .883806931033158D0,.878010620604707D0,.872082599995488D0,
     *  .866023758466555D0,.859835004903376D0,.853517267679503D0/
      DATA T(173),T(174),T(175),T(176),T(177),T(178),T(179),T(180),
     *     T(181),T(182),T(183),T(184),T(185),T(186),T(187)/
     *  .847071494517296D0,.840498652345763D0,.833799727155505D0,
     *  .826975723850813D0,.820027666098917D0,.812956596176432D0,
     *  .805763574812999D0,.798449681032171D0,.791016011989546D0,
     *  .783463682808184D0,.775793826411326D0,.768007593352446D0,
     *  .760106151642655D0,.752090686575492D0,.743962400549112D0/
      DATA T(188),T(189),T(190),T(191),T(192),T(193),T(194),T(195),
     *     T(196),T(197),T(198),T(199),T(200),T(201),T(202)/
     *  .735722512885918D0,.727372259649652D0,.718912893459971D0,
     *  .710345683304543D0,.701671914348685D0,.692892887742577D0,
     *  .684009920426076D0,.675024344931163D0,.665937509182049D0,
     *  .656750776292973D0,.647465524363725D0,.638083146272911D0,
     *  .628605049469015D0,.619032655759261D0,.609367401096334D0/
      DATA T(203),T(204),T(205),T(206),T(207),T(208),T(209),T(210),
     *     T(211),T(212),T(213),T(214),T(215),T(216),T(217)/
     *  .599610735362968D0,.589764122154454D0,.579829038559083D0,
     *  .569806974936569D0,.559699434694481D0,.549507934062719D0,
     *  .539234001866059D0,.528879179294822D0,.518445019673674D0,
     *  .507933088228616D0,.497344961852181D0,.486682228866890D0,
     *  .475946488786983D0,.465139352078479D0,.454262439917590D0/
      DATA T(218),T(219),T(220),T(221),T(222),T(223),T(224),T(225),
     *     T(226),T(227),T(228),T(229),T(230),T(231),T(232)/
     *  .443317383947527D0,.432305826033741D0,.421229418017624D0,
     *  .410089821468717D0,.398888707435459D0,.387627756194516D0,
     *  .376308656998716D0,.364933107823654D0,.353502815112970D0,
     *  .342019493522372D0,.330484865662417D0,.318900661840106D0,
     *  .307268619799319D0,.295590484460136D0,.283868007657082D0/
      DATA T(233),T(234),T(235),T(236),T(237),T(238),T(239),T(240),
     *     T(241),T(242),T(243),T(244),T(245),T(246),T(247)/
     *  .272102947876337D0,.260297069991943D0,.248452145001057D0,
     *  .236569949758284D0,.224652266709132D0,.212700883622626D0,
     *  .200717593323127D0,.188704193421389D0,.176662486044902D0,
     *  .164594277567554D0,.152501378338656D0,.140385602411376D0,
     *  .128248767270607D0,.116092693560333D0,.103919204810509D0/
      DATA   T(248),T(249),T(250),T(251),T(252),T(253),T(254),T(255)/
     *  .917301271635196D-1,.795272891002330D-1,.673125211657164D-1,
     *  .550876556946340D-1,.428545265363791D-1,.306149687799790D-1,
     *  .183708184788137D-1,.612391237518953D-2/
      DATA W(1),W(2),W(3),W(4),W(5),W(6),W(7),W(8),W(9),W(10),W(11),
     *     W(12),W(13),W(14),W(15)/1.0D0,.347854845137454D0,
     *  .652145154862546D0,.101228536290376D0,.222381034453374D0,
     *  .313706645877887D0,.362683783378362D0,.271524594117541D-1,
     *  .622535239386479D-1,.951585116824928D-1,.124628971255534D0,
     *  .149595988816577D0,.169156519395003D0,.182603415044924D0,
     *  .189450610455068D0/
      DATA W(16),W(17),W(18),W(19),W(20),W(21),W(22),W(23),W(24),W(25),
     *     W(26),W(27),W(28),W(29),W(30),W(31)/.701861000947010D-2,
     *  .162743947309057D-1,.253920653092621D-1,.342738629130214D-1,
     *  .428358980222267D-1,.509980592623762D-1,.586840934785355D-1,
     *  .658222227763618D-1,.723457941088485D-1,.781938957870703D-1,
     *  .833119242269468D-1,.876520930044038D-1,.911738786957639D-1,
     *  .938443990808046D-1,.956387200792749D-1,.965400885147278D-1/
      DATA W(32),W(33),W(34),W(35),W(36),W(37),W(38),W(39),W(40),W(41),
     *     W(42),W(43),W(44),W(45),W(46),W(47)/.178328072169643D-2,
     *  .414703326056247D-2,.650445796897836D-2,.884675982636395D-2,
     *  .111681394601311D-1,.134630478967186D-1,.157260304760247D-1,
     *  .179517157756973D-1,.201348231535302D-1,.222701738083833D-1,
     *  .243527025687109D-1,.263774697150547D-1,.283396726142595D-1,
     *  .302346570724025D-1,.320579283548516D-1,.338051618371416D-1/
      DATA W(48),W(49),W(50),W(51),W(52),W(53),W(54),W(55),W(56),W(57),
     *     W(58),W(59),W(60),W(61),W(62),W(63)/ .354722132568824D-1,
     *  .370551285402400D-1,.385501531786156D-1,.399537411327203D-1,
     *  .412625632426235D-1,.424735151236536D-1,.435837245293235D-1,
     *  .445905581637566D-1,.454916279274181D-1, .462847965813144D-1,
     *  .469681828162100D-1,.475401657148303D-1,.479993885964583D-1,
     *  .483447622348030D-1,.485754674415034D-1,.486909570091397D-1/
      DATA W(64),W(65),W(66),W(67),W(68),W(69),W(70),W(71),W(72),W(73),
     *     W(74),W(75),W(76),W(77),W(78),W(79)/ .449380960292090D-3,
     *  .104581267934035D-2,.164250301866903D-2,.223828843096262D-2,
     *  .283275147145799D-2,.342552604091022D-2,.401625498373864D-2,
     *  .460458425670296D-2,.519016183267633D-2,.577263754286570D-2,
     *  .635166316170719D-2,.692689256689881D-2,.749798192563473D-2,
     *  .806458989048606D-2,.862637779861675D-2,.918300987166087D-2/
      DATA W(80),W(81),W(82),W(83),W(84),W(85),W(86),W(87),W(88),
     *   W(89),W(90),W(91),W(92),W(93),W(94),W(95)/.973415341500681D-2,
     *  .102794790158322D-1,.108186607395031D-1,.113513763240804D-1,
     *  .118773073727403D-1,.123961395439509D-1,.129075627392673D-1,
     *  .134112712886163D-1,.139069641329520D-1,.143943450041668D-1,
     *  .148731226021473D-1,.153430107688651D-1,.158037286593993D-1,
     *  .162550009097852D-1,.166965578015892D-1,.171281354231114D-1/
      DATA W(96),W(97),W(98),W(99),W(100),W(101),W(102),W(103),W(104),
     *     W(105),W(106),W(107),W(108),W(109),W(110),W(111),W(112)/
     *  .175494758271177D-1,.179603271850087D-1,.183604439373313D-1,
     *  .187495869405447D-1,.191275236099509D-1,.194940280587066D-1,
     *  .198488812328309D-1,.201918710421300D-1,.205227924869601D-1,
     *  .208414477807511D-1,.211476464682213D-1,.214412055392085D-1,
     *  .217219495380521D-1,.219897106684605D-1,.222443288937998D-1,
     *  .224856520327450D-1,.227135358502365D-1/
      DATA W(113),W(114),W(115),W(116),W(117),W(118),W(119),W(120),
     *     W(121),W(122),W(123),W(124),W(125),W(126),W(127)/
     *  .229278441436868D-1,.231284488243870D-1,.233152299940628D-1,
     *  .234880760165359D-1,.236468835844476D-1,.237915577810034D-1,
     *  .239220121367035D-1,.240381686810241D-1,.241399579890193D-1,
     *  .242273192228152D-1,.243002001679719D-1,.243585572646906D-1,
     *  .244023556338496D-1,.244315690978500D-1,.244461801962625D-1/
      DATA W(128),W(129),W(130),W(131),W(132),W(133),W(134),W(135),
     *     W(136),W(137),W(138),W(139),W(140),W(141),W(142)/
     *  .112789017822272D-3,.262534944296446D-3,.412463254426176D-3,
     *  .562348954031410D-3,.712154163473321D-3,.861853701420089D-3,
     *  .101142439320844D-2,.116084355756772D-2,.131008868190250D-2,
     *  .145913733331073D-2,.160796713074933D-2,.175655573633073D-2,
     *  .190488085349972D-2,.205292022796614D-2,.220065164983991D-2/
      DATA W(143),W(144),W(145),W(146),W(147),W(148),W(149),W(150),
     *     W(151),W(152),W(153),W(154),W(155),W(156),W(157)/
     *  .234805295632731D-2,.249510203470371D-2,.264177682542749D-2,
     *  .278805532532771D-2,.293391559082972D-2,.307933574119934D-2,
     *  .322429396179420D-2,.336876850731555D-2,.351273770505631D-2,
     *  .365617995814250D-2,.379907374876626D-2,.394139764140883D-2,
     *  .408313028605267D-2,.422425042138154D-2,.436473687796806D-2/
      DATA W(158),W(159),W(160),W(161),W(162),W(163),W(164),W(165),
     *     W(166),W(167),W(168),W(169),W(170),W(171),W(172)/
     *  .450456858144790D-2,.464372455568006D-2,.478218392589269D-2,
     *  .491992592181387D-2,.505692988078684D-2,.519317525086928D-2,
     *  .532864159391593D-2,.546330858864431D-2,.559715603368291D-2,
     *  .573016385060144D-2,.586231208692265D-2,.599358091911534D-2,
     *  .612395065556793D-2,.625340173954240D-2,.638191475210788D-2/
      DATA W(173),W(174),W(175),W(176),W(177),W(178),W(179),W(180),
     *     W(181),W(182),W(183),W(184),W(185),W(186),W(187)/
     *  .650947041505366D-2,.663604959378107D-2,.676163330017380D-2,
     *  .688620269544632D-2,.700973909296982D-2,.713222396107539D-2,
     *  .725363892583391D-2,.737396577381235D-2,.749318645480588D-2,
     *  .761128308454566D-2,.772823794738156D-2,.784403349893971D-2,
     *  .795865236875435D-2,.807207736287350D-2,.818429146643827D-2/
      DATA W(188),W(189),W(190),W(191),W(192),W(193),W(194),W(195),
     *     W(196),W(197),W(198),W(199),W(200),W(201),W(202)/
     *  .829527784623523D-2,.840501985322154D-2,.851350102502249D-2,
     *  .862070508840101D-2,.872661596169881D-2,.883121775724875D-2,
     *  .893449478375821D-2,.903643154866287D-2,.913701276045081D-2,
     *  .923622333095630D-2,.933404837762327D-2,.943047322573775D-2,
     *  .952548341062928D-2,.961906467984073D-2,.971120299526628D-2/
      DATA W(203),W(204),W(205),W(206),W(207),W(208),W(209),W(210),
     *     W(211),W(212),W(213),W(214),W(215),W(216),W(217)/
     *  .980188453525733D-2,.989109569669583D-2,.997882309703491D-2,
     *  .100650535763064D-1,.101497741990949D-1,.102329722564782D-1,
     *  .103146352679340D-1,.103947509832117D-1,.104733073841704D-1,
     *  .105502926865815D-1,.106256953418966D-1,.106995040389798D-1,
     *  .107717077058046D-1,.108422955111148D-1,.109112568660490D-1/
      DATA W(218),W(219),W(220),W(221),W(222),W(223),W(224),W(225),
     *     W(226),W(227),W(228),W(229),W(230),W(231),W(232)/
     *  .109785814257296D-1,.110442590908139D-1,.111082800090098D-1,
     *  .111706345765534D-1,.112313134396497D-1,.112903074958755D-1,
     *  .113476078955455D-1,.114032060430392D-1,.114570935980906D-1,
     *  .115092624770395D-1,.115597048540436D-1,.116084131622531D-1,
     *  .116553800949452D-1,.117005986066207D-1,.117440619140606D-1/
      DATA W(233),W(234),W(235),W(236),W(237),W(238),W(239),W(240),
     *     W(241),W(242),W(243),W(244),W(245),W(246),W(247)/
     *  .117857634973434D-1,.118256971008240D-1,.118638567340711D-1,
     *  .119002366727665D-1,.119348314595636D-1,.119676359049059D-1,
     *  .119986450878058D-1,.120278543565826D-1,.120552593295601D-1,
     *  .120808558957245D-1,.121046402153405D-1,.121266087205273D-1,
     *  .121467581157945D-1,.121650853785355D-1,.121815877594818D-1/
      DATA W(248),W(249),W(250),W(251),W(252),W(253),W(254),W(255)/
     *  .121962627831147D-1,.122091082480372D-1,.122201222273040D-1,
     *  .122293030687103D-1,.122366493950402D-1,.122421601042728D-1,
     *  .122458343697479D-1,.122476716402898D-1/
C
C-------WE NOW STORE THE REQUIRED NUMBER OF ZEROS AND WEIGHTS .
C
         NG  =  2**NG2
         NN2  =  NG/2
         DO 10 I = 1,NN2
            NN2I  =  NN2 - 1 + I
            NNI  =  NG + 1 - I
            XG(I)  =   - T(NN2I)
            XG(NNI)  =   - XG(I)
            WG(I)  =  W(NN2I)
            WG(NNI)  =  WG(I)
   10    CONTINUE
C
C-------IF NPTG  =  1 PRINT THE ZEROS AND WEIGHTS .
C
      IF(NPTG .NE. 1) GO TO 20
      WRITE(6,100) NG
      WRITE(6,200) (IG,XG(IG),IG,WG(IG),IG = 1,NG)
 20   CONTINUE
C
 100  FORMAT(/' FOR ',I3,' POINT GAUSS INTEGRATION ZEROS AND WEIGHTS FOL
     .LOW:'//)
 200  FORMAT(2(1X,'X(',I3,') = ',D20.10,1X,'W(',I3,') = ',D20.10))
      RETURN
      END
      SUBROUTINE MATINC (A,NDM,N)
C
C---------------------------------------------------------------------
C
C        SUBROUTINE FOR THE INVERSION OF A GENERAL COMPLEX MATRIX
C        A  -   COMPLEX MATRIX TO BE INVERTED
C        NDM-   MAXIMUM DIMENSION OF A IN THE CALLING PROGRAM
C        N  -   ACTUAL DIMENSION OF THE MATRIX TO BE INVERTED
C
C        NOTE: MAXIMUM DIMENSION OF A MATRIX FOR INVERSION IS 1000x1000
C
C-----------------------------------------------------------------------
C
      COMPLEX*16 A,SWAP,T
      REAL*8 AMAX,TEMP
      DIMENSION A(NDM,NDM),PIVOT(1000),INDEX(1000)
C
C     INITIALIZE AND PIVOT ELEMENT ARRAY
C
      DO 20 I=1,N
      PIVOT(I)=0.0
   20 CONTINUE
C
C     PERFORM SUCCESSIVE PIVOT OPERATIONS (GRAND LOOP)
C
      DO 550 I=1,N
C
C     SEARCH FOR PIVOT ELEMENT
C
      AMAX=0.0
      DO 105 J=1,N
      IF (PIVOT(J).NE.0.0) GO TO 105
      DO 100 K=1,N
      IF (PIVOT(K).NE.0.0) GO TO 100
      TEMP=CDABS(A(J,K))
      IF (TEMP.LT.AMAX) GO TO 100
      IROW=J
      ICOLUM=K
      AMAX=TEMP
  100 CONTINUE
  105 CONTINUE
      INDEX(I)=4096*IROW+ICOLUM
      J=IROW
      T=A(J,ICOLUM)
      PIVOT(ICOLUM)=AMAX
C
C     INTERCHANGE ROWS TO PUT PIVOT ELEMENT ON DIAGONAL
C
      IF (IROW.EQ.ICOLUM) GO TO 260
      DO 200 K=1,N
      SWAP=A(J,K)
      A(J,K)=A(ICOLUM,K)
      A(ICOLUM,K)=SWAP
  200 CONTINUE
C
C     DIVIDE PIVOT ROW BY PIVOT ELEMENT
C
  260 DO 350 K=1,N
      IF (K.EQ.ICOLUM) A(ICOLUM,K)=1.0
      A(ICOLUM,K)=A(ICOLUM,K)/T
  350 CONTINUE
C
C     REDUCE NON-PIVOT ROWS
C
      DO 550 J=1,N
      IF (J.EQ.ICOLUM) GO TO 550
      T=A( J,ICOLUM)
      A( J,ICOLUM)=0.0
      DO 450 K=1,N
      A( J,K)=A( J,K)-A(ICOLUM,K)*T
  450 CONTINUE
  550 CONTINUE
C
C     INTERCHANGE COLUMNS AFTER ALL PIVOT OPERATIONS HAVE BEEN PERFORMED
C
      DO 710 I=1,N
      I1=N+1-I
      K=INDEX(I1)/4096
      ICOLUM=INDEX(I1)-4096*K
      IF (K.EQ.ICOLUM) GO TO 710
      DO 705 J=1,N
      SWAP=A(J,K)
      A(J,K)=A(J,ICOLUM)
      A(J,ICOLUM)=SWAP
  705 CONTINUE
  710 CONTINUE
      RETURN
      END
      SUBROUTINE ZREAL1 (F,EPS,EPS2,ETA,NSIG,N,X,ITMAX)                 
C                                                                       
C-----------------------------------------------------------------------
C                                                                       
C   PURPOSE             - THE REAL ZEROS OF A REAL FUNCTION - TO BE     
C                           USED WHEN INITIAL GUESSES ARE POOR          
C                                                                       
C   USAGE               - CALL ZREAL1 (F,EPS,EPS2,ETA,NSIG,N,X,ITMAX)   
C                                                                       
C   ARGUMENTS    F      - A SINGLE-ARGUMENT REAL FUNCTION SUBPROGRAM    
C                           SUPPLIED BY THE USER. (INPUT)               
C                           F MUST BE DECLARED EXTERNAL IN THE CALLING  
C                           PROGRAM. F DEFINES THE FUNCTION FOR WHICH   
C                           THE ROOTS ARE TO BE FOUND.                  
C                EPS    - CONVERGENCE CRITERION. (INPUT)                
C                           A ROOT, X(I), IS ACCEPTED IF                
C                           ABS(F(X(I)) .LE. EPS.                       
C                EPS2   - SPREAD CRITERIA FOR MULTIPLE ROOTS. (INPUT)   
C                ETA        IF THE ROOT X(I) HAS BEEN COMPUTED AND IT IS
C                           FOUND THAT                                  
C                           ABS(X(I)-X(J)) .LT. EPS2                    
C                           WHERE X(J) IS A PREVIOUSLY COMPUTED ROOT,   
C                           THEN THE COMPUTATION IS RESTARTED WITH A    
C                           GUESS EQUAL TO X(I) + ETA.                  
C                NSIG   - CONVERGENCE CRITERION. (INPUT)                
C                           A ROOT IS ACCEPTED IF TWO SUCCESSIVE        
C                           APPROXIMATIONS TO A GIVEN ROOT AGREE        
C                           IN THE FIRST NSIG DIGITS.                   
C                         NOTE THAT IF EITHER CONVERGENCE CRITERION     
C                         IS SATISFIED, THE ROOT IS ACCEPTED.           
C                N      - THE NUMBER OF ROOTS TO BE FOUND. (INPUT)      
C                X      - VECTOR OF LENGTH N. (INPUT/OUTPUT)            
C                           ON INPUT, X CONTAINS THE INITIAL GUESSES    
C                           FOR THE ROOTS.                              
C                           ON OUTPUT, X CONTAINS THE COMPUTED ROOTS.   
C                ITMAX  - ITERATION INDICATOR. (INPUT/OUTPUT)           
C                           ON INPUT, ITMAX IS THE MAXIMUM NUMBER OF    
C                           ITERATIONS TO BE TAKEN PER ROOT.            
C                           ON OUTPUT, ITMAX IS THE NUMBER OF ITERATIONS
C                           USED IN FINDING THE LAST ROOT.              
C                                                                       
C   REMARKS  1.  ZREAL1 ASSUMES THAT THERE EXIST N DISTINCT REAL ROOTS  
C                FOR THE FUNCTION F AND THAT THEY CAN BE REACHED FROM   
C                THE INITIAL GUESSES SUPPLIED. THE ROUTINE IS DESIGNED  
C                SO THAT CONVERGENCE TO ANY SINGLE ROOT CANNOT BE       
C                OBTAINED FROM TWO DIFFERENT INITIAL GUESSES.           
C            2.  SCALING THE X VECTOR IN THE FUNCTION F MAY BE REQUIRED 
C                IF ANY OF THE ROOTS ARE KNOWN TO BE LESS THAN ONE.     
C                                                                       
C-----------------------------------------------------------------------
C                                                                       
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION X(N)                                                    
      DATA               TEN,ONE,ZERO,P9/10.0D0,1.0D0,0.0D0,0.90D0/     
      DATA               P11,HALF,PP1,F4/1.10D0,0.50D0,0.10D0,4.0D0/    
C                                                                       
      DIGT = TEN**(-NSIG)                                               
      P = -ONE                                                          
      P1 = ONE                                                          
      P2 = ZERO                                                         
      H = ZERO                                                          
      DO 95 L = 1,N                                                     
         JK = 0                                                         
         IF (X(L) .EQ. ZERO) GO TO 5                                    
         P = P9*X(L)                                                    
         P1 = P11*X(L)                                                  
         P2 = X(L)                                                      
    5    RT = P                                                         
         GO TO 65                                                       
   10    IF (JK .NE. 1) GO TO 15                                        
         RT = P1                                                        
         X0 = FPRT                                                      
         GO TO 65                                                       
   15    IF (JK .NE. 2) GO TO 20                                        
         RT = P2                                                        
         X1 = FPRT                                                      
         GO TO 65                                                       
   20    IF (JK .NE. 3) GO TO 55                                        
         X2 = FPRT                                                      
         D = -HALF                                                      
         IF (X(L) .EQ. ZERO) GO TO 25                                   
         H =-PP1*X(L)                                                   
         GO TO 30                                                       
   25    H = -ONE                                                       
   30    DD = ONE+D                                                     
         BI = X0*D**2-X1*DD**2+X2*(DD+D)                                
         DEN = BI**2 -F4*X2*D*DD*(X0*D-(X1*DD)+X2)                      
         IF (DEN .LE. ZERO) GO TO 35                                    
         DEN =DSQRT(DEN)                                                
         GO TO 40                                                       
   35    DEN = ZERO                                                     
   40    DN = BI + DEN                                                  
         DM = BI - DEN                                                  
         IF (DABS(DN) .LE. DABS(DM)) GO TO 45                           
         DEN = DN                                                       
         GO TO 50                                                       
   45    DEN = DM                                                       
   50    IF (DEN .EQ. ZERO) DEN = ONE                                   
         DI=-DD*(X2+X2)/DEN                                             
         H = DI * H                                                     
         RT = RT + H                                                    
C                                  TEST FOR CONVERGENCE                 
         IF (DABS(H) .LT. DABS(RT)*DIGT) GO TO 90                       
         GO TO 65                                                       
   55    IF (DABS(FPRT) .GE. DABS(X2*10.0))  GO TO 60                   
         X0 = X1                                                        
         X1 = X2                                                        
         X2 = FPRT                                                      
         D = DI                                                         
         GO TO 30                                                       
   60    DI = DI * HALF                                                 
         H = H * HALF                                                   
         RT = RT - H                                                    
   65    JK = JK + 1                                                    
         IF (JK .LT. ITMAX)  GO TO 75                                   
C                                  WARNING  ERROR ITERATIONS = MAXIMUM  
         WRITE(6,*) ' ZREAL1: NO CONVERGENCE IN ITMAX ITERATIONS'       
         WRITE(6,*) ' ROOT ',L,'  IS SET TO 111111.0'
         X(L)=111111.                                                   
         GO TO 95                                                       
   75    FRT = F(RT)                                                    
         FPRT = FRT                                                     
         IF (L .LT. 2) GO TO 81                                         
         DO 80 I = 2,L                                                  
            TEM = RT - X(I-1)                                           
            IF (DABS(TEM) .LT. EPS2)  GO TO 85                          
            FPRT = FPRT/TEM                                             
   80    CONTINUE                                                       
C                                  TEST FOR CONVERGENCE                 
   81    IF ((DABS(FRT) .LT. EPS) .AND. (DABS(FPRT) .LT. EPS))  GO TO 90
         GO TO 10                                                       
   85    RT = RT + ETA                                                  
         JK = JK - 1                                                    
         GO TO 65                                                       
   90    X(L) = RT                                                      
   95 CONTINUE                                                          
      ITMAX = JK                                                        
      RETURN                                                            
      END                                                               
      SUBROUTINE EIGRF (A,N,IA,IJOB,W,Z,IZ,WK)
C
C-----------------------------------------------------------------------
C
C   PURPOSE             - EIGENVALUES AND (OPTIONALLY) EIGENVECTORS OF
C                           A REAL GENERAL MATRIX IN FULL STORAGE MODE
C
C   USAGE               - CALL EIGRF (A,N,IA,IJOB,W,Z,IZ,WK)
C
C   ARGUMENTS    A      - THE INPUT REAL GENERAL MATRIX OF ORDER N
C                           WHOSE EIGENVALUES AND EIGENVECTORS ARE
C                           TO BE COMPUTED. INPUT A IS DESTROYED IF
C                           IJOB IS EQUAL TO 0 OR 1.
C                N      - THE INPUT ORDER OF THE MATRIX A.
C                IA     - THE INPUT ROW DIMENSION OF MATRIX A EXACTLY
C                           AS SPECIFIED IN THE DIMENSION STATEMENT IN
C                           THE CALLING PROGRAM.
C                IJOB   - THE INPUT OPTION PARAMETER. WHEN
C                           IJOB = 0, COMPUTE EIGENVALUES ONLY
C                           IJOB = 1, COMPUTE EIGENVALUES AND EIGEN-
C                             VECTORS.
C                           IJOB = 2, COMPUTE EIGENVALUES, EIGENVECTORS
C                             AND PERFORMANCE INDEX.
C                           IJOB = 3, COMPUTE PERFORMANCE INDEX ONLY.
C                           IF THE PERFORMANCE INDEX IS COMPUTED, IT IS
C                           RETURNED IN WK(1). THE ROUTINES HAVE
C                           PERFORMED (WELL, SATISFACTORILY, POORLY) IF
C                           WK(1) IS (LESS THAN 1, BETWEEN 1 AND 100,
C                           GREATER THAN 100).
C                W      - THE OUTPUT COMPLEX VECTOR OF LENGTH N,
C                           CONTAINING THE EIGENVALUES OF A.
C                         NOTE - THE ROUTINE TREATS W AS A REAL VECTOR
C                           OF LENGTH 2*N. AN APPROPRIATE EQUIVALENCE
C                           STATEMENT MAY BE REQUIRED. SEE DOCUMENT
C                           EXAMPLE.
C                Z      - THE OUTPUT N BY N COMPLEX MATRIX CONTAINING
C                           THE EIGENVECTORS OF A.
C                           THE EIGENVECTOR IN COLUMN J OF Z CORRES-
C                           PONDS TO THE EIGENVALUE W(J).
C                           IF IJOB = 0, Z IS NOT USED.
C                         NOTE - THE ROUTINE TREATS Z AS A REAL VECTOR
C                           OF LENGTH 2*N*N. AN APPROPRIATE EQUIVALENCE
C                           STATEMENT MAY BE REQUIRED. SEE DOCUMENT
C                           EXAMPLE.
C                IZ     - THE INPUT ROW DIMENSION OF MATRIX Z EXACTLY
C                           AS SPECIFIED IN THE DIMENSION STATEMENT IN
C                           THE CALLING PROGRAM. IZ MUST BE GREATER
C                           THAN OR EQUAL TO N IF IJOB IS NOT EQUAL TO
C                           ZERO.
C                WK     - WORK AREA, THE LENGTH OF WK DEPENDS
C                           ON THE VALUE OF IJOB, WHEN
C                           IJOB = 0, THE LENGTH OF WK IS AT LEAST N.
C                           IJOB = 1, THE LENGTH OF WK IS AT LEAST 2N.
C                           IJOB = 2, THE LENGTH OF WK IS AT LEAST
C                             (2+N)N.
C                           IJOB = 3, THE LENGTH OF WK IS AT LEAST 1.
C
C-----------------------------------------------------------------------
C
C                                  SPECIFICATIONS FOR ARGUMENTS
      INTEGER            N,IA,IJOB,IZ,IER
      DOUBLE PRECISION   A(IA,1),WK(N,1),W(1),Z(1)
C                                  SPECIFICATIONS FOR LOCAL VARIABLES
      INTEGER            JER,IZ2,K,L,I,N1,N2,II,JJ,NP1,IIZ,NPI,JW,J,
     *                   IS,IG,IGZ,LW,LLZ,KKZ,LZ,KZ
      DOUBLE PRECISION   ANORM,ASUM,PI,SUMZ,SUMR,SUMI,S,TEN,RDELP,
     *                   ZERO,ONE,THOUS,AN,Z11
      DATA               RDELP/0.222045D-15/
      DATA               ZERO,ONE/0.0D0,1.0D0/,TEN/10.0D0/,THOUS/1000.0D
     *0/
C                                  INITIALIZE ERROR PARAMETERS
C                                  FIRST EXECUTABLE STATEMENT
      JER = 0
      IZ2 = IZ+IZ
      IF (IJOB .GE. 0 .AND. IJOB .LE. 3) GO TO 5
C                                  WARNING ERROR - IJOB IS NOT IN THE
C                                    RANGE
      WRITE(6,*) ' EIGRF: IJOB NOT IN RANGE, IJOB=1'
      IJOB = 1
      GO TO 10
    5 IF (IJOB .EQ. 0) GO TO 16
   10 IF (IZ .GE. N) GO TO 15
C                                  WARNING ERROR - IZ IS LESS THAN N
C                                    EIGENVECTORS CAN NOT BE COMPUTED,
C                                    IJOB SET TO ZERO
      WRITE(6,*) ' EIGRF: IZ .LT. N. IJOB SET TO ZERO'
      IJOB = 0
   15 IF (IJOB .EQ. 3) GO TO 95
C                                  PACK A INTO AN N BY N ARRAY
   16 K = 1
      L = 1
      DO 20 J=1,N
         DO 20 I=1,N
            A(K,L) = A(I,J)
C                                  SAVE INPUT A IF IJOB = 2
            IF (IJOB .EQ. 2) WK(I,J)=A(I,J)
            K = K+1
            IF (K .GT. IA) K = 1
            IF (K .EQ. 1) L = L+1
   20 CONTINUE
      N1 = 1
      IF (IJOB .EQ. 2) N1 = N+1
      N2 = N1+1
      IF (IJOB .EQ. 0) N2 = 1
C                                  BALANCE THE INPUT A
      CALL EBALAF (A,N,N,WK(1,N1),K,L)
      IF (IJOB .EQ. 0 .AND. L .EQ. 0) GO TO 35
C                                  IF L = 0, A IS ALREADY IN HESSENBERG
C                                    FORM
      CALL EHESSF (A,K,L,N,N,WK(1,N2))
      IF (IJOB .EQ. 0) GO TO 35
C                                  SET Z IDENTITY MATRIX
      II = 1
      JJ = 1
      NP1 = N+1
      DO 30 I=1,N
         DO 25 J=1,N
            Z(II) = ZERO
            II = II+1
   25    CONTINUE
         Z(JJ) = ONE
         JJ = JJ+NP1
   30 CONTINUE
      CALL EHBCKF (Z,A,WK(1,N2),N,N,N,K,L)
      IIZ = N
   35 IF (IJOB .EQ. 0) IIZ = 1
      IF(IJOB .EQ. 0 .AND. N .EQ. 1) Z11 = Z(1)
      CALL EQRH3F (A,N,N,K,L,W(1),W(N+1),Z,IIZ,JER)
      IF(IJOB .EQ. 0 .AND. N .EQ. 1) Z(1) = Z11
      IF (JER .GT. 0 .OR. IJOB .EQ. 0) GO TO 40
      CALL EBBCKF (WK(1,N1),Z,K,L,N,N,N)
C                                  CONVERT W (EIGENVALUES) TO COMPLEX
C                                    FORMAT
   40 DO 45 I=1,N
         NPI = N+I
         WK(I,N1) = W(NPI)
   45 CONTINUE
      JW = N+N
      J = N
      DO 50 I=1,N
         W(JW-1) = W(J)
         W(JW) = WK(J,N1)
         JW = JW-2
         J = J-1
   50 CONTINUE
      IF (IJOB .EQ. 0) GO TO 9000
C                                  CONVERT Z (EIGENVECTORS) TO COMPLEX
C                                    FORMAT Z(IZ,N)
      J = N
   60 IF (J .LT. 1) GO TO 85
      IF (W(J+J) .EQ. ZERO) GO TO 75
C                                  MOVE PAIR OF COMPLEX CONJUGATE
C                                    EIGENVECTORS
      IS = IZ2*(J-1)+1
      IG = N*(J-2)+1
      IGZ = IG+N
C                                  MOVE COMPLEX CONJUGATE EIGENVECTOR
      DO 65 I=1,N
         Z(IS) = Z(IG)
         Z(IS+1) = -Z(IGZ)
         IS = IS+2
         IG = IG+1
         IGZ = IGZ+1
   65 CONTINUE
C                                  MOVE COMPLEX EIGENVECTOR
      IS = IZ2*(J-2)+1
      IG = IS+IZ2
      DO 70 I=1,N
         Z(IS) = Z(IG)
         Z(IS+1) = -Z(IG+1)
         IS = IS+2
         IG = IG+2
   70 CONTINUE
      J = J-2
      GO TO 60
C                                  MOVE REAL EIGENVECTOR
   75 IS = IZ2*(J-1)+N+N
      IG = N*J
      DO 80 I=1,N
         Z(IS-1) = Z(IG)
         Z(IS) = ZERO
         IS = IS-2
         IG = IG-1
   80 CONTINUE
      J = J-1
      GO TO 60
C                                  Z IS NOW IN COMPLEX FORMAT Z(IZ,N).
C                                    NEXT, MOVE ORIGINAL MATRIX BACK
C                                    TO A
   85 IF (IJOB .LE. 1) GO TO 9000
      DO 90 I=1,N
         DO 90 J=1,N
            A(I,J) = WK(I,J)
   90 CONTINUE
      WK(1,1) = THOUS
      IF (JER .NE. 0) GO TO 9000
C                                  COMPUTE 1-NORM OF A
   95 ANORM = ZERO
      DO 105 J=1,N
         ASUM = ZERO
         DO 100 I=1,N
            ASUM = ASUM+DABS(A(I,J))
  100    CONTINUE
         ANORM = DMAX1(ANORM,ASUM)
  105 CONTINUE
      IF (ANORM .EQ. ZERO) ANORM = ONE
C                                  COMPUTE PERFORMANCE INDEX
      PI = ZERO
      LW = 1
      LLZ = 0
      KKZ = 0
      DO 120 J=1,N
         S = ZERO
         SUMZ = ZERO
         LZ = LLZ+1
         KZ = KKZ+1
         LW = J+J-1
         DO 115 L=1,N
            SUMZ = SUMZ+CDABS(DCMPLX(Z(LZ),Z(LZ+1)))
            KZ = KKZ+1
            SUMR = -W(LW)*Z(LZ)+W(LW+1)*Z(LZ+1)
            SUMI = -W(LW)*Z(LZ+1)-W(LW+1)*Z(LZ)
            DO 110 K=1,N
               SUMR =SUMR+A(L,K)*Z(KZ)
               SUMI = SUMI+A(L,K)*Z(KZ+1)
               KZ = KZ+2
  110       CONTINUE
            S = S+CDABS(DCMPLX(SUMR,SUMI))
            LZ = LZ+2
  115    CONTINUE
         PI = DMAX1(PI,S/SUMZ)
         KKZ = KKZ+IZ2
         LLZ = LLZ+IZ2
  120 CONTINUE
      AN = N
      PI = PI/(ANORM*TEN*AN*RDELP)
      WK(1,1) = PI
 9000 CONTINUE
      IF (JER .EQ. 0) GO TO 9005
      WRITE(6,*) ' EIGRF: EIGENVALUES .LT. ',JER,'TH EIGENVALUE DID NOT
     XCONVERGE'
 9005 RETURN
      END
      SUBROUTINE EBALAF (A,N,IA,D,K,L)
C                                  SPECIFICATIONS FOR ARGUMENTS
      INTEGER            N,IA,K,L
      DOUBLE PRECISION   A(IA,1),D(1)
C                                  SPECIFICATIONS FOR LOCAL VARIABLES
      INTEGER            L1,K1,K1P1,K11,JJ,J,I,LL,NOCONV
      DOUBLE PRECISION   R,C,F,G,B,S,B2,ONE,ZERO,P95
      DATA               B/16.0D0/,B2/256.0D0/
      DATA               ZERO/0.0D0/,ONE/1.0D0/,P95/.95D0/
C                                  REDUCE NORM A BY DIAGONAL SIMILARITY
C                                  TRANSFORMATION STORED IN D
C                                  FIRST EXECUTABLE STATEMENT
      L1 = 1
      K1 = N
C                                  SEARCH FOR ROWS ISOLATING AN EIGEN-
C                                    VALUE AND PUSH THEM DOWN
    5 K1P1 = K1+1
      IF (K1.LT.1) GO TO 35
      K11=K1
      DO 30 JJ=1,K11
         J = K1P1-JJ
         R = ZERO
         DO 10 I=1,K1
            IF (I.EQ.J) GO TO 10
            R=R+DABS(A(J,I))
   10    CONTINUE
         IF (R.NE.ZERO) GO TO 30
         D(K1) = J
         IF (J.EQ.K1) GO TO 25
         DO 15 I=1,K1
            F = A(I,J)
            A(I,J) = A(I,K1)
            A(I,K1) = F
   15    CONTINUE
         DO 20 I=L1,N
            F = A(J,I)
            A(J,I) = A(K1,I)
            A(K1,I) = F
   20    CONTINUE
   25    K1 = K1-1
         GO TO 5
   30 CONTINUE
C                                  SEARCH FOR COLUMNS ISOLATING AN
C                                    EIGENVALUE AND PUSH THEM LEFT
   35 IF (K1.LT.L1) GO TO 65
      LL = L1
      DO 60 J=LL,K1
         C = ZERO
         DO 40 I=L1,K1
            IF (I.EQ.J) GO TO 40
            C = C+DABS(A(I,J))
   40    CONTINUE
         IF (C.NE.ZERO) GO TO 60
         D(L1) = J
         IF (J.EQ.L1) GO TO 55
         DO 45 I=1,K1
            F = A(I,J)
            A(I,J) = A(I,L1)
            A(I,L1) = F
   45    CONTINUE
         DO 50  I=L1,N
            F = A(J,I)
            A(J,I) = A(L1,I)
            A(L1,I) = F
   50    CONTINUE
   55    L1 = L1+1
         GO TO 35
   60 CONTINUE
C                                  NOW BALANCE THE SUBMATRIX IN ROWS
C                                    L1 THROUGH K1
   65 K = L1
      L = K1
      IF (K1.LT.L1) GO TO 75
      DO 70  I=L1,K1
         D(I) = ONE
   70 CONTINUE
   75 NOCONV = 0
      IF (K1.LT.L1) GO TO 120
      DO 115 I=L1,K1
         C = ZERO
         R = ZERO
         DO 80 J=L1,K1
            IF (J.EQ.I) GO TO 80
            C = C+DABS(A(J,I))
            R = R+DABS(A(I,J))
   80    CONTINUE
         G = R/B
         F = ONE
         S = C+R
   85    IF (C.GE.G) GO TO 90
         F = F * B
         C = C*B2
         GO TO 85
   90    G = R*B
   95    IF (C.LT.G) GO TO 100
         F = F/B
         C = C/B2
         GO TO 95
C                                  NOW BALANCE
  100    IF ((C+R)/F.GE.P95*S) GO TO 115
         G = ONE/F
         D(I) = D(I)*F
         NOCONV = 1
         DO 105 J=L1,N
            A(I,J) = A(I,J)*G
  105    CONTINUE
         DO 110 J=1,K1
            A(J,I) = A(J,I)*F
  110    CONTINUE
  115 CONTINUE
  120 IF (NOCONV.EQ.1) GO TO 75
      RETURN
      END
      SUBROUTINE EBBCKF (D,Z,K,L,MM,N,IZ)
C                                  SPECIFICATIONS FOR ARGUMENTS
      INTEGER            K,L,MM,N,IZ
      DOUBLE PRECISION   Z(IZ,1),D(1)
C                                  SPECIFICATIONS FOR LOCAL VARIABLES
      INTEGER            I,J,KM1,II,JJ,LP1
      DOUBLE PRECISION   S
C                                  COLUMN SCALE Z BY APPROPRIATE D VALUE
C                                  FIRST EXECUTABLE STATEMENT
      IF (L.EQ.0) GO TO 15
      DO 10 I=K,L
         S = D(I)
         DO 5 J=1,MM
            Z(I,J) = Z(I,J)*S
    5    CONTINUE
   10 CONTINUE
C                                  INTERCHANGE ROWS IF PERMUTATIONS
C                                    OCCURRED IN EBALAF
   15 IF (K.EQ.1) GO TO 30
      KM1 = K-1
      DO 25 I=1,KM1
         II = K-I
         JJ = D(II)
         IF (II.EQ.JJ) GO TO 25
         DO 20 J=1,MM
            S = Z(II,J)
            Z(II,J) = Z(JJ,J)
            Z(JJ,J) = S
   20    CONTINUE
   25 CONTINUE
   30 IF (L.EQ.N) GO TO 45
      LP1 = L+1
      DO 40 II=LP1,N
         JJ = D(II)
         IF (II.EQ.JJ) GO TO 40
         DO 35 J=1,MM
            S = Z(II,J)
            Z(II,J) = Z(JJ,J)
            Z(JJ,J) = S
   35    CONTINUE
   40 CONTINUE
   45 RETURN
      END
      SUBROUTINE EHBCKF (Z,H,D,N,MM,IZH,K,L)
C                                  SPECIFICATIONS FOR ARGUMENTS
      INTEGER            N,MM,IZH,K,L
      DOUBLE PRECISION   Z(IZH,1),H(IZH,1),D(1)
C                                  SPECIFICATIONS FOR LOCAL VARIABLES
      INTEGER            LM2,KI,LTEMP,M,MA,MP2,I,J
      DOUBLE PRECISION   G,ZERO
      DATA               ZERO/0.0D0/
C                                  ADAPTED FROM EISPACK ROUTINE ORTBAK
C                                  FIRST EXECUTABLE STATEMENT
      LM2=L-2
      IF(LM2.LT.K) GO TO 9005
      LTEMP=LM2+K
      DO 30 KI=K,LM2
         M=LTEMP-KI
         MA=M+1
         IF(H(MA,M).EQ.ZERO) GO TO 30
         MP2=M+2
         IF(MP2.GT.L) GO TO 10
         DO 5 I=MP2,L
            D(I)=H(I,M)
    5    CONTINUE
   10    IF(MA.GT.L) GO TO 30
         DO 25 J=1,MM
            G=ZERO
            DO 15 I=MA,L
               G=G+D(I)*Z(I,J)
   15       CONTINUE
C                                  DOUBLE DIVISION AVOIDS POSSIBLE
C                                  UNDERFLOW
            G = (G/D(MA))/H(MA,M)
            DO 20 I=MA,L
               Z(I,J)=Z(I,J)+G*D(I)
   20       CONTINUE
   25    CONTINUE
   30 CONTINUE
 9005 RETURN
      END
      SUBROUTINE EHESSF (A,K,L,N,IA,D)
C                                  SPECIFICATIONS FOR ARGUMENTS
      INTEGER            K,L,N,IA
      DOUBLE PRECISION   A(IA,N),D(N)
C                                  SPECIFICATIONS FOR LOCAL VARIABLES
      INTEGER            LA,KP1,M,I,MP,II,J,JJ
      DOUBLE PRECISION   F,G,H,SCALE,ZERO
      DATA               ZERO/0.0D0/
C                                  FIRST EXECUTABLE STATEMENT
      LA = L - 1
      KP1 = K + 1
      IF (LA .LT. KP1) GO TO 50
      DO 45 M = KP1, LA
         H = ZERO
         D(M) = ZERO
         SCALE = ZERO
C                                  SCALE COLUMN
         DO 5 I = M, L
            SCALE = SCALE + DABS(A(I,M-1))
    5    CONTINUE
         IF (SCALE .EQ. ZERO ) GO TO 45
         MP = M + L
C                                  DO 10 I=L,M,-1
         DO 10 II = M, L
            I = MP - II
            D(I) = A(I,M-1) / SCALE
            H = H + D(I) * D(I)
   10    CONTINUE
         G = -DSIGN(DSQRT(H),D(M))
         H = H - D(M) * G
         D(M) = D(M) - G
C                                  FORM (I-(U*UT)/H) * A
         DO 25 J = M,N
            F = ZERO
C                                  DO 15 I=L,M,-1
            DO 15 II = M, L
               I = MP - II
               F = F + D(I) * A(I,J)
   15       CONTINUE
            F = F / H
            DO 20 I = M, L
               A(I,J) = A(I,J) - F * D(I)
   20       CONTINUE
   25    CONTINUE
C                                  FORM (I-(U*UT)/H)*A*(I-(U*UT)/H)
         DO 40 I = 1,L
            F = ZERO
C                                  DO 30 J=L,M,-1
            DO 30 JJ = M, L
               J = MP - JJ
               F = F + D(J) * A(I,J)
   30       CONTINUE
            F = F / H
            DO 35 J = M, L
               A(I,J) = A(I,J) - F * D(J)
   35       CONTINUE
   40    CONTINUE
         D(M) = SCALE * D(M)
         A(M,M-1) = SCALE * G
   45 CONTINUE
   50 RETURN
      END
      SUBROUTINE EQRH3F (H,N,IH,K,L,WR,WI,Z,IZ,IER)
C                                  SPECIFICATIONS FOR ARGUMENTS
      INTEGER            N,IH,K,L,IZ,IER
      DOUBLE PRECISION   H(IH,N),WR(N),WI(N),Z(IZ,N)
C                                  SPECIFICATIONS FOR LOCAL VARIABLES
      INTEGER            I,IEN,ITS,IENM2,NPL,LL,LB,NAML,MM,M,MP2,KA,NA,
     *                   J,JJ
      DOUBLE PRECISION   T3(2),RDELP,P4,P5,P7,ZERO,ONE,T,X,Y,W,S,ZZ,R,P,
     *                   Q,RNORM,RA,SA,SCALE,VR,VI
      COMPLEX*16         Z3
      LOGICAL            NOTLAS
      EQUIVALENCE        (Z3,T3(1))
      DATA               RDELP/0.222045D-15/
      DATA               P4 /0.4375D0/,P5 /0.5D0/,P7 /0.75D0/,ZERO /0.0D
     *0/,ONE
     *                   /1.0D0/
C                                  FIRST EXECUTABLE STATEMENT
      IER = 0
C                                  STORE ROOTS ISOLATED BY EBALAF
      RNORM = 0.0D0
      KA = 1
      DO 10 I=1,N
         DO 5 J=KA,N
    5    RNORM = RNORM+DABS(H(I,J))
         KA = I
         IF (I.GE.K .AND. I.LE.L) GO TO 10
         WR(I) = H(I,I)
         WI(I) = ZERO
   10 CONTINUE
      IEN = L
      T = ZERO
C                                  SEARCH FOR NEXT EIGENVALUES
   15 IF (IEN.LT.K) GO TO 145
      ITS = 0
      NA = IEN-1
      IENM2 = NA-1
C                                  LOOK FOR SINGLE SMALL SUB-DIAGONAL
C                                  ELEMENT
   20 NPL = IEN+K
      DO 25 LL=K,IEN
         LB = NPL-LL
         IF (LB.EQ.K) GO TO 30
         S = DABS(H(LB-1,LB-1))+DABS(H(LB,LB))
         IF (S.EQ.0.0D0) S = RNORM
         IF (DABS(H(LB,LB-1)).LE.RDELP*S) GO TO 30
   25 CONTINUE
C
   30 X = H(IEN,IEN)
      IF (LB.EQ.IEN) GO TO 110
      Y = H(NA,NA)
      W = H(IEN,NA)*H(NA,IEN)
      IF (LB.EQ.NA) GO TO 115
      IF (ITS.EQ.30) GO TO 250
C                                  FORM SHIFT
      IF (ITS.NE.10 .AND. ITS.NE.20) GO TO 40
      T = T+X
      DO 35 I=K,IEN
         H(I,I) = H(I,I)-X
   35 CONTINUE
      S = DABS(H(IEN,NA))+DABS(H(NA,IENM2))
      X = P7*S
      Y = X
      W = -P4*S*S
   40 ITS = ITS+1
C                                  LOOK FOR TWO CONSECUTIVE SMALL
C                                  SUB-DIAGONAL ELEMENTS
      NAML = IENM2+LB
      DO 45 MM=LB,IENM2
         M = NAML-MM
         ZZ = H(M,M)
         R = X-ZZ
         S = Y-ZZ
         P = (R*S-W)/H(M+1,M)+H(M,M+1)
         Q = H(M+1,M+1)-ZZ-R-S
         R = H(M+2,M+1)
         S = DABS(P)+DABS(Q)+DABS(R)
         P = P/S
         Q = Q/S
         R = R/S
         IF (M.EQ.LB) GO TO 50
         IF (DABS(H(M,M-1))*(DABS(Q)+DABS(R)).LE.RDELP*DABS(P)*(DABS(H(M
     *-1,
     *   M-1))+DABS(ZZ)+DABS(H(M+1,M+1)))) GO TO 50
   45 CONTINUE
   50 MP2 = M+2
      DO 55 I=MP2,IEN
         H(I,I-2) = ZERO
         IF (I.EQ.MP2) GO TO 55
         H(I,I-3) = ZERO
   55 CONTINUE
C                                  DOUBLE QR STEP INVOLVING ROWS
C                                  L TO EN AND COLUMNS M TO EN
      DO 105 KA=M,NA
         NOTLAS = KA.NE.NA
         IF (KA.EQ.M) GO TO 60
         P = H(KA,KA-1)
         Q = H(KA+1,KA-1)
         R = ZERO
         IF (NOTLAS) R = H(KA+2,KA-1)
         X = DABS(P)+DABS(Q)+DABS(R)
         IF (X.EQ.ZERO) GO TO 105
         P = P/X
         Q = Q/X
         R = R/X
   60    CONTINUE
         S = DSIGN(DSQRT(P*P+Q*Q+R*R),P)
         IF (KA.EQ.M) GO TO 65
         H(KA,KA-1) = -S*X
         GO TO 70
   65    IF (LB.NE.M) H(KA,KA-1) = -H(KA,KA-1)
   70    P = P+S
         X = P/S
         Y = Q/S
         ZZ = R/S
         Q = Q/P
         R = R/P
C                                  ROW MODIFICATION
         DO 80 J=KA,N
            P = H(KA,J)+Q*H(KA+1,J)
            IF (.NOT.NOTLAS) GO TO 75
            P = P+R*H(KA+2,J)
            H(KA+2,J) = H(KA+2,J)-P*ZZ
   75       H(KA+1,J) = H(KA+1,J)-P*Y
            H(KA,J) = H(KA,J)-P*X
   80    CONTINUE
         J = MIN0(IEN,KA+3)
C                                  COLUMN MODIFICATION
         DO 90 I=1,J
            P = X*H(I,KA)+Y*H(I,KA+1)
            IF (.NOT.NOTLAS) GO TO 85
            P = P+ZZ*H(I,KA+2)
            H(I,KA+2) = H(I,KA+2)-P*R
   85       H(I,KA+1) = H(I,KA+1)-P*Q
            H(I,KA) = H(I,KA)-P
   90    CONTINUE
         IF (IZ.LT.N) GO TO 105
C                                  ACCUMULATE TRANSFORMATIONS
         DO 100 I=K,L
            P = X*Z(I,KA)+Y*Z(I,KA+1)
            IF (.NOT.NOTLAS) GO TO 95
            P = P+ZZ*Z(I,KA+2)
            Z(I,KA+2) = Z(I,KA+2)-P*R
   95       Z(I,KA+1) = Z(I,KA+1)-P*Q
            Z(I,KA) = Z(I,KA)-P
  100    CONTINUE
  105 CONTINUE
      GO TO 20
C                                  ONE ROOT FOUND
  110 H(IEN,IEN) = X+T
      WR(IEN) = H(IEN,IEN)
      WI(IEN) = ZERO
      IEN = NA
      GO TO 15
C                                  TWO ROOTS FOUND
  115 P = (Y-X)*P5
      Q = P*P+W
      ZZ = DSQRT(DABS(Q))
      H(IEN,IEN) = X+T
      X = H(IEN,IEN)
      H(NA,NA) = Y+T
      IF (Q.LT.ZERO) GO TO 135
C                                  REAL PAIR
      ZZ = P+DSIGN(ZZ,P)
      WR(NA) = X+ZZ
      WR(IEN) = WR(NA)
      IF (ZZ.NE.ZERO) WR(IEN) = X-W/ZZ
      WI(NA) = ZERO
      WI(IEN) = ZERO
      X = H(IEN,NA)
C                                  EMPLOY SCALE FACTOR IN CASE X AND
C                                  ZZ ARE VERY SMALL
      SCALE = DABS(X) + DABS(ZZ)
      R = SCALE * DSQRT( (X/SCALE)**2 + (ZZ/SCALE)**2 )
      P = X/R
      Q = ZZ/R
C                                  ROW MODIFICATION
      DO 120 J=NA,N
         ZZ = H(NA,J)
         H(NA,J) = Q*ZZ+P*H(IEN,J)
         H(IEN,J) = Q*H(IEN,J)-P*ZZ
  120 CONTINUE
C                                  COLUMN MODIFICATION
      DO 125 I=1,IEN
         ZZ = H(I,NA)
         H(I,NA) = Q*ZZ+P*H(I,IEN)
         H(I,IEN) = Q*H(I,IEN)-P*ZZ
  125 CONTINUE
      IF (IZ.LT.N) GO TO 140
C                                  ACCUMULATE TRANSFORMATIONS
      DO 130 I=K,L
         ZZ = Z(I,NA)
         Z(I,NA) = Q*ZZ+P*Z(I,IEN)
         Z(I,IEN) = Q*Z(I,IEN)-P*ZZ
  130 CONTINUE
      GO TO 140
C                                  COMPLEX PAIR
  135 WR(NA) = X+P
      WR(IEN) = X+P
      WI(NA) = ZZ
      WI(IEN) = -ZZ
  140 IEN = IENM2
      GO TO 15
C                                  ALL ROOTS FOUND, NOW
C                                  BACKSUBSTITUTE
  145 IF (IZ.LT.N) GO TO 9005
      IF (RNORM.EQ.ZERO) GO TO 9005
      DO 220 NN=1,N
         IEN = N+1-NN
         P = WR(IEN)
         Q = WI(IEN)
         NA = IEN-1
         IF (Q.GT.ZERO) GO TO 220
         IF (Q.LT.ZERO) GO TO 180
C                                  REAL VECTOR
         M = IEN
         H(IEN,IEN) = ONE
         IF (NA.EQ.0) GO TO 220
         DO 175 II=1,NA
            I = IEN-II
            W = H(I,I)-P
            R = H(I,IEN)
            IF (M.GT.NA) GO TO 155
            DO 150 J=M,NA
               R = R+H(I,J)*H(J,IEN)
  150       CONTINUE
  155       IF (WI(I).GE.ZERO) GO TO 160
            ZZ = W
            S = R
            GO TO 175
  160       M = I
            IF (WI(I).NE.ZERO) GO TO 165
            T = W
            IF (W.EQ.ZERO) T = RDELP*RNORM
            H(I,IEN) = -R/T
            GO TO 175
C                                  SOLVE REAL EQUATIONS
  165       X = H(I,I+1)
            Y = H(I+1,I)
            Q = (WR(I)-P)*(WR(I)-P)+WI(I)*WI(I)
            T = (X*S-ZZ*R)/Q
            H(I,IEN) = T
            IF (DABS(X).LE.DABS(ZZ)) GO TO 170
            H(I+1,IEN) = (-R-W*T)/X
            GO TO 175
  170       H(I+1,IEN) = (-S-Y*T)/ZZ
  175    CONTINUE
C                                  END REAL VECTOR
         GO TO 220
C                                  LAST VECTOR COMPONENT CHOSEN
C                                    IMAGINARY SO THAT EIGENVECTOR
C                                    MATRIX IS TRIANGULAR
  180    M = NA
C                                  COMPLEX VECTOR
         IF (DABS(H(IEN,NA)).LE.DABS(H(NA,IEN))) GO TO 185
         H(NA,NA) = Q/H(IEN,NA)
         H(NA,IEN) = -(H(IEN,IEN)-P)/H(IEN,NA)
         GO TO 190
  185    CONTINUE
         Z3 = DCMPLX(ZERO,-H(NA,IEN))/DCMPLX(H(NA,NA)-P,Q)
         H(NA,NA) = T3(1)
         H(NA,IEN) = T3(2)
  190    H(IEN,NA) = ZERO
         H(IEN,IEN) = ONE
         IENM2 = NA-1
         IF (IENM2.EQ.0) GO TO 220
         DO 215 II=1,IENM2
            I = NA-II
            W = H(I,I)-P
            RA = ZERO
            SA = H(I,IEN)
            DO 195 J=M,NA
               RA = RA+H(I,J)*H(J,NA)
               SA = SA+H(I,J)*H(J,IEN)
  195       CONTINUE
            IF (WI(I).GE.ZERO) GO TO 200
            ZZ = W
            R = RA
            S = SA
            GO TO 215
  200       M = I
            IF (WI(I).NE.ZERO) GO TO 205
            Z3 = DCMPLX(-RA,-SA)/DCMPLX(W,Q)
            H(I,NA) = T3(1)
            H(I,IEN) = T3(2)
            GO TO 215
C                                  SOLVE COMPLEX EQUATIONS
  205       X = H(I,I+1)
            Y = H(I+1,I)
            VR = (WR(I)-P)*(WR(I)-P)+WI(I)*WI(I)-Q*Q
            VI = (WR(I)-P)*Q
            VI = VI+VI
            IF (VR.EQ.ZERO .AND. VI.EQ.ZERO) VR = RDELP*RNORM*(DABS(W)
     *      +DABS(Q)+DABS(X)+DABS(Y)+DABS(ZZ))
            Z3 = DCMPLX(X*R-ZZ*RA+Q*SA,X*S-ZZ*SA-Q*RA)/DCMPLX(VR,VI)
            H(I,NA) = T3(1)
            H(I,IEN) = T3(2)
            IF (DABS(X).LE.DABS(ZZ)+DABS(Q)) GO TO 210
            H(I+1,NA) = (-RA-W*H(I,NA)+Q*H(I,IEN))/X
            H(I+1,IEN) = (-SA-W*H(I,IEN)-Q*H(I,NA))/X
            GO TO 215
  210       CONTINUE
            Z3 = DCMPLX(-R-Y*H(I,NA),-S-Y*H(I,IEN))/DCMPLX(ZZ,Q)
            H(I+1,NA) = T3(1)
            H(I+1,IEN) = T3(2)
  215    CONTINUE
C                                  END COMPLEX VECTOR
  220 CONTINUE
C                                  END BACKSUBSTITUTION
C                                  VECTORS OF ISOLATED ROOTS
      DO 230 I=1,N
         IF (I.GE.K .AND. I.LE.L) GO TO 230
         DO 225 J=I,N
            Z(I,J) = H(I,J)
  225    CONTINUE
  230 CONTINUE
      IF (L.EQ.0) GO TO 9005
C                                  MULTIPLY BY TRANSFORMATION MATRIX
      DO 245 JJ=K,N
         J = N+K-JJ
         M = MIN0(J,L)
         DO 240 I=K,L
            ZZ = ZERO
            DO 235 KA=K,M
               ZZ = ZZ+Z(I,KA)*H(KA,J)
  235       CONTINUE
            Z(I,J) = ZZ
  240    CONTINUE
  245 CONTINUE
      GO TO 9005
C                                  NO CONVERGENCE AFTER 30 ITERATIONS
C                                  SET ERROR INDICATOR  TO THE INDEX
C                                  OF THE CURRENT EIGENVALUE
  250 IER = IEN
      DO 255 I=1,IEN
         WR(I) = ZERO
         WI(I) = ZERO
  255 CONTINUE
      IF (IZ.LT.N) GO TO 9000
      DO 265 I=1,N
         DO 260 J=1,N
            Z(I,J) = ZERO
  260    CONTINUE
  265 CONTINUE
 9000 CONTINUE
 9005 RETURN
      END
      SUBROUTINE EIGCH  (A,N,JOBN,D,Z,IZ,WK)
C
C-----------------------------------------------------------------------
C
C   PURPOSE             - EIGENVALUES AND (OPTIONALLY) EIGENVECTORS OF
C                           A COMPLEX HERMITIAN MATRIX
C
C   USAGE               - CALL EIGCH (A,N,JOBN,D,Z,IZ,WK,IER)
C
C   ARGUMENTS    A      - INPUT COMPLEX HERMITIAN MATRIX OF ORDER N
C                           WHOSE EIGENVALUES AND EIGENVECTORS ARE
C                           TO BE COMPUTED. INPUT A IS DESTROYED IF
C                           IJOB IS EQUAL TO 0 OR 1.
C                         NOTE - THE ROUTINE TREATS A AS A REAL VECTOR.
C                           AN EQUIVALENCE STATEMENT MAY BE REQUIRED-
C                           SEE DOCUMENT EXAMPLE.
C                N      - INPUT ORDER OF THE MATRIX A AND MATRIX Z.
C                JOBN   - INPUT OPTION PARAMETER. IF JOBN.GE.10
C                         A IS ASSUMED TO BE IN FULL COMPLEX STORAGE
C                         MODE (MUST BE DIMENSIONED EXACTLY N BY N).
C                         IF JOBN.LT.10 THEN A IS ASSUMED TO BE IN
C                         HERMITIAN STORAGE MODE.  DEFINE
C                         IJOB=MOD(JOBN,10).  THEN WHEN
C                           IJOB = 0, COMPUTE EIGENVALUES ONLY.
C                           IJOB = 1, COMPUTE EIGENVALUES AND EIGEN-
C                             VECTORS.
C                           IJOB = 2, COMPUTE EIGENVALUES, EIGENVECTORS
C                             AND PERFORMANCE INDEX.
C                           IJOB = 3, COMPUTE PERFORMANCE INDEX ONLY.
C                           IF THE PERFORMANCE INDEX IS COMPUTED, IT IS
C                           RETURNED IN WK(1). THE ROUTINES HAVE
C                           PERFORMED (WELL, SATISFACTORILY, POORLY) IF
C                           WK(1) IS (LESS THAN 1, BETWEEN 1 AND 100,
C                           GREATER THAN 100).
C                D      - OUTPUT VECTOR OF LENGTH N CONTAINING THE
C                           EIGENVALUES OF A.
C                Z      - OUTPUT N BY N COMPLEX MATRIX CONTAINING
C                           THE EIGENVECTORS OF A.
C                           THE EIGENVECTOR IN COLUMN J OF Z CORRES-
C                           PONDS TO THE EIGENVALUE D(J).
C                           IF IJOB = 0, Z IS NOT USED.
C                         NOTE - THE ROUTINE TREATS Z AS A REAL VECTOR
C                           OF LENGTH 2*N*N. AN APPROPRIATE EQUIVALENCE
C                           STATEMENT MAY BE REQUIRED. SEE DOCUMENT
C                           EXAMPLE.
C                IZ     - INPUT ROW DIMENSION OF MATRIX Z EXACTLY AS
C                           SPECIFIED IN THE DIMENSION STATEMENT IN THE
C                           CALLING PROGRAM. IZ MUST BE GREATER THAN
C                           OR EQUAL TO N IF IJOB IS NOT EQUAL TO ZERO.
C                WK     - WORK AREA, THE LENGTH OF WK DEPENDS
C                           ON THE VALUE OF IJOB, WHEN
C                           IJOB = 0, THE LENGTH OF WK IS AT LEAST 3N.
C                           IJOB = 1, THE LENGTH OF WK IS AT LEAST 3N.
C                           IJOB = 2, THE LENGTH OF WK IS AT LEAST
C                             N*N+4N.
C                           IJOB = 3, THE LENGTH OF WK IS AT LEAST 1.
C
C-----------------------------------------------------------------------
C
C                                  SPECIFICATIONS FOR ARGUMENTS
      INTEGER            N,JOBN,IZ,IER
      DOUBLE PRECISION   A(1),D(N),Z(1),WK(1)
C                                  SPECIFICATIONS FOR LOCAL VARIABLES
      INTEGER            JER,K,I,NE,NTAU,NA,NI,NI2,IM1,J,IIZ,NZ,IIZ1,
     1                   IJOB,JR,IR,IJ,JI,NP1,
     2                   JZ,JZI,L,M,II,IL,KK,LZ,MZ,LK,KZ
      DOUBLE PRECISION   ANORM,ASUM,PI,SUMZ,SUMR,SUMI,S,TEN,RDELP,
     1                   ZERO,ONE,THOUS,AN,SIGNA
      DATA               RDELP/0.222045D-15/
      DATA               ZERO,ONE/0.0D0,1.0D0/,TEN/10.0D0/,THOUS/1000.0D
     *0/
      IF (JOBN.LT.10) GO TO 15
C                                  CONVERT TO HERMETIAN STORAGE MODE
      JR = N + N - 2
      IJ = 2
      K = 2
      DO 10 J=1,N
         DO 5 I=1,J
            A(K-1) = A(IJ-1)
            A(K) = -A(IJ)
            K = K+2
            IJ = IJ + 2
    5    CONTINUE
         IJ = IJ + JR
         JR = JR - 2
   10 CONTINUE
   15 IJOB = MOD(JOBN,10)
      IF (IJOB.GE.0.AND.IJOB.LE.3) GO TO 20
C                                  WARNING ERROR - IJOB IS NOT IN THE
C                                    RANGE
      WRITE(6,*) ' EIGCH: IJOB OUT OF RANGE, IJOB=1'
      IJOB = 1
      GO TO 25
   20 IF (IJOB.EQ.0) GO TO 45
   25 IF (IZ.GE.N) GO TO 30
C                                  WARNING ERROR - IZ IS LESS THAN N
C                                    EIGENVECTORS CAN NOT BE COMPUTED,
C                                    IJOB SET TO ZERO
      WRITE(6,*) ' EIGCH: INCORRECT IZ. ONLY EIGENVALUES CALCULATED'
      IJOB = 0
   30 K = 2
      DO 40 I=1,N
         IF (A(K).EQ.ZERO) GO TO 35
         A(K) = ZERO
C                                  WARNING ERROR - SOME DIAGONAL
C                                    ELEMENT(S) NOT REAL
         WRITE(6,*) ' EIGCH: SOME DIAGONAL ELEMENTS ARE NOT REAL!'
   35    K = K+I+I+2
   40 CONTINUE
      IF (IJOB.EQ.3) GO TO 110
   45 NE = 1
      NTAU = NE+N
      NA = NTAU+N+N
      NI = (N*(N+1))/2
      NI2 = NI+NI
      IF (IJOB.NE.2) GO TO 55
C                                  SAVE INPUT A IF IJOB = 2
      K = NA
      DO 50 I=1,NI2
         WK(K) = A(I)
         K = K+1
   50 CONTINUE
C                                  SEPARATE A INTO REAL AND IMAGINARY
C                                    PARTS
   55 IF (NI.LT.2) GO TO 70
      IM1 = 1
      DO 65 I=2,NI
         K = IM1+I
         PI = A(K)
         DO 60 J=1,IM1
            A(K) = A(K-1)
            K = K-1
   60    CONTINUE
         A(I) = PI
         IM1 = I
   65 CONTINUE
C                                  REDUCE HERMITIAN MATRIX TO A REAL
C                                    SYMMETRIC TRIDIAGONAL MATRIX
   70 CALL EHOUSH (A(1),A(NI+1),N,D,WK(NE),WK(NTAU))
      IIZ = 1
      IF (IJOB.NE.0) IIZ = IZ+IZ
      IF (IIZ.EQ.1) GO TO 85
C                                  SET Z TO AN IDENTITY MATRIX
      NZ = (IZ+IZ)*N
      DO 75 I=1,NZ
         Z(I) = ZERO
   75 CONTINUE
      K = 1
      IIZ1 = IIZ+1
      DO 80 I=1,N
         Z(K) = ONE
         K = K+IIZ1
   80 CONTINUE
C                                  COMPUTE EIGENVALUES AND EIGENVECTORS
   85 CALL EQRT2S (D,WK(NE),N,Z(1),IIZ)
      IF (IJOB.EQ.0) GO TO 9000
C                                  BACKTRANSFORM THE EIGENVECTORS
      CALL EHBCKH (A(1),A(NI+1),WK(NTAU),N,Z(1),Z(IZ+1),IIZ)
C                                  CONVERT Z (EIGENVECTORS) TO COMPLEX
C                                    FORMAT Z(IZ,N)
      JZ = 0
      DO 100 J=1,N
         JZI = JZ+IZ
         DO 90 I=1,N
            K = JZI+I
            WK(I) = Z(K)
   90    CONTINUE
         K = JZ+N
         L = K+N-1
         M = N
         DO 95 I=1,N
            Z(L) = Z(K)
            Z(L+1) = WK(M)
            K = K-1
            L = L-2
            M = M-1
   95    CONTINUE
         JZ = JZ+IZ+IZ
  100 CONTINUE
C                                  Z IS NOW IN COMPLEX FORMAT Z(IZ,N).
      IF (IJOB.NE.2) GO TO 9000
C                                  MOVE ORIGINAL MATRIX BACK TO A
      K = NA
      DO 105 I=1,NI2
         A(I) = WK(K)
         K = K+1
  105 CONTINUE
      WK(1) = THOUS
C                                  COMPUTE 1-NORM OF A
  110 ANORM = ZERO
      II = 1
      DO 120 I=1,N
         ASUM = ZERO
         IL = II
         KK = 2
         DO 115 L=1,N
            ASUM = ASUM+CDABS(DCMPLX(A(IL),A(IL+1)))
            IF (L.GE.I) KK = L+L
            IL = IL+KK
  115    CONTINUE
         ANORM = DMAX1(ANORM,ASUM)
         II = II+I+I
  120 CONTINUE
      IF (ANORM.EQ.ZERO) ANORM = ONE
C                                  COMPUTE PERFORMANCE INDEX
      PI = ZERO
      DO 135 I=1,N
         II = 1
         S = ZERO
         SUMZ = ZERO
         LZ = (IZ+IZ)*(I-1)+1
         LZ = IZ*(I-1)*2+1
         MZ = LZ
         DO 130 L=1,N
            LK = II
            KK = 2
            SUMZ = SUMZ+CDABS(DCMPLX(Z(LZ),Z(LZ+1)))
            SUMR = -D(I)*Z(LZ)
            SUMI = -D(I)*Z(LZ+1)
            KZ = MZ
            DO 125 K=1,N
               SIGNA = ONE
               IF (K.GT.L) SIGNA = -ONE
               SUMR = SUMR+A(LK)*Z(KZ)-SIGNA*A(LK+1)*Z(KZ+1)
               SUMI = SUMI+A(LK)*Z(KZ+1)+SIGNA*A(LK+1)*Z(KZ)
               IF (K.GE.L) KK = K+K
               LK = LK+KK
               KZ = KZ+2
  125       CONTINUE
            S = S+CDABS(DCMPLX(SUMR,SUMI))
            LZ = LZ+2
            II = II+L+L
  130    CONTINUE
         IF (SUMZ.EQ.ZERO) GO TO 135
         PI = DMAX1(PI,S/SUMZ)
  135 CONTINUE
      AN = N
      PI = PI/(ANORM*TEN*AN*RDELP)
      WK(1) = PI
      IF (JOBN.LT.10) GO TO 9000
C                                  CONVERT BACK TO FULL COMPLEX MODE
      NP1 = N + 1
      IJ = (N-1) * NP1
      IJ = IJ + IJ + 2
      K = N * NP1
      DO 145 JR=1,N
         J = N+1-JR
         DO 140 IR=1,J
            A(IJ-1) = A(K-1)
            A(IJ) = -A(K)
            K = K-2
            IJ = IJ - 2
  140    CONTINUE
         IJ = IJ - JR - JR
  145 CONTINUE
      JR = N + N
      II = 2
      JI = 2
      DO 155 I=1,N
         IJ = II
         DO 150 J=1,I
            A(IJ-1) = A(JI-1)
            A(IJ) = -A(JI)
            JI = JI+2
            IJ = IJ+JR
  150    CONTINUE
         JI = JI + JR - I - I
         II = II + 2
  155 CONTINUE
 9000 CONTINUE
 9005 RETURN
      END
      SUBROUTINE EHBCKH (AR,AI,TAU,N,ZR,ZI,IZ)
C                                  SPECIFICATIONS FOR ARGUMENTS
      INTEGER            N,IZ
      DOUBLE PRECISION   AR(1),AI(1),TAU(2,1),ZR(IZ,1),ZI(IZ,1)
C                                  SPECIFICATIONS FOR LOCAL VARIABLES
      INTEGER            J,K,NR,L,NRM1,INX1,INX2,K1
      DOUBLE PRECISION   DELTA,ZERO,ALPHA1,ALPHA2
      DATA               ZERO/0.0D0/
C                                  TRANSFORM THE EIGENVECTORS OF THE
C                                    REAL SYMMETRIC TRIDIAGONAL MATRIX
C                                    TO THOSE OF THE HERMITIAN TRIDIA-
C                                    GONAL MATRIX
C                                  FIRST EXECUTABLE STATEMENT
      DO 5 J=1,N
         DO 5 K=1,N
            ZI(J,K)=-ZR(J,K)*TAU(2,J)
            ZR(J,K)=ZR(J,K)*TAU(1,J)
    5 CONTINUE
      IF (N .LE. 2) GO TO 30
C                                  RECOVER THE HOUSEHOLDER MATRICES IN
C                                    REVERSE ORDER
      DO 25 L=3,N
         NR=N-L+2
         NRM1=NR-1
         INX1=(NR*(NRM1))/2+NR
         INX2=INX1-1
         IF (AI(INX1) .EQ. ZERO) GO TO 25
         DELTA=AI(INX1)* DSQRT(AR(INX2)**2+AI(INX2)**2)
         DO 20 J=1,N
            ALPHA1=ZERO
            ALPHA2=ZERO
            DO 10 K=NR,N
               K1=(K*(K-1))/2+NRM1
               ALPHA1=ALPHA1+AR(K1)*ZR(K,J)+AI(K1)*ZI(K,J)
               ALPHA2=ALPHA2-AI(K1)*ZR(K,J)+AR(K1)*ZI(K,J)
   10       CONTINUE
            ALPHA1=ALPHA1/DELTA
            ALPHA2=ALPHA2/DELTA
            DO 15 K=NR,N
               K1=(K*(K-1))/2+NRM1
               ZR(K,J)=ZR(K,J)-AR(K1)*ALPHA1+AI(K1)*ALPHA2
               ZI(K,J)=ZI(K,J)-AR(K1)*ALPHA2-AI(K1)*ALPHA1
   15       CONTINUE
   20    CONTINUE
   25 CONTINUE
   30 RETURN
      END
      SUBROUTINE EHOUSH (AR,AI,N,D,E,TAU)
C                                  SPECIFICATIONS FOR ARGUMENTS
      INTEGER            N
      DOUBLE PRECISION   AR(1),AI(1),D(1),E(1),TAU(2,1)
C                                  SPECIFICATIONS FOR LOCAL VARIABLES
      INTEGER            NM1,NN,I,NR,NRM1,L,INDX,J,JJ,INX1,INX2,JP1,KK,
     *                   IX,IM1
      DOUBLE PRECISION   RHO,TOLER,ZERO,ONE,T1,T2,TESTBB,VR,ROOT,DELTA,
     *                   RATIO,RDELP,Q1,Q2,X1,X2,TT1,TT2,BB
      DATA               ZERO/0.0D0/,ONE/1.0D0/
      DATA               RDELP/0.222045D-15/
C                                  FIRST EXECUTABLE STATEMENT
      NM1=N-1
      TOLER=ZERO
      NN=(N*(N+1))/2
      DO 5 I=1,NN
         T1=DABS(AR(I))
         T2=DABS(AI(I))
         IF(T2.GT.T1) T1=T2
         IF (T1.GT.TOLER) TOLER=T1
    5 CONTINUE
      TESTBB=RDELP*TOLER
      IF (N.LE.2) GO TO 65
C                                  PERFORM N - 2 SIMILARITY
C                                    TRANSFORMATIONS
      DO 60 NR=2,NM1
         NRM1=NR-1
         VR=ZERO
         TAU(1,NR)=ZERO
         TAU(2,NR)=ZERO
         TAU(2,1)=ZERO
         DO 10 L=NR,N
            INDX=(L*(L-1))/2+NRM1
            VR=AR(INDX)**2+AI(INDX)**2+VR
   10    CONTINUE
         INDX=(NR*NRM1)/2+NRM1
         IF ((TESTBB)**2 .GE. VR) GO TO 60
         ROOT = CDABS(DCMPLX(AR(INDX),AI(INDX)))*DSQRT(VR)
         IF(ROOT.NE.ZERO) GO TO 15
         AR(INDX)=DSQRT(VR)
         DELTA=VR
         TAU(1,1)=-AR(INDX)
         GO TO 20
   15    DELTA=VR+ROOT
         RATIO=VR/ROOT
         TAU(1,1)=-RATIO*AR(INDX)
         TAU(2,1)= RATIO*AI(INDX)
         AR(INDX)=(RATIO+ONE)*AR(INDX)
         AI(INDX)=(RATIO+ONE)*AI(INDX)
C                                  THE MATRIX TO BE USED IN THE
C                                    SIMILARITY TRANSFORMATION HAS
C                                    BEEN DETERMINED. THE TRANSFOR-
C                                    MATION FOLLOWS
   20    DO 35 J=NR,N
            JJ=(J*(J-1))/2
            INDX=JJ+NRM1
            TAU(1,J)=AR(INDX)/DELTA
            TAU(2,J)=AI(INDX)/DELTA
            D(J)=ZERO
            E(J)=ZERO
            DO 25 L=NR,J
               INX1=(L*(L-1))/2+NRM1
               INX2=JJ+L
               D(J)= D(J)+AR(INX2)*AR(INX1)-AI(INX2)*AI(INX1)
               E(J)= E(J)+AR(INX2)*AI(INX1)+AI(INX2)*AR(INX1)
   25       CONTINUE
            JP1=J+1
            IF (JP1 .GT. N) GO TO 40
            DO 30 L=JP1,N
               KK=(L*(L-1))/2
               INX1=KK+NRM1
               INX2=KK+J
               D(J)=D(J)+AR(INX2)*AR(INX1)+AI(INX2)*AI(INX1)
               E(J)=E(J)+AR(INX2)*AI(INX1)-AI(INX2)*AR(INX1)
   30       CONTINUE
   35    CONTINUE
   40    RHO=ZERO
         DO 45 L=NR,N
            RHO=RHO+D(L)*TAU(1,L)+E(L)*TAU(2,L)
   45    CONTINUE
         IX=(NRM1*(NR-2))/2
         DO 55 I=NR,N
            IX=IX+I-1
            INX2=IX+NRM1
            DO 50 J=NR,I
               INX1=IX+J
               X1=TAU(1,I)*D(J)+TAU(2,I)*E(J)
               X2=TAU(2,I)*D(J)-TAU(1,I)*E(J)
               Q1=D(I)-RHO*AR(INX2)
               Q2=E(I)-RHO*AI(INX2)
               T1=Q1*TAU(1,J)+Q2*TAU(2,J)
               T2=Q2*TAU(1,J)-Q1*TAU(2,J)
               AR(INX1)=AR(INX1)-X1-T1
               AI(INX1)=AI(INX1)-X2-T2
   50       CONTINUE
   55    CONTINUE
         TAU(1,NR)=TAU(1,1)
         TAU(2,NR)=TAU(2,1)
   60 CONTINUE
C                                  THE MATRIX HAS BEEN REDUCED TO TRI-
C                                    DIAGONAL HERMITIAN FORM. THE SUB-
C                                    DIAGONAL HAS BEEN TEMPORARILY
C                                    STORED IN VECTOR TAU. STORE THE
C                                    DIAGONAL OF THE REDUCED MATRIX IN D
   65 INDX=0
      DO 70 I=1,N
         INDX=INDX+I
         D(I)=AR(INDX)
   70 CONTINUE
C                                  PERFORM THE DIAGONAL UNITARY SIMILA-
C                                    RITY TRANSFORMATION
      TAU(1,1)=ONE
      TAU(2,1)=ZERO
      E(1)=ZERO
      IF (N .EQ. 1) GO TO 85
      INDX=(N*NM1)/2+NM1
      TAU(1,N)=AR(INDX)
      TAU(2,N)=-AI(INDX)
C                                  CALCULATE SUBDIAGONAL E OF THE REAL
C                                    SYMMETRIC TRIDIAGONAL MATRIX. CAL-
C                                    CULATE TAU, THE DIAGONAL OF THE
C                                    DIAGONAL UNITARY MATRIX
      INDX=1
      DO 80 I=2,N
         INDX=INDX+I
         IM1=I-1
         BB= DSQRT(TAU(1,I)**2+TAU(2,I)**2)
         E(I)=BB
         AI(INDX)=BB
         IF (TESTBB .LT. BB) GO TO 75
         TAU(1,I)=ONE
         TAU(2,I)=ZERO
         BB=ONE
   75    TT1=TAU(1,I)*TAU(1,IM1)-TAU(2,I)*TAU(2,IM1)
         TT2=TAU(1,I)*TAU(2,IM1)+TAU(2,I)*TAU(1,IM1)
         TAU(1,I)=TT1/BB
         TAU(2,I)=TT2/BB
   80 CONTINUE
   85 RETURN
      END
