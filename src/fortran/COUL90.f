C----------------------------------------------------------------------
      SUBROUTINE COUL90(X, ETA, XLMIN,LRANGE, FC,GC,FCP,GCP, KFN,IFAIL)
C----------------------------------------------------------------------
C
C  COULOMB & BESSEL FUNCTION PROGRAM-- COUL90 -- USING STEED'S METHOD
C
C  COUL90 RETURNS ARRAYS FC = F, GC = G, FCP = (D/DX) F, GCP = (D/DX) G
C   FOR REAL X .GT. 0. ,REAL ETA (INCLUDING 0.), AND REAL XLMIN .GT.-1.
C   FOR (LRANGE+1) INTEGER-SPACED LAMBDA VALUES.
C   IT HENCE GIVES POSITIVE-ENERGY SOLUTIONS TO THE COULOMB SCHRODINGER
C   EQUATION, TO THE KLEIN-GORDON EQUATION AND TO SUITABLE FORMS OF
C   THE DIRAC EQUATION.    BY SETTING ETA = 0.0 AND RENORMALISING
C   SPHERICAL & CYLINDRICAL BESSEL FUNCTIONS ARE COMPUTED INSTEAD.
C----------------------------------------------------------------------
C   CALLING VARIABLES; ALL REALS ARE DOUBLE PRECISION (REAL*8)
C
C   X      - REAL ARGUMENT FOR COULOMB FUNCTIONS > 0.0
C      [ X > SQRT(ACCUR) : ACCUR IS TARGET ACCURACY 1.0D-14 ]
C   ETA     - REAL SOMMERFELD PARAMETER, UNRESTRICTED > = < 0.0
C   XLMIN   - REAL MINIMUM LAMBDA-VALUE (L-VALUE OR ORDER),
C      GENERALLY IN RANGE 0.0 - 1.0 AND MOST USUALLY 0.0
C   LRANGE  - INTEGER NUMBER OF ADDITIONAL L-VALUES : RESULTS RETURNED
C      FOR L-VALUES XLMIN TO XLMIN + LRANGE INCLUSIVE
C   FC ,GC  - REAL VECTORS F,G OF REGULAR, IRREGULAR COULOMB FUNCTIONS
C   FCP,GCP - REAL VECTORS FOR THE X-DERIVATIVES OF  F,G
C      THESE VECTORS TO BE OF LENGTH AT LEAST MINL + LRANGE
C      STARTING ELEMENT MINL = MAX0( IDINT(XLMIN+ACCUR),0 )
C   KFN     - INTEGER CHOICE OF FUNCTIONS TO BE COMPUTED :
C      = 0      REAL COULOMB FUNCTIONS AND DERIVATIVES F & G
C      = 1    SPHERICAL BESSEL      "      "     "      J & Y
C      = 2  CYLINDRICAL BESSEL      "      "     "      J & Y
C
C   PRECISION:  RESULTS TO WITHIN 2-3 DECIMALS OF "MACHINE ACCURACY"
C   IN OSCILLATING REGION X .GE. [ETA + SQRT{ETA**2 + XLM*(XLM+1)}]
C   I.E. THE TARGET ACCURACY ACCUR SHOULD BE 100 * ACC8 WHERE ACC8 IS
C   THE SMALLEST NUMBER WITH 1.+ACC8.NE.1. FOR OUR WORKING PRECISION.
C   THIS RANGES BETWEEN 4E-15 AND 2D-17 ON CRAY, VAX, SUN, PC FORTRANS
C   SO CHOOSE A SENSIBLE  ACCUR = 1.0D-14
C   IF X IS SMALLER THAN [ ] ABOVE THE ACCURACY BECOMES STEADLY WORSE :
C   THE VARIABLE PACCQ IN COMMON /STEED/ HAS AN ESTIMATE OF ACCURACY.
C----------------------------------------------------------------------
C   ERROR RETURNS      THE USER SHOULD TEST IFAIL ON EXIT
C
C   IFAIL ON INPUT IS SET TO 0      LIMIT = 20000
C   IFAIL IN OUTPUT =  0 : CALCULATIONS SATISFACTORY
C      =  1 : CF1 DID NOT CONVERGE AFTER LIMIT ITERATIONS
C      =  2 : CF2 DID NOT CONVERGE AFTER LIMIT ITERATIONS
C      = -1 : X < 1D-7 = SQRT(ACCUR)
C      = -2 : INCONSISTENCY IN ORDER VALUES (L-VALUES)
C----------------------------------------------------------------------
C  MACHINE-DEPENDENT PARAMETERS:    ACCUR - SEE ABOVE
C      SMALL - OFFSET FOR RECURSION = APPROX SQRT(MIN REAL NO.)
C      IE 1D-30 FOR IBM REAL*8,    1D-150 FOR DOUBLE PRECISION
C----------------------------------------------------------------------
C  PROGRAMMING HISTORY AND BACKGROUND: CPC IS COMPUTER PHYSICS COMMUN.
C  ORIGINAL PROGRAM  RCWFN       IN    CPC  8 (1974) 377-395
C                 +  RCWFF       IN    CPC 11 (1976) 141-142
C  FULL DESCRIPTION OF ALGORITHM IN    CPC 21 (1981) 297-314
C  REVISED STANDARD  COULFG      IN    CPC 27 (1982) 147-166
C  BACKGROUND MATERIAL IN J. COMP. PHYSICS 46 (1982) 171-188
C  CURRENT PROGRAM   COUL90  (FORTRAN77) SUPERCEDES THESE EARLIER ONES
C  (WHICH WERE WRITTEN IN FORTRAN 4) AND ALSO BY INCORPORATING THE NEW
C  LENTZ-THOMPSON ALGORITHM FOR EVALUATING THE FIRST CONTINUED FRACTION
C  ..SEE ALSO APPENDIX TO J. COMP. PHYSICS 64 (1986) 490-509
C----------------------------------------------------------------------
C  AUTHOR: A. R. BARNETT      MANCHESTER  MARCH   1981
C                             AUCKLAND    MARCH   1991
C----------------------------------------------------------------------
      INTEGER          LRANGE, KFN, IFAIL
      DOUBLE PRECISION X, ETA, XLMIN
      DOUBLE PRECISION FC (0:*), GC (0:*), FCP(0:*), GCP(0:*)
C---- ARRAYS INDEXED FROM 0 INSIDE SUBROUTINE: STORAGE FROM MINL
      DOUBLE PRECISION ACCUR,ACCH,SMALL, ONE,ZERO,HALF,TWO,TEN2, RT2DPI
      DOUBLE PRECISION XINV,PK,CF1,C,D,PK1,ETAK,RK2,TK,DCF1,DEN,XLM,XLL
      DOUBLE PRECISION EL,XL,RL,SL, F,FCMAXL,FCMINL,GCMINL,OMEGA,WRONSK
      DOUBLE PRECISION WI, A,B, AR,AI,BR,BI,DR,DI,DP,DQ, ALPHA,BETA
      DOUBLE PRECISION E2MM1, FJWKB,GJWKB, P,Q,PACCQ, GAMMA,GAMMAI
      INTEGER      IEXP, NFP, NPQ, L, MINL,MAXL, LIMIT
      LOGICAL      ETANE0, XLTURN
      PARAMETER      ( LIMIT = 900000, SMALL = 1.0D-150 )
      COMMON  /STEED/  PACCQ,NFP,NPQ,IEXP,MINL   
      COMMON  /DESET/  CF1,P,Q,F,GAMMA,WRONSK     
C----------------------------------------------------------------------
C     COUL90 HAS CALLS TO: DSQRT,DABS,MAX0,IDINT,DSIGN,DBLE,MIN
C----------------------------------------------------------------------
      DATA ZERO,ONE,TWO,TEN2,HALF /0.0D0, 1.0D0, 2.0D0, 1.0D2, 0.5D0/
      DATA RT2DPI /0.79788 45608 02865  D0/
CQ    DATA RT2DPI /0.79788 45608 02865 35587 98921 19868 76373 Q0/
C---- THIS CONSTANT IS  SQRT(TWO / PI):
C---- USE Q0 FOR IBM REAL*16: D0 FOR REAL*8 AND DOUBLE PRECISION
C---- CHANGE ACCUR TO SUIT MACHINE AND PRECISION REQUIRED
      ACCUR = 1.0D-14
      IFAIL = 0
      IEXP  = 1
      NPQ   = 0
      GJWKB = ZERO
      PACCQ = ONE
      IF(KFN .NE. 0) ETA = ZERO
      ETANE0 = ETA .NE. ZERO
      ACCH = SQRT(ACCUR)
C---- TEST RANGE OF X, EXIT IF.LE.SQRT(ACCUR) OR IF NEGATIVE
      IF( X .LE. ACCH ) GOTO 70
      IF( KFN .EQ. 2 ) THEN
       XLM = XLMIN - HALF
      ELSE
       XLM = XLMIN
      ENDIF
      IF( XLM .LE. -ONE .OR. LRANGE .LT. 0 ) GOTO 80
      E2MM1 = XLM * XLM + XLM
      XLTURN = X * (X - TWO * ETA) .LT. E2MM1
      E2MM1 = E2MM1 + ETA * ETA
      XLL = XLM + DBLE(LRANGE)
C---- LRANGE IS NUMBER OF ADDITIONAL LAMBDA VALUES TO BE COMPUTED
C---- XLL    IS MAX LAMBDA VALUE [ OR 0.5 SMALLER FOR J,Y BESSELS ]
C---- DETERMINE STARTING ARRAY ELEMENT (MINL) FROM XLMIN
      MINL = MAX0( IDINT(XLMIN + ACCUR),0 ) 
      MAXL = MINL + LRANGE
C---- EVALUATE CF1  =  F   =  DF(L,ETA,X)/DX   /   F(L,ETA,X)
      XINV = ONE / X
C---- UNNORMALISED F(MAXL,ETA,X)
      DEN = ONE 
      PK = XLL + ONE
      CF1 = ETA / PK + PK * XINV
      IF( ABS(CF1) .LT. SMALL ) CF1 = SMALL
      RK2 = ONE
      D = ZERO
      C = CF1
C---- BEGIN CF1 LOOP ON PK = K STARTING AT LAMBDA + 1: LENTZ-THOMPSON
      DO 10 L = 1 , LIMIT 
       PK1 = PK + ONE
       IF( ETANE0 ) THEN
        ETAK = ETA / PK
        RK2 = ONE + ETAK * ETAK
        TK = (PK + PK1) * (XINV + ETAK / PK1)
       ELSE
        TK = (PK + PK1) * XINV
       ENDIF
C---- DIRECT RATIO OF B CONVERGENTS
       D = TK - RK2 * D 
C---- INVERSE RATIO OF A CONVERGENTS
       C = TK - RK2 / C 
       IF( ABS(C) .LT. SMALL ) C = SMALL
       IF( ABS(D) .LT. SMALL ) D = SMALL
       D = ONE / D
       DCF1= D * C
       CF1 = CF1 * DCF1
       IF( D .LT. ZERO ) DEN = -DEN
       PK = PK1
C----  PROPER EXIT
       IF( ABS(DCF1-ONE) .LT. ACCUR ) GOTO 20     
10    CONTINUE
C---- ERROR EXIT
      GOTO 90                                      
20    NFP = PK - XLL - 1 
      F = CF1 
C---- DOWNWARD RECURRENCE TO LAMBDA = XLM; ARRAYS GC, GCP STORE RL, SL
      IF( LRANGE .GT. 0 ) THEN
       FCMAXL = SMALL * DEN
       FCP(MAXL) = FCMAXL * CF1
       FC (MAXL) = FCMAXL
       XL = XLL
       RL = ONE
       DO 30 L = MAXL, MINL+1, -1
        IF( ETANE0 ) THEN
         EL = ETA / XL
         RL = SQRT( ONE + EL * EL )
         SL = XL * XINV + EL
         GC (L) = RL 
         GCP(L) = SL
        ELSE
         SL = XL * XINV
        ENDIF
        FC (L-1) = ( FC(L) * SL + FCP(L) ) / RL
        FCP(L-1) = FC(L-1) * SL - FC (L) * RL
C---- END VALUE IS XLM
        XL = XL - ONE 
30     CONTINUE
       IF( ABS(FC(MINL)) .LT. ACCUR*SMALL ) FC(MINL) = ACCUR * SMALL
       F = FCP(MINL) / FC(MINL) 
       DEN = FC (MINL)
      ENDIF
C---------------------------------------------------------------------
C---- NOW WE HAVE REACHED LAMBDA = XLMIN = XLM
C---- EVALUATE CF2 = P + I.Q  USING STEED'S ALGORITHM (NO ZEROS)
C---------------------------------------------------------------------
      IF( XLTURN ) CALL JWKB(X,ETA,MAX(XLM,ZERO),FJWKB,GJWKB,IEXP)
      IF( IEXP .GT. 1 .OR. GJWKB .GT. ONE / (ACCH*TEN2)) THEN
       OMEGA = FJWKB
       GAMMA = GJWKB * OMEGA
       P = F
       Q = ONE
       ELSE 
       XLTURN = .FALSE.
       PK = ZERO
       WI = ETA + ETA
       P = ZERO
       Q = ONE - ETA * XINV
       AR = -E2MM1
       AI = ETA
       BR = TWO * (X - ETA)
       BI = TWO
       DR = BR / (BR * BR + BI * BI)
       DI = -BI / (BR * BR + BI * BI)
       DP = -XINV * (AR * DI + AI * DR)
       DQ = XINV * (AR * DR - AI * DI)
       DO 40 L = 1, LIMIT
        P = P + DP
        Q = Q + DQ
        PK = PK + TWO
        AR = AR + PK
        AI = AI + WI
        BI = BI + TWO
        D = AR * DR - AI * DI + BR
        DI = AI * DR + AR * DI + BI
        C = ONE / (D * D + DI * DI)
        DR = C * D
        DI = -C * DI
        A = BR * DR - BI * DI - ONE
        B = BI * DR + BR * DI
        C = DP * A - DQ * B
        DQ = DP * B + DQ * A
        DP = C
        IF( ABS(DP)+ABS(DQ) .LT. (ABS(P)+ABS(Q)) * ACCUR ) GOTO 50
40     CONTINUE
C---- ERROR EXIT
       GOTO 100                                   
50     NPQ = PK / TWO 
       PACCQ = HALF * ACCUR / MIN( ABS(Q),ONE )
       IF( ABS(P) .GT. ABS(Q) ) PACCQ = PACCQ * ABS(P)
C---------------------------------------------------------------------
C    SOLVE FOR FCMINL = F AT LAMBDA = XLM AND NORMALISING FACTOR OMEGA
C---------------------------------------------------------------------
       GAMMA = (F - P) / Q
       GAMMAI = ONE / GAMMA
       IF( ABS(GAMMA) .LE. ONE ) THEN
        OMEGA = SQRT( ONE + GAMMA * GAMMA )
       ELSE
        OMEGA = SQRT( ONE + GAMMAI* GAMMAI) * ABS(GAMMA)
       ENDIF
       OMEGA = ONE / ( OMEGA * SQRT(Q) )
       WRONSK = OMEGA
      ENDIF
C---------------------------------------------------------------------
C---- RENORMALISE IF SPHERICAL OR CYLINDRICAL BESSEL FUNCTIONS
C---------------------------------------------------------------------
      IF(KFN .EQ. 1) THEN 
C---- SPHERICAL
       ALPHA = XINV
       BETA = XINV
      ELSE
     - IF(KFN .EQ. 2) THEN
C----   CYLINDRICAL 
        ALPHA = HALF * XINV
        BETA = SQRT( XINV ) * RT2DPI
        ELSE 
C----   KFN = 0, COULOMB FUNCTIONS
        ALPHA = ZERO
        BETA = ONE
      ENDIF
      FCMINL = SIGN( OMEGA,DEN ) * BETA
      IF( XLTURN ) THEN
       GCMINL = GJWKB * BETA
      ELSE
       GCMINL = FCMINL * GAMMA
      ENDIF
C---- BESSEL SIGN
      IF( KFN .NE. 0 ) GCMINL = - GCMINL 
      FC (MINL) = FCMINL
      GC (MINL) = GCMINL
      GCP(MINL) = GCMINL * (P - Q * GAMMAI - ALPHA)
      FCP(MINL) = FCMINL * (F - ALPHA)
      IF( LRANGE .EQ. 0 ) RETURN
C---------------------------------------------------------------------
C---- UPWARD RECURRENCE FROM GC(MINL),GCP(MINL) STORED VALUES ARE RL,SL
C---- RENORMALISE FC,FCP AT EACH LAMBDA AND CORRECT REGULAR DERIVATIVE
C---- XL   = XLM HERE  AND RL = ONE , EL = ZERO FOR BESSELS
C---------------------------------------------------------------------
      OMEGA = BETA * OMEGA / ABS(DEN)
      XL = XLM
      RL = ONE
      DO 60 L = MINL+1 , MAXL 
       XL = XL + ONE
       IF( ETANE0 ) THEN
        RL = GC (L)
        SL = GCP(L)
       ELSE
        SL = XL * XINV
       ENDIF
       GC (L) = ( (SL - ALPHA) * GC(L-1) - GCP(L-1) ) / RL
       GCP(L) = RL * GC(L-1) - (SL + ALPHA) * GC(L)
       FCP(L) = OMEGA * ( FCP(L) - ALPHA * FC(L) )
       FC (L) = OMEGA * FC (L)
60    CONTINUE
      RETURN
C---- ERROR MESSAGES
70    IFAIL = -1
      WRITE(6,1000) X,ACCH
1000  FORMAT(' FOR X = ',1P,D12.3,'     TRY SMALL-X  SOLUTIONS',                 
     *' OR X IS NEGATIVE'/ ,' SQUARE ROOT (ACCURACY) =  ',D12.3/)               
      RETURN
80    IFAIL = -2
      WRITE (6,1001) LRANGE,XLMIN,XLM
1001  FORMAT(/' PROBLEM WITH INPUT ORDER VALUES: LRANGE, XLMIN, XLM = ',        
     *I10,1P,2D15.6/)                                                            
      RETURN
90    IFAIL = 1
      WRITE (99,1002) LIMIT, CF1,DCF1, PK,ACCUR
1002  FORMAT(' CF1 HAS FAILED TO CONVERGE AFTER ',I10,' ITERATIONS',/           
     *' CF1,DCF1,PK,ACCUR =  ',1P,4D12.3/)                                       
      RETURN
100   IFAIL = 2
      WRITE (99,1003) LIMIT,P,Q,DP,DQ,ACCUR
1003  FORMAT(' CF2 HAS FAILED TO CONVERGE AFTER ',I7,' ITERATIONS',/            
     *' P,Q,DP,DQ,ACCUR =  ',1P,4D17.7,D12.3/)                                   
      RETURN
      END
