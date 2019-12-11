C---------------------------------------------------------------------
      SUBROUTINE  JWKB   (X,ETA,XL, FJWKB,GJWKB, IEXP)
      DOUBLE PRECISION    X,ETA,XL, FJWKB,GJWKB, DZERO, PHI, PHI10
C----------------------------------------------------------------------
C---- COMPUTES JWKB APPROXIMATIONS TO COULOMB FUNCTIONS  FOR XL .GE. 0.
C---- AS MODIFIED BY BIEDENHARN ET AL. PHYS REV 97 (1955) 542-554
C---- CALCULATED IN SINGLE, RETURNED IN DOUBLE PRECISION VARIABLES
C---- CALLS MAX, SQRT, LOG, EXP, ATAN2, REAL, INT
C     AUTHOR:    A.R.BARNETT   FEB 1981    LAST UPDATE MARCH 1991
C----------------------------------------------------------------------
      REAL    ZERO,HALF,ONE,SIX,TEN,RL35,ALOGE
      REAL    GH2,XLL1,HLL,HL,SL,RL2,GH
      INTEGER IEXP, MAXEXP
      PARAMETER  ( MAXEXP = 300 )
      DATA  ZERO,HALF,ONE,SIX,TEN  /0.0E0, 0.5E0, 1.0E0, 6.0E0, 1.0E1/
      DATA DZERO,RL35,ALOGE /0.0D0, 35.0E0, 0.43429 45 E0 /
C----------------------------------------------------------------------
C---- CHOOSE MAXEXP NEAR MAX EXPONENT RANGE 
C---- E.G. 1.D300 FOR DOUBLE PRECISION
C----------------------------------------------------------------------
      GH2   = X * (ETA + ETA - X)
      XLL1  = MAX( XL * XL + XL, DZERO )
      IF( GH2 + XLL1 .LE. ZERO ) RETURN
      HLL = XLL1 + SIX / RL35
      HL = SQRT(HLL)
      SL = ETA / HL + HL / X
      RL2 = ONE + ETA * ETA / HLL
      GH = SQRT(GH2 + HLL) / X
      PHI = X*GH - HALF*( HL*LOG((GH + SL)**2 / RL2) - LOG(GH) )
      IF( ETA .NE. ZERO ) PHI = PHI - ETA * ATAN2(X*GH,X - ETA)
      PHI10 = -PHI * ALOGE
      IEXP = INT(PHI10)
      IF( IEXP .GT. MAXEXP ) THEN
       GJWKB = TEN**(PHI10 - REAL(IEXP))
      ELSE
       GJWKB = EXP(-PHI)
       IEXP = 0
      ENDIF
      FJWKB = HALF / (GH * GJWKB)
      RETURN
      END
