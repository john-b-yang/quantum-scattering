FUNCTION V_pot(z)
!
! Time-independent potential
!
IMPLICIT NONE
INTEGER, PARAMETER :: DBL = Selected_Real_Kind(p=13,r=200)
COMPLEX(kind=DBL) :: V_pot, z
REAL(kind=DBL) :: Znuc, lang
!
lang = 1.d0
Znuc = 1.d0
! Noro Taylor model problem with one resonance
  V_pot = 7.5d0*(z**2)*Exp(-z)
! Coulomb   with centrifugal potential  
!  V_pot = lang*(lang+1.d0)/(2.d0*z**2) - Znuc/z
! exponential  potential 
!  V_pot = 2/(2.d0*z**2) - 3.d0*Exp(-z) 
!  resonance (narrow)
!     V_pot  = -3.d0*zexp(-z) + 0.4*zexp(-0.5*(z-5.d0)**2)
!  resonance, one narrow one broad -- modified Coulomb potential 
!  V_pot = lang*(lang+1.d0)/(2.d0*z**2) - Znuc/z + 0.1d0*zexp(-0.05*(z-20.d0)**2)
!  V_pot = lang*(lang+1.d0)/(2.d0*z**2) - Znuc/z + 0.05d0*zexp(-0.025*(z-20.d0)**2)
RETURN
END

