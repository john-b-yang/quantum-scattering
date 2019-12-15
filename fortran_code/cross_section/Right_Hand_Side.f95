SUBROUTINE Right_Hand_Side()
!
! right hand side of driven Schroedinger Equation 
!
USE parameters 
USE declare
!
IMPLICIT NONE
!
REAL(kind=DBL) :: gaunter
INTEGER :: l1, l2
INTEGER :: m1, m2, m3
!
ALLOCATE(rhs(nbas))
l2 = 1 !Photon is a boson with integer=1 spin
    l1=1 
    m1=0
    m2=0
    m3=0
    CALL GauntCoef(l1,l2,m1,m2,m3,gaunter) !Get appropriate gaunt coefficient
rhs = 0.d0
DO i=1,nbas
    rhs(i) = gaunter*2.d0*SQRT(REAL(Znuc))**3*x(i+1)**2*exp(-Znuc*x(i+1))*sqrt(wt(i+1))
            WRITE(unit_wave,'(3d25.15)')  REAL(x(i+1)), REAL(rhs(i)/sqrt(wt(i+1)))
ENDDO
RETURN
END
