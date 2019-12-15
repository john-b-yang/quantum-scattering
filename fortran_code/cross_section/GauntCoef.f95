SUBROUTINE GauntCoef(l1,l2,m1,m2,m3,gaunter)
!
USE parameters
USE declare

  REAL(kind=DBL) :: tj,tjm1,tjm2,norm
  INTEGER :: const,new_m1
  INTEGER, INTENT(IN) :: l1, l2, m1, m2, m3
  REAL(kind=DBL), INTENT(OUT) :: gaunter
!
 !DO m1=-l1,l1
  !DO m2=-l2,l2
   !DO m3=-l3,l3
   !IF (m1+m2+m3==0) THEN !Spin multiplicity cannot change selection rule
!
        !Complex Conjugation Relation CONJ(Y_l^m) = (-1)^m*Y_l^(-m)
        const = (-1)**m1
        new_m1 = -m1
!
        CALL setupfls(90) !setup factorial rfl used in the subroutines below
        CALL threej(l1, l2, l_start, new_m1, m2, m3, tj)
        tjm1 = tj
        CALL threej0(l1, l2, l_start, tj)
        tjm2 = tj
        DEALLOCATE(rfl)
!
        norm = SQRT((2*l1+1)*(2*l2+1)*(2*l_start+1)/(4*pi)) !Connects to spherical harmonics
        IF (l1 ==2) THEN
            gaunter = SQRT(2.5)*SQRT(4.d0*pi/3.d0)*const*norm*tjm1*tjm2 
        ELSE
            gaunter = SQRT(4.d0*pi/3.d0)*const*norm*tjm1*tjm2 
        ENDIF
    !ENDDO
   !ENDDO
  !ENDDO
END SUBROUTINE GauntCoef
