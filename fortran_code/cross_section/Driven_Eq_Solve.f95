SUBROUTINE Driven_Eq_Solve()
!
!  Solves linear equations for (ES -H) Psi_scattered = rhs
!
USE parameters, ONLY: unit_stdout,DBL
USE declare
!
IMPLICIT NONE
!
COMPLEX(kind=DBL), DIMENSION(nbas,nbas) :: ES_minus_H ! matrix for ES-H 
COMPLEX(kind=DBL), DIMENSION(nbas) :: Psi_scattered_scratch 
INTEGER :: lda, ldb, nrhs
INTEGER, ALLOCATABLE, DIMENSION(:) :: ipiv
!
IF (l_start>0) THEN !Ionizing from orbital other than S 
!**
!
 ALLOCATE(Psi_scattered(2*nbas,l_start-1:l_start+1))
 !load the different hamiltonians for each l
 DO l=l_start-1,l_start+1,2 
   !load rhs into solution vector for passing to linear equation routine
   DO n=1,nbas
    Psi_scattered_scratch(n) = r_eig_vec(n,l)
   ENDDO
   DO n=1,nbas
    DO m=1,nbas
       ES_minus_H(n,m) = E_tot*S_mat(n,m) - H_mat(n,m,l)
    ENDDO
   ENDDO
   ! linear equations
   lda = nbas
   ldb = nbas 
   nrhs = 1
   ALLOCATE(ipiv(nbas))
   CALL zgesv(nbas,nrhs,ES_minus_H,lda,ipiv,Psi_scattered_scratch,ldb,info)
   IF(info/=0) THEN
    WRITE(unit_stdout,*) " problem in zgesv solving lin. eqs. ", info
    STOP
   ENDIF
   DEALLOCATE(ipiv)
!
  DO n=1,nbas
    Psi_scattered(n,l) = Psi_scattered_scratch(n)
  ENDDO
!
!Checking****
  IF (E_tot==5.d0) THEN
   DO j=1,nbas
    WRITE(unit_Psi_scattered,'(3d25.15)')  REAL(x(j+1)), REAL(Psi_scattered(j,l))!/sqrt(wt(j+1)))  
   ENDDO
    WRITE(unit_Psi_scattered,*)  
    WRITE(unit_Psi_scattered,*)  
   ENDIF
!
 ENDDO
!**
ELSE !Ionizing from s orbital
!**
 ALLOCATE(Psi_scattered(nbas,1))
 l=l_start+1 !Ionize from s to p state
!
! load rhs into solution vector for passing to linear equation routine
   DO n = 1, nbas
    Psi_scattered_scratch(n) = r_eig_vec(n,l)
    !Psi_scattered_scratch(n) = rhs(n)
   ENDDO
!
   DO n=1,nbas
    DO m=1,nbas
         ES_minus_H(n,m) = E_tot*S_mat(n,m) - H_mat(n,m,l)
    ENDDO
   ENDDO
 ! linear equations
 lda = nbas
 ldb = nbas 
 nrhs = 1
 ALLOCATE(ipiv(nbas))
 CALL zgesv(nbas,nrhs,ES_minus_H,lda,ipiv,Psi_scattered_scratch,ldb,info)
 IF(info/=0) THEN
   WRITE(unit_stdout,*) " problem in zgesv solving lin. eqs. ", info
   STOP
 ENDIF
 DEALLOCATE(ipiv)
!
  DO n=1,nbas
    Psi_scattered(n,l) = Psi_scattered_scratch(n)
  ENDDO
!
!Checking******
 IF (E_tot==5.d0) THEN
  DO j=1,nbas
    WRITE(unit_Psi_scattered,'(3d25.15)')  REAL(x(j+1)), REAL(Psi_scattered(j,l))!/sqrt(wt(j+1))) 
  ENDDO
 ENDIF
!
!**
ENDIF
RETURN
END SUBROUTINE  
