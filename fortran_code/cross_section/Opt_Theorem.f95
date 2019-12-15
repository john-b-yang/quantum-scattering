SUBROUTINE Opt_Theorem()
!
! Calculates Cross Section from Optical Theorem: (4*Pi*omega)*Im(<Psi_scattered|\mu|Psi_bound>)/c
USE parameters
USE declare
!
IMPLICIT NONE
!
COMPLEX(kind=DBL) :: Full_optmat
REAL(kind=DBL) :: Ang_integral, k_momentum, dipole_mat_imaginary_part
COMPLEX(kind=DBL), ALLOCATABLE, DIMENSION(:) :: right_vec,left_vec,Opt_ther_mat
!
IF (l_start>0) THEN !Ionizing from orbital other than s
!**
!
 Full_optmat = (0.d0,0.d0)
 ALLOCATE(Opt_ther_mat(l_start-1:l_start+1))
 DO l=l_start-1,l_start+1
  IF (l/=l_start) THEN
   ALLOCATE(right_vec(nbas))
   ALLOCATE(left_vec(nbas))
   DO n=1,nbas
    right_vec(n) = r_eig_vec(n,l)
    left_vec(n) = Psi_scattered(n,l)
   ENDDO
   Opt_ther_mat(l) = dot_product(left_vec,right_vec)
   DEALLOCATE(right_vec)
   DEALLOCATE(left_vec)
  ENDIF
  Full_optmat = Full_optmat + Opt_ther_mat(l)
 ENDDO
 DEALLOCATE(Opt_ther_mat)
!**
ELSE !Ionizing from s orbital
!**
 l=1
  ALLOCATE(Opt_ther_mat(l))
  ALLOCATE(right_vec(nbas))
  ALLOCATE(left_vec(nbas))
  DO n=1,nbas
   right_vec(n) = r_eig_vec(n,l)
   !right_vec(n) = rhs(n)
   left_vec(n) = Psi_scattered(n,l)
  ENDDO
!
  Opt_ther_mat(l) = dot_product(left_vec,right_vec)
  Full_optmat = Opt_ther_mat(l)
  DEALLOCATE(right_vec)
  DEALLOCATE(left_vec)
  DEALLOCATE(Opt_ther_mat)
!
!**
ENDIF
!
  dipole_mat_imaginary_part = AIMAG(Full_optmat)
  Ang_integral = 1.d0/Sqrt(3.d0)! Average orientation 
!
IF (l_start>0) THEN
  Cross_section = a0_sqr_to_mb*4.d0*pi*omega*alpha_fine_structure*dipole_mat_imaginary_part*Ang_integral**2
ELSE
  Cross_section = a0_sqr_to_mb*4.d0*pi*omega*alpha_fine_structure*dipole_mat_imaginary_part
ENDIF
!
WRITE(unit_Cross_SectionvsEnergy_Opt,*) E_tot, Cross_Section
!
!**************Possibly Useful Stuff*************************
!SQRT(4.d0*Pi/3.d0)*Cross_section = Cross_section*(1.d0/(2.d0))
!Cross_section = Cross_section*(1.d0/(k_momentum))
!k_momentum = Sqrt(2.d0*(omega)) !Dispersion relationship omega(k) = k**2*hbar/(2*m)
!
RETURN
END
