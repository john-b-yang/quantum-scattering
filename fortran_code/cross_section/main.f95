PROGRAM Single_Photoionization 
!
!*******MODULES********
USE parameters
USE declare
!**********************
!
IMPLICIT NONE
!
CALL FEM_DVR() !Building grid
CALL Hamiltonian_build() !Building Hamiltonian
!CALL Right_Hand_Side() !Analytic wavefunctions for checking
CALL cdiag() !Get Eigenvectors and Eigenvalues
!****Bound Energies***** Energy = -Znuc**2/(2*n**2)
E_1 = eig(1) ! These depend on the l value, just like orbital  
E_2 = eig(2) ! These depend on the l value, just like orbital  
!print *, 'H_1s eigenvalue 1',E_1
!print *, 'H_1s eigenvalue 2',E_2
!************************
!Set up for Cross Section calculation
energy_min=0
energy_max=2
step=0.01
n=nint((energy_max-energy_min)/step)
DO i=0,n+1
    !E_tot = i*step
    E_tot = 5.d0 !Check 
    omega = E_tot - E_1
    CALL Driven_Eq_Solve()
    CALL Opt_Theorem()
    CALL Flux()
    DEALLOCATE(Psi_scattered)
ENDDO
!
END PROGRAM
