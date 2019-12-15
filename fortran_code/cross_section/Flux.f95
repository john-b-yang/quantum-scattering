SUBROUTINE Flux()
!
! Outgoing amplitude from simple flux formula
!
USE parameters
USE declare
!
IMPLICIT NONE
!
INTEGER :: ipt, i_R0, l_ang, g 
REAL(kind=DBL) :: Ang_integral, k_momentum, Fluxer
COMPLEX(kind=DBL) ::  psival, psider, F, zpt, step_back_from_R0
COMPLEX(kind=DBL) :: Psi_evaluate, Psi_derivative_evaluate
COMPLEX(kind=DBL), DIMENSION(nbas) :: Psi_scattered_scratch !temp  
!
! find R_0 on grid
i_R0 = 1
DO ipt=1,nbas
IF(IMAG(x(ipt))<=1.d-4) THEN
  i_R0 = ipt
ENDIF
ENDDO
!
! evaluate function and its derivative at a point just inside R0
step_back_from_R0 = 1.d0
zpt = x(i_R0) - step_back_from_R0
Fluxer = (0.d0,0.d0)
IF (l_start==0) THEN
    DO g=1,nbas
        Psi_scattered_scratch(g)=Psi_scattered(g,1)
    ENDDO
    psival =  Psi_evaluate(zpt,Psi_scattered_scratch)
    psider =  Psi_derivative_evaluate(zpt,Psi_scattered_scratch)
    F = (CONJG(psival)*psider - psival*CONJG(psider))/(0.d0,2.d0) 
    Fluxer = REAL(F)
ELSE
    DO l_ang=l_start-1,l_start+1,2
        DO g=1,nbas
         Psi_scattered_scratch(g)=Psi_scattered(g,l_ang)
        ENDDO
        psival =  Psi_evaluate(zpt,Psi_scattered_scratch)
        psider =  Psi_derivative_evaluate(zpt,Psi_scattered_scratch)
        F = (CONJG(psival)*psider - psival*CONJG(psider))/(0.d0,2.d0) 
        Fluxer = Fluxer + REAL(F)
    ENDDO
ENDIF
!
! General flux formula
!F = (CONJG(psival)*psider - psival*CONJG(psider))/(0.d0,2.d0) 
!WRITE(UNIT_stdout,*) "flux1 = ", F !Check
!
Ang_integral = 1.d0/Sqrt(3.d0)
IF (l_start/=0) THEN
    Cross_section = a0_sqr_to_mb*2.d0*Pi*alpha_fine_structure*omega*Fluxer*Ang_integral**2
ELSE
    Cross_section = a0_sqr_to_mb*2.d0*Pi*alpha_fine_structure*omega*Fluxer
ENDIF
WRITE(unit_Cross_SectionvsEnergy_Flux,*) E_tot, Cross_Section

!************Possible Useful Stuff***********************
!k_momentum = Sqrt(2.d0*omega)
!Cross_section = Cross_section*(1.d0/(2.d0))
!Cross_section = Cross_section*(2.d0/(Pi*k_momentum))

RETURN
END
