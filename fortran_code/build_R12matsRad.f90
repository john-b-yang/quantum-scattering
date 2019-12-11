!-----------------------------------!
subroutine build_R12matsRadau()

  use dvrecsRad
  implicit none

  integer lam, i, j
  complex*16 fval, faci, facj, b0

  if(mq.eq.0) then
     b0=az(nel)
  else
     b0=(2.d0*dble(mq)*eit)/alphaRad+R0 !!!effective surface term of Radau quad
  end if
  
  do lam = 0, lamax
     do i=1, nbas
        faci=1.d0/(xz(i)*sqrt(wz(i)))
        if ((i>(rbas+1)).and.(mq.ne.0)) faci=faci*exp(-alphaRad*conjg(eit)*(xz(i)-R0))
        do j=1, nbas
           facj=1.d0/(xz(j)*sqrt(wz(j)))
           if ((j>(rbas+1)).and.(mq.ne.0)) facj=facj*exp(-alphaRad*conjg(eit)*(xz(j)-R0))

           fval=dble(2*lam+1)*TIXX(j, i, lam)*faci*facj
           fval=fval+((xz(j)*xz(i))**(lam))/(b0**(2*lam+1))
           TIXX(j, i, lam) = fval
        end do
     end do
  end do
  
end subroutine build_R12matsRadau
!-----------------------------------!
