!----------------------------------!
subroutine setupfls(n)
  
  use numbers

  implicit none
  integer n, i
  
  nfls=n
  allocate(rfl(0:nfls))
  rfl(0)=1.d0
  do i=1, nfls
     rfl(i)=i*rfl(i-1)
  end do
  
end subroutine setupfls
!----------------------------------!
