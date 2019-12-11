!--------------------------------------!
subroutine cdiag(n, AA, eig)
  
  implicit none

  integer n, i, j
  complex*16 AA(n, n), eig(n)
  
  integer lwork, info
  real*8 rwork(2*n)
  complex*16 vl(1, n), evec(n, n), work(4*n), fval, alph(n)


  lwork = 4*n
  call zgeev('n', 'v', n, AA, n, eig, vl, 1, evec, n, work, lwork, rwork, info)


  do i=1, n
     do j=i+1, n
        if(dble(eig(j))<=dble(eig(i)))then
           fval=eig(i)
           eig(i)=eig(j)
           eig(j)=fval
           alph=evec(:, i)
           evec(:, i)=evec(:, j)
           evec(:, j)=alph           
        endif
     enddo
  enddo

  AA = evec

end subroutine cdiag
!--------------------------------------!
