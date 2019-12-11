!-----------------------------------!
subroutine build_KEmats()
  
  use dvrecs
  implicit none

  integer i, j, l, ii, jj
  integer i1, i2, o, m, lam, info
  complex*16 de

  integer ipiv(nbas)
  real*8 Dq(nq, nq)
  complex*16 D2(nbas, nbas)
  complex*16 zImat(nbas, nbas), zR1(nbas, nbas)


  call build_dq_nw(nq, Dq)
  D2=0.d0
  do l=1, nel
     i1=1
     i2=nq
     if(l==1) i1=2
     if(l==nel) i2=nq-1
     o=(nq-1)*(l-1)
     do i=i1, i2
        ii=(i-1)+o
        do j=i1, i2
           jj=(j-1)+o
           de=0.d0
           do m=1, nq
              de=de-dq(m, i)*dq(m, j)*2.d0/(az(l)-az(l-1))*(wq(m)/sqrt(wz(ii)*wz(jj)))
           end do
           D2(ii, jj)=D2(ii, jj)-de
        end do
     end do
  end do
  allocate(TXX(nbas, nbas, 0:lamax))
  TXX=0.d0
  do lam=0, lamax
     TXX(:, :, lam)=D2
     do i=1, nbas
        TXX(i, i, lam)=TXX(i, i, lam)+dble(lam*(lam+1))/(xz(i)**2)
     end do
  end do
  zImat=0.d0
  do i=1, nbas
     zImat(i, i)=1.d0
  end do
  allocate(TIXX(nbas, nbas, 0:lamax))
  do lam=0, lamax
     D2=TXX(:, :, lam)
     zR1=zImat
     call zgesv(nbas, nbas, D2, nbas, ipiv, zR1, nbas, info)
     TIXX(:, :, lam)=zR1
  end do
  
end subroutine Build_KEmats
!-----------------------------------!
