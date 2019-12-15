!-----------------------------------!
subroutine build_KEmats()
  
  use dvrecs
  use numbers, only: pi, eye
  implicit none

  integer i, j, l, ii, jj
  integer i1, i2, o, m, lam, info
  complex*16 de

  integer ipiv(nbas)
  real*8 dq_l(nq, nq), dq_r(mq, mq), lag_alpha
  complex*16 D2(nbas, nbas), scal_w
  complex*16 zImat(nbas, nbas), zR1(nbas, nbas)

  lag_alpha=0.3d0 !Scrinzi
  eit2=exp(eye*(pi*theta*(-2)/180.d0))
  eit_conj=conjg(eit)

  call build_dq_nw(nq, dq_l, 1)
  call build_dq_nw(mq, dq_r, 2)
  
  D2=0.d0
    do l=1, nel
        i1=1
        i2=nq
        o=(nq-1)*(l-1)
        if(l==1) i1=2
        if(rel==nel)then
            if(l==nel) i2 = nq-1 
        else
            if(l==nel) i2 = mq-1 
        end if
        do i=i1, i2
            ii=(i-1)+o
            do j=i1, i2
                jj=(j-1)+o
                de=0.d0
                if (l<=rel) then !Labatto body
                    scal_w=2.d0/(az(l)-az(l-1))/sqrt(wz(ii)*wz(jj))
                    do m=1, nq
                        de=de+dq_l(m, i)*dq_l(m, j)*scal_w*wq_l(m)
                    end do
                    D2(ii, jj)=D2(ii, jj)+de
                else !Radau tail
                    scal_w=2*lag_alpha/(eit*sqrt(wz(ii)*wz(jj))) 
                    do m=1, mq
                        de=de+dq_r(m, i)*dq_r(m, j)*wq_r(m)*scal_w
                    end do
                    if (i==j) then
                        D2(ii,jj)=D2(ii,jj)+de- &
                        lag_alpha*eit_conj*scal_w*(wz(jj)*dq_r(j,i)+wz(ii)*dq_r(i,j))+ & 
                        lag_alpha**2*eit2
                    else
                        D2(ii,jj)=D2(ii,jj)+de- &
                        lag_alpha*eit_conj*scal_w*(wz(jj)*dq_r(j,i)+wz(ii)*dq_r(i,j))
                    end if
                end if
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
