!----------------------------------------------!
subroutine Gcoeff(j1, m1, j2, m2, j3, m3, Cval)

  use numbers
  implicit none
  integer j1, j2, j3, m1, m2, m3
  real*8 Cval, tjay
  
  ! Aug.24,2015 - I think m1 is being conjugated in this formula
  ! by making m2 and m3 negative, but it might differ in sign change
  ! from by Gaunt function below 
  
  Cval=dble((-1)**(m3+m2))*sqrt(0.25d0*dble((2*j1+1)*(2*j2+1)*(2*j3+1))/pi)
  call threej0(j3, j2, j1, tjay)
  Cval=Cval*tjay
  call threej(j3, j2, j1, -m3, -m2, m1, tjay)
  Cval=Cval*tjay


end subroutine Gcoeff
!----------------------------------------------!
real*8 function Gaunt(j1, m1, j2, m2, j3, m3)

  use numbers
  implicit none
  integer j1, j2, j3, m1, m2, m3
  real*8 Cval, tjay

  ! does the conjugation of j1,m1 Ylm (integration of Y1*,dipole,Y1p)
  ! also written as a function, not a subroutine
  ! Aug.24,2015 - in any set of codes for a given problem, 
  ! should exclusively use Gaunt or Gcoeff to avoid getting 
  ! sign change inconsistencies

  Cval=sqrt(0.25d0*dble((2*j1+1)*(2*j2+1)*(2*j3+1))/pi)
  call threej0(j1, j2, j3, tjay)
  Cval=Cval*tjay
  ! Condon-Shortley phase Ylm*= ((-1)**m)Yl-m 
  call threej(j1, j2, j3, -m1, m2, m3, tjay)
  Gaunt=Cval*tjay*(dble(-1)**m1)

end function Gaunt
!----------------------------------------------!
subroutine threej(j1, j2, j3, m1, m2, m3, tj)

  use numbers
  implicit none
  integer j1, j2, j3, m1, m2, m3
  integer is, k, k1, k2, if1, if2
  real*8 tj
  real*8 rsum, rden, rf1, rf2


  !This routine computes the Wigner 3j symbol.
  !j1, j2, j3, m1, m2, m3 are integers
  !returns real*8 tj


  tj=0.d0
  if((m1+m2+m3).ne.0) goto 999
  if(abs(m1)>j1) goto 999
  if(abs(m2)>j2) goto 999
  if(abs(m3)>j3) goto 999
  if((j1-j2-j3)>0) goto 999
  if((abs(j2-j3)-j1)>0) goto 999
  
  k1=max(-j3+j2-m1, -j3+j1+m2, 0)
  k2=min(j1+j2-j3, j1-m1, j2+m2)
  is=(-1)**(k1)
  rsum=0
  do k=k1, k2
     rden=1.d0/rfl(k)
     rden=rden/rfl(j1+j2-j3-k)
     rden=rden/rfl(j1-m1-k)
     rden=rden/rfl(j2+m2-k)
     rden=rden/rfl(j3-j2+m1+k)
     rden=rden/rfl(j3-j1-m2+k)

     rsum=rsum+dble(is)*rden
     
     is=-is
  end do
  rf1=sqrt(rfl(j1+j2-j3)*rfl(j1-j2+j3)*rfl(-j1+j2+j3)/rfl(j1+j2+j3+1))
  rf2=sqrt(rfl(j1+m1)*rfl(j1-m1)*rfl(j2+m2)*rfl(j2-m2)*rfl(j3+m3)*rfl(j3-m3))
  is=(-1)**(j1-j2-m3)
  
  tj=rf1*rf2*dble(is)*rsum
  
999 continue

end subroutine threej
!----------------------------------------------!
!----------------------------------------------!
subroutine threej0(j1, j2, j3, tj)

  use numbers
  implicit none
  integer j1, j2, j3
  integer is, k, k1, k2, if1, if2
  real*8 tj
  real*8 rsum, rden, rf1, rf2


  !This routine computes the Wigner 3j symbol.
  !j1, j2, j3, are integers all m's are 0
  !returns real*8 tj


  tj=0.d0
  if((j1-j2-j3)>0) goto 999
  if((abs(j2-j3)-j1)>0) goto 999
  
  k1=max(-j3+j2, -j3+j1, 0)
  k2=min(j1+j2-j3, j1, j2)
  is=(-1)**(k1)
  rsum=0
  do k=k1, k2
     rden=1.d0/rfl(k)
     rden=rden/rfl(j1+j2-j3-k)
     rden=rden/rfl(j1-k)
     rden=rden/rfl(j2-k)
     rden=rden/rfl(j3-j2+k)
     rden=rden/rfl(j3-j1+k)

     rsum=rsum+dble(is)*rden
     
     is=-is
  end do
  rf1=sqrt(rfl(j1+j2-j3)*rfl(j1-j2+j3)*rfl(-j1+j2+j3)/rfl(j1+j2+j3+1))
  rf2=rfl(j1)*rfl(j2)*rfl(j3)
  is=(-1)**(j1-j2)
  
  tj=rf1*rf2*dble(is)*rsum
  
999 continue

end subroutine threej0
!----------------------------------------------!
real*8 function CG(j1, m1, j2, m2, L, M)
  implicit none

  integer j1, m1, j2, m2, L, M
  real*8 Cval
  
  call threej(j1, j2, L, m1, m2, -M, Cval)
  CG=(-1)**(M+j1-j2)*sqrt(dble(2*L+1))*Cval
end function CG
!----------------------------------------------!  
