!-----------------------------------!
subroutine gl_quad(n, t, w)

  implicit none
  real*8 t(n), w(n), scr(n)
  real*8 alpha, beta, endpts(2)
  integer n, kind, kpts
  
 !generate Gauss-Lobatto quadrature points and weights
  alpha=0.d0
  beta=0.d0
  endpts(1)=-1.d0
  endpts(2)=1.d0
  kind=1
  kpts=2
  call gaussq(kind, n, alpha, beta, kpts, endpts, scr, t, w)
  t(1)=-1.d0
  t(n)=1.d0

end subroutine gl_quad
!-----------------------------------!
!-----------------------------------!
subroutine gr_quad(n, t, w)

  implicit none
  real*8 t(n), w(n), scr(n)
  real*8 alpha, beta, endpts(2)
  integer n, kind, kpts
  
 !generate Gauss-Lobatto quadrature points and weights
  alpha=0.d0
  beta=0.d0
  endpts(1)=0.d0
  endpts(2)=1.d0
  kind=6 !Laguerre [0,inf]
  kpts=1 !left endpoint is fixed
  call gaussq(kind, n, alpha, beta, kpts, endpts, scr, t, w)

end subroutine gr_quad
!-----------------------------------!
