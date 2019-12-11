!-----------------------------------!
subroutine build_dq_nw(n, d, quad)
  
  use dvrecs
  implicit none
  integer n, i, j, k, quad
  real*8 d(n, n), x(n), w(n), pr

  
  if (quad==1) then !Evaluate Labatto derivative
    call gl_quad(n, x, w)
  else
    call gr_quad(n, x, w) !Evaluate Radau derivative
  end if

        do i=1, n
            d(i, i)=0.d0
            do k=1, n
                if(i.ne.k) d(i, i)=d(i, i)+1.d0/(x(i)-x(k))
            end do
            do j=1, n
                if(i.ne.j)then
                pr=1.d0/(x(j)-x(i))
                do k=1, n
                    if((k.ne.i).and.(k.ne.j)) pr=pr*(x(i)-x(k))/(x(j)-x(k))
                end do
                d(i, j)=pr
                end if
            end do
        end do

end subroutine build_dq_nw
!-----------------------------------!
