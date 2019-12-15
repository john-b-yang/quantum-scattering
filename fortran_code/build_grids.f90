!-----------------------------------!
subroutine build_grids()

  use dvrecs
  implicit none
  integer,parameter:: debug=2 !For debugging
  complex*16 xs_l(nq), ws_l(nq), xs_r(mq), ws_r(mq)!Temps 
  real*8 lag_alpha
  complex*16 lscal
  integer i, j, ii, p

  allocate(xz(nbas), wz(nbas))

  !create nth order DVR on [-1, 1]
  allocate(xq_l(nq), wq_l(nq))!Raw points and weights
  call gl_quad(nq, xq_l, wq_l)
  allocate(xqz(nq))
  xqz=xq_l
  !Write raw points and weights if desired
  if(debug==1) then
      write(unit_stdout,*) 'Raw DVR points and weights on -1,1'
      write(unit_stdout,'(6f11.6)')(xq_l(i),i=1,nq)
      write(unit_stdout,*)
      write(unit_stdout,'(6f11.6)')(wq_l(i),i=1,nq)
      write(unit_stdout,*)
    endif

  !create nth order Radau DVR on [0, Inf]
  allocate(xq_r(mq), wq_r(mq))!Raw points and weights
  call gr_quad(mq, xq_r, wq_r)
  !Write raw points and weights if desired
  if(debug==2) then
    !200 format (F10.40)
    !open(701, file='Laguerre_Radau', form='formatted', status='unknown')
    open(701, file='Laguerre_Radau', status='unknown')
      !write(unit_stdout,*) 'Raw DVR points and weights on 0,inf'
      !write(701,'(6f11.6)')(xq_r(i),i=1,mq)
      do i=1,mq
        write(701,*) xq_r(i)
        !write(701,200) xq_r(i)
      enddo
      !write(701,'(6f11.6)')(wq_r(i),i=1,mq)
        write(701,*)
      do i=1,mq
        write(701,*) wq_r(i)
        !write(701,200) wq_r(i)
      enddo
    endif

    lag_alpha=0.3d0 !Scrinzi says this was a good choice.
    lscal=eit/(2*lag_alpha) !Transformation to Laguerre boundaries on a complex ray

!*******************************************
!Creating Grid from raw points--------------
!*******************************************
if(rel==nel)then
    p=rel-1
else
    p=rel
endif
ii=0
do i=1, p !Labatto body
    xs_l=az(i-1)+0.5d0*(xq_l+1)*(az(i)-az(i-1))
    ws_l=0.5*wq_l*(az(i)-az(i-1))
    if(ii.ne.0) wz(ii)=wz(ii)+ws_l(1)
    do j=2, nq
        ii=ii+1
        xz(ii)=xs_l(j)  
        wz(ii)=ws_l(j)
    end do
end do

if(rel==nel)then
    xs_l=az(i-1)+0.5d0*(xq_l+1)*(az(i)-az(i-1))
    ws_l=0.5*wq_l*(az(i)-az(i-1))
    if(ii.ne.0) wz(ii)=wz(ii)+ws_l(1)
    do j=2, nq-1
        ii=ii+1
        xz(ii)=xs_l(j)  
        wz(ii)=ws_l(j)
    end do
else !Radau tail
    xs_r=lscal*xq_r+R0
    ws_r=lscal*wq_r
    wz(ii)=wz(ii)+ws_r(1)
    do j=2, mq-1
        ii=ii+1
        xz(ii)=xs_r(j)
        wz(ii)=ws_r(j)
    end do
end if

!do i=1,mq
    !print*,xs_r(i)
!enddo
do i=1,nbas-1
    print*,xz(i)
enddo

end subroutine build_grids
!-----------------------------------!
