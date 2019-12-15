!+++++++++++++++++++++++++++++++++++!
function DVReval1d(x, az, nq, nel, vq, xq, wq)

  implicit none
  complex*16 DVReval1d, x
  complex*16 vq(nel*(nq-1)-1)
  complex*16 xq(nel*(nq-1)-1), wq(nel*(nq-1)-1)
  complex*16 az(0:nel), xl(nq)
  complex*16 lsx, xm
  complex*16 accum
  integer nel, nq, nbas
  integer i, ii, iel, i1, i2

  if((x.ne.az(0)).and.(x.ne.az(nel)))then

     nbas=nel*(nq-1)-1

     iel=count((dble(az)<dble(x)))

     if((iel.ne.1).and.(iel.ne.nel))then
        i1=1
        i2=nq
        xl=xq((iel-1)*(nq-1):iel*(nq-1))
     else
        if(iel==1)then
           i1=2
           i2=nq
           xl(1)=az(0)
           xl(2:nq)=xq(1:nq-1)
        else
           i1=1
           i2=nq-1
           xl(1:nq-1)=xq(nbas-nq+2:nbas)
           xl(nq)=az(nel)     
        end if
     end if


     accum=0.d0
     ii=max((iel-1)*(nq-1), 1)
     do i=i1, i2
        lsx=product(((x-xl)/(xl(i)-xl)), mask= xl.ne.xl(i))
        accum=accum+lsx*vq(ii)/sqrt(wq(ii))
        ii=ii+1
     end do

     DVReval1d=accum

  else
     DVReval1d=0.d0
  end if

99 return
end function DVReval1d
!+++++++++++++++++++++++++++++++++++!
!+++++++++++++++++++++++++++++++++++!
function DVReval2d(x, y, az, nq, nel, vq, xq, wq)

  implicit none
  complex*16 DVReval2d, x, y, z
  complex*16 vq(nel*(nq-1)-1, nel*(nq-1)-1)
  complex*16 xq(nel*(nq-1)-1), wq(nel*(nq-1)-1)
  complex*16 az(0:nel), xl(nq), yl(nq)
  complex*16 lsx, lsy
  complex*16 accum
  integer nel, nq, nbas
  integer i, ii, iel, i1, i2
  integer j, jj, jel, j1, j2

  if((x.ne.az(0)).and.(x.ne.az(nel)).and.(y.ne.az(0)).and.(y.ne.az(nel)))then


     nbas=nel*(nq-1)-1

     iel=count((dble(az)<dble(x)))
     jel=count((dble(az)<dble(y)))

     if((iel.ne.1).and.(iel.ne.nel))then
        i1=1
        i2=nq
        xl=xq((iel-1)*(nq-1):iel*(nq-1))
     else
        if(iel==1)then
           i1=2
           i2=nq
           xl(1)=az(0)
           xl(2:nq)=xq(1:nq-1)
        else
           i1=1
           i2=nq-1
           xl(1:nq-1)=xq(nbas-nq+2:nbas)
           xl(nq)=az(nel)     
        end if
     end if
     if((jel.ne.1).and.(jel.ne.nel))then
        j1=1
        j2=nq
        yl=xq((jel-1)*(nq-1):jel*(nq-1))
     else
        if(jel==1)then
           j1=2
           j2=nq
           yl(1)=az(0)
           yl(2:nq)=xq(1:nq-1)
        else
           j1=1
           j2=nq-1
           yl(1:nq-1)=xq(nbas-nq+2:nbas)
           yl(nq)=az(nel)     
        end if
     end if


     accum=0.d0
     ii=max((iel-1)*(nq-1), 1)
     do i=i1, i2
        lsx=product(((x-xl)/(xl(i)-xl)), mask= xl.ne.xl(i))
        jj=max((jel-1)*(nq-1), 1)
        do j=j1, j2
           lsy=product(((y-yl)/(yl(j)-yl)), mask= yl.ne.yl(j))
           accum=accum+lsx*lsy*vq(ii, jj)/sqrt(wq(ii)*wq(jj))
           jj=jj+1
        end do
        ii=ii+1
     end do
     DVReval2d=accum

  else
     DVReval2d=0.d0
  end if

99 return
end function DVReval2d
!+++++++++++++++++++++++++++++++++++!
!+++++++++++++++++++++++++++++++++++!
function DVReval3d(x, y, z, az, nq, nel, vq, xq, wq)

  implicit none
  complex*16 DVReval3d, x, y, z
  complex*16 vq(nel*(nq-1)-1, nel*(nq-1)-1,  nel*(nq-1)-1)
  complex*16 xq(nel*(nq-1)-1), wq(nel*(nq-1)-1)
  complex*16 az(0:nel), xl(nq), yl(nq), zl(nq)
  complex*16 lsx, lsy, lsz
  complex*16 accum
  integer nel, nq, nbas
  integer i, ii, iel, i1, i2
  integer j, jj, jel, j1, j2
  integer k, kk, kel, k1, k2


  if((x.ne.az(0)).and.(x.ne.az(nel)).and.(y.ne.az(0)).or.(y.ne.az(nel)).and.(z.ne.az(0)).and.(z.ne.az(nel)))then

     nbas=nel*(nq-1)-1

     iel=count((dble(az)<dble(x)))
     jel=count((dble(az)<dble(y)))
     kel=count((dble(az)<dble(z)))

     if((iel.ne.1).and.(iel.ne.nel))then
        i1=1
        i2=nq
        xl=xq((iel-1)*(nq-1):iel*(nq-1))
     else
        if(iel==1)then
           i1=2
           i2=nq
           xl(1)=az(0)
           xl(2:nq)=xq(1:nq-1)
        else
           i1=1
           i2=nq-1
           xl(1:nq-1)=xq(nbas-nq+2:nbas)
           xl(nq)=az(nel)     
        end if
     end if
     if((jel.ne.1).and.(jel.ne.nel))then
        j1=1
        j2=nq
        yl=xq((jel-1)*(nq-1):jel*(nq-1))
     else
        if(jel==1)then
           j1=2
           j2=nq
           yl(1)=az(0)
           yl(2:nq)=xq(1:nq-1)
        else
           j1=1
           j2=nq-1
           yl(1:nq-1)=xq(nbas-nq+2:nbas)
           yl(nq)=az(nel)     
        end if
     end if
     if((kel.ne.1).and.(kel.ne.nel))then
        k1=1
        k2=nq
        zl=xq((kel-1)*(nq-1):kel*(nq-1))
     else
        if(kel==1)then
           k1=2
           k2=nq
           zl(1)=az(0)
           zl(2:nq)=xq(1:nq-1)
        else
           k1=1
           k2=nq-1
           zl(1:nq-1)=xq(nbas-nq+2:nbas)
           zl(nq)=az(nel)     
        end if
     end if

     accum=0.d0
     ii=max((iel-1)*(nq-1), 1)
     do i=i1, i2
        lsx=product(((x-xl)/(xl(i)-xl)), mask= xl.ne.xl(i))
        jj=max((jel-1)*(nq-1), 1)
        do j=j1, j2
           lsy=product(((y-yl)/(yl(j)-yl)), mask= yl.ne.yl(j))
           kk=max((kel-1)*(nq-1), 1)
           do k=k1, k2
              lsz=product(((z-zl)/(zl(k)-zl)), mask= zl.ne.zl(k))
              accum=accum+lsx*lsy*lsz*vq(ii, jj, kk)/sqrt(wq(ii)*wq(jj)*wq(kk))
              kk=kk+1
           end do
           jj=jj+1
        end do
        ii=ii+1
     end do
     DVReval3d=accum

  else
     DVReval3d=0.d0
  end if


99 return
end function DVReval3d
!+++++++++++++++++++++++++++++++++++!
!-----------------------------------!
subroutine lsf_z(x, m, n, xdvr, fls)

  !f_m(x)
  implicit none
  integer m, n
  complex*16 x, xdvr(n), fls, xm, xi
  integer i

  xm=xdvr(m)
  fls=product(((x-xdvr)/(xm-xdvr)), mask= xdvr.ne.xm)

end subroutine lsf_z
!-----------------------------------!
!-----------------------------------!
subroutine lsfp_z(x, m, n, xdvr, flsp)

  implicit none
  integer m, n
  complex*16 x, xdvr(n), flsp, pr
  integer p, i


  flsp=0.d0
  do p=1, m-1
     pr=1.d0
     do i=1, p-1
        pr=pr*((x-xdvr(i))/(xdvr(m)-xdvr(i)))
     end do
     do i=p+1, m-1
        pr=pr*((x-xdvr(i))/(xdvr(m)-xdvr(i)))
     end do
     do i=m+1, n
        pr=pr*((x-xdvr(i))/(xdvr(m)-xdvr(i)))
     end do
     flsp=flsp+pr/(xdvr(m)-xdvr(p))
  end do
  do p=m+1, n
     pr=1.d0
     do i=1, m-1
        pr=pr*((x-xdvr(i))/(xdvr(m)-xdvr(i)))
     end do
     do i=m+1, p-1
        pr=pr*((x-xdvr(i))/(xdvr(m)-xdvr(i)))
     end do
     do i=p+1, n
        pr=pr*((x-xdvr(i))/(xdvr(m)-xdvr(i)))
     end do
     flsp=flsp+pr/(xdvr(m)-xdvr(p))
  end do

end subroutine lsfp_z
!-----------------------------------!
