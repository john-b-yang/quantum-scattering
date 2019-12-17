!+++++++++++++++++++++++++++++++++++!
complex*16 function DVReval1dR(x, vq)

  ! produces the function value at x of a DVR-Radau vector as defined in
  ! Bill's notes (Sept. 2017) in Equation 15 (and Lobatto FEM-DVR otherwise)

  ! note::"fac" on the Radau tail must be combatted with the weight (e^-x)
  ! in a quadrature, where x=(r-R0)*2alphaRad and r is a quadarature point 

  ! note::az runs up to R0 at (nel-1) and az(nel) is last(complex) Radau point

  use dvrecsRad, only: nel, rel, nq, mq, nbas, az, xz, wz, eit, R0, alphaRad
  implicit none
  complex*16, intent(in) ::  vq(nbas), x
  !      xl will hold either Lobatto or Radau pts
  complex*16 xl(nq+mq), lsx, accum, fac 
  integer i, ii, iel, i1, i2, iq

  DVReval1dR=0.d0

  if((dble(x)>dble(az(0))).and.(dble(x)<=dble(az(nel)))) then

     iel=count((dble(az)<dble(x))) ! tells which FEM the point x lies in
     ii=(iel-1)*(nq-1)
     fac=1.d0 ! fac accounts for the weight of Lobatto (1.) or Radau quad (e^-x)

!***     xl=1.d0 ! need to define this way for the product when using mask below
     if((iel.ne.1).and.(iel.ne.nel))then !interior FEM (Lobatto)
        i1=1
        i2=nq
        xl(1:nq)=xz((iel-1)*(nq-1):iel*(nq-1))
     else
        if(iel==1)then  ! first FEM (Lobatto always) 
           i1=2
           i2=nq
           xl(1)=az(0)
           xl(2:nq)=xz(1:nq-1)
           ii=1
        else if (mq.eq.0) then   ! last FEM on Lobatto-only grid
           i1=1
           i2=nq-1
           xl(1:nq-1)=xz((iel-1)*(nq-1):iel*(nq-1))
           xl(nq)=az(nel)  
        else            ! last FEM (Radau) 
           i1=1
           i2=mq
           xl(1:mq)=xz(nbas-mq+1:nbas)
           ii=nbas-mq+1
           ! Radau weight factor on ECS (Eq. 15 of Bill's notes)
           fac=exp(-alphaRad*conjg(eit)*(x-R0)) 
        end if
     end if

     accum=0.d0
     do i=i1, i2
!*** alternative way of finding the product with mask (xl initialization above)
!        lsx=product(((x-xl)/(xl(i)-xl)), mask= xl.ne.xl(i))
!*** loop below may be more efficient (only runs over nonzero functions in FEM)
        lsx=1.d0
        do iq=i1, i2
           if(iq.eq.i) cycle
           lsx=lsx*((x-xl(iq))/(xl(i)-xl(iq)))
        end do
        accum=accum+lsx*vq(ii)/sqrt(wz(ii))
        ii=ii+1
     end do

     DVReval1dR=accum*fac

  end if

99 return
end function DVReval1dR
!+++++++++++++++++++++++++++++++++++!
!+++++++++++++++++++++++++++++++++++!
complex*16 function DVReval2dR(x, y, vq)

  ! produces the function value at x of a DVR-Radau vector as defined in
  ! Bill's notes (Sept. 2017) in Equation 15 (and Lobatto FEM-DVR otherwise)

  ! note::"fac" on the Radau tail must be combatted with the weight (e^-x)
  ! in a quadrature, where x=(r-R0)*2alphaRad and r is a quadarature point 

  ! note::az runs up to R0 at (nel-1) and az(nel) is last(complex) Radau point

  use dvrecsRad, only: nel, rel, nq, mq, nbas, az, xz, wz, eit, R0, alphaRad
  implicit none
  complex*16, intent(in) ::  vq(nbas, nbas), x, y
  !      xl and yl will hold either Lobatto or Radau pts
  complex*16 xl(nq+mq), yl(nq+mq), lsx, lsy, accum, faci, facj
  integer i, ii, iel, i1, i2, iq, iis
  integer j, jj, jel, j1, j2, jq, jjs

  DVReval2dR=0.d0

  if((dble(x)>dble(az(0))).and.(dble(x)<=dble(az(nel))).and.&
       (dble(y)>dble(az(0))).and.(dble(y)<=dble(az(nel)))) then

     iel=count((dble(az)<dble(x))) ! tells which FEM the points x, y lie in
     jel=count((dble(az)<dble(y)))
     iis=(iel-1)*(nq-1)
     jjs=(jel-1)*(nq-1)
     faci=1.d0 !fac accounts for the weight of Lobatto (1.) or Radau quad (e^-x)
     facj=1.d0 

     if((iel.ne.1).and.(iel.ne.nel))then !interior FEM (Lobatto)
        i1=1
        i2=nq
        xl(1:nq)=xz((iel-1)*(nq-1):iel*(nq-1))
     else
        if(iel==1)then ! first FEM (Lobatto) 
           i1=2
           i2=nq
           xl(1)=az(0)
           xl(2:nq)=xz(1:nq-1)
           iis=1
        else if (mq.eq.0) then   ! last FEM on Lobatto-only grid  
           i1=1
           i2=nq-1
           xl(1:nq-1)=xz((iel-1)*(nq-1):iel*(nq-1))
           xl(nq)=az(nel)
        else            ! last FEM (Radau) 
           i1=1
           i2=mq
           xl(1:mq)=xz(nbas-mq+1:nbas)
           iis=nbas-mq+1
           ! Radau weight factor on ECS (Eq. 15 of Bill's notes)
           faci=exp(-alphaRad*conjg(eit)*(x-R0)) 
        end if
     end if
     if((jel.ne.1).and.(jel.ne.nel))then  !interior FEM (Lobatto)
        j1=1
        j2=nq
        yl(1:nq)=xz((jel-1)*(nq-1):jel*(nq-1))
     else
        if(jel==1)then ! first FEM (Lobatto) 
           j1=2
           j2=nq
           yl(1)=az(0)
           yl(2:nq)=xz(1:nq-1)
           jjs=1
        else if (mq.eq.0) then   ! last FEM on Lobatto-only grid
           j1=1
           j2=nq-1
           yl(1:nq-1)=xz((jel-1)*(nq-1):jel*(nq-1))
           yl(nq)=az(nel)
        else            ! last FEM (Radau) 
           j1=1
           j2=mq
           yl(1:mq)=xz(nbas-mq+1:nbas)
           jjs=nbas-mq+1
           ! Radau weight factor on ECS (Eq. 15 of Bill's notes)
           facj=exp(-alphaRad*conjg(eit)*(y-R0))            
        end if
     end if

     accum=0.d0
     ii=iis
     do i=i1, i2
        !lsx=product(((x-xl)/(xl(i)-xl)), mask= xl.ne.xl(i))
        lsx=1.d0
        do iq=i1, i2
           if(iq.eq.i) cycle
           lsx=lsx*((x-xl(iq))/(xl(i)-xl(iq)))
        end do
        jj=jjs
        do j=j1, j2
           !lsy=product(((y-yl)/(yl(j)-yl)), mask= yl.ne.yl(j))
           lsy=1.d0
           do jq=j1, j2
              if(jq.eq.j) cycle
              lsy=lsy*((y-yl(jq))/(yl(j)-yl(jq)))
           end do
           accum=accum+lsx*lsy*vq(ii, jj)/sqrt(wz(ii)*wz(jj))
           jj=jj+1
        end do
        ii=ii+1
     end do

     DVReval2dR=accum*faci*facj

  end if

99 return
end function DVReval2dR
!+++++++++++++++++++++++++++++++++++!
!+++++++++++++++++++++++++++++++++++!
!!$!!! UNTESTED ...!!! UNTESTED ...!!! UNTESTED ...!!! UNTESTED ..
!!$!!! 3D-DVR function HASN'T BEEN TESTED SINCE written (June 2018)
!!$!!! delete UNTESTED ...!!! AFTER CHECKING THE BELOW CAREFULLY (NEVER DONE!!)
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!!$complex*16 function DVReval3dR(x, y, z, vq)
!!$
!!$!!! UNTESTED ...!!! UNTESTED ...!!! UNTESTED ...!!! UNTESTED ...
!!$!!! UNTESTED ...!!! UNTESTED ...!!! UNTESTED ...!!! UNTESTED (written) June2018
!!$!!! UNTESTED ...!!! UNTESTED ...!!! UNTESTED ...!!! UNTESTED ...
!!$
!!$  ! produces the function value at x of a DVR-Radau vector as defined in
!!$  ! Bill's notes (Sept. 2017) in Equation 15 (and Lobatto FEM-DVR otherwise)
!!$
!!$  ! note::"fac" on the Radau tail must be combatted with the weight (e^-x)
!!$  ! in a quadrature, where x=(r-R0)*2alphaRad and r is a quadarature point 
!!$
!!$  ! note::az runs up to R0 at (nel-1) and az(nel) is last(complex) Radau point
!!$
!!$  use dvrecsRad, only: nel, rel, nq, mq, nbas, az, xz, wz, eit, R0, alphaRad
!!$  implicit none
!!$  complex*16, intent(in) ::  vq(nbas, nbas, nbas), x, y, z
!!$  !      xl yl zl will hold either Lobatto or Radau pts
!!$  complex*16 xl(nq+mq),yl(nq+mq),zl(nq+mq), lsx,lsy,lsz, accum, faci, facj, fack
!!$  integer i, ii, iel, i1, i2, iq, iis
!!$  integer j, jj, jel, j1, j2, jq, jjs
!!$  integer k, kk, kel, k1, k2, kq, kks
!!$
!!$  DVReval3dR=0.d0
!!$  if((dble(x)>dble(az(0))).and.(dble(x)<=dble(az(nel))).and.&
!!$       (dble(y)>dble(az(0))).and.(dble(y)<=dble(az(nel))).and.&
!!$       (dble(z)>dble(az(0))).and.(dble(z)<=dble(az(nel)))) then
!!$
!!$     iel=count((dble(az)<dble(x)))! tells which FEM the points x, y, z lie in
!!$     jel=count((dble(az)<dble(y)))
!!$     kel=count((dble(az)<dble(z)))
!!$     iis=(iel-1)*(nq-1)
!!$     jjs=(jel-1)*(nq-1)
!!$     kks=(kel-1)*(nq-1)
!!$     faci=1.d0 !fac accounts for the weight of Lobatto (1.) or Radau quad (e^-x)
!!$     facj=1.d0
!!$     fack=1.d0
!!$     
!!$     if((iel.ne.1).and.(iel.ne.nel))then !interior FEM (Lobatto)
!!$        i1=1
!!$        i2=nq
!!$        xl(1:nq)=xz((iel-1)*(nq-1):iel*(nq-1))
!!$     else
!!$        if(iel==1)then ! first FEM (Lobatto) 
!!$           i1=2
!!$           i2=nq
!!$           xl(1)=az(0)
!!$           xl(2:nq)=xz(1:nq-1)
!!$           iis=1
!!$        else if (mq.eq.0) then   ! last FEM on Lobatto-only grid
!!$           i1=1
!!$           i2=nq-1
!!$           xl(1:nq-1)=xz((iel-1)*(nq-1):iel*(nq-1))
!!$           xl(nq)=az(nel)
!!$        else            ! last FEM (Radau) 
!!$           i1=1
!!$           i2=mq
!!$           xl(1:mq)=xz(nbas-mq+1:nbas)
!!$           iis=nbas-mq+1
!!$           ! Radau weight factor on ECS (Eq. 15 of Bill's notes)
!!$           faci=exp(-alphaRad*conjg(eit)*(x-R0)) 
!!$        end if
!!$     end if
!!$     if((jel.ne.1).and.(jel.ne.nel))then  !interior FEM (Lobatto)
!!$        j1=1
!!$        j2=nq
!!$        yl(1:nq)=xz((jel-1)*(nq-1):jel*(nq-1))
!!$     else
!!$        if(jel==1)then ! first FEM (Lobatto) 
!!$           j1=2
!!$           j2=nq
!!$           yl(1)=az(0)
!!$           yl(2:nq)=xz(1:nq-1)
!!$           jjs=1
!!$        else if (mq.eq.0) then   ! last FEM on Lobatto-only grid
!!$           j1=1
!!$           j2=nq-1
!!$           yl(1:nq-1)=xz((jel-1)*(nq-1):jel*(nq-1))
!!$           yl(nq)=az(nel)
!!$        else            ! last FEM (Radau) 
!!$           j1=1
!!$           j2=mq
!!$           yl(1:mq)=xz(nbas-mq+1:nbas)
!!$           jjs=nbas-mq+1
!!$           ! Radau weight factor on ECS (Eq. 15 of Bill's notes)
!!$           facj=exp(-alphaRad*conjg(eit)*(y-R0))            
!!$        end if
!!$     end if
!!$     if((kel.ne.1).and.(kel.ne.nel))then  !interior FEM (Lobatto)
!!$        k1=1
!!$        k2=nq
!!$        zl(1:nq)=xz((kel-1)*(nq-1):kel*(nq-1))
!!$     else
!!$        if(kel==1)then ! first FEM (Lobatto) 
!!$           k1=2
!!$           k2=nq
!!$           zl(1)=az(0)
!!$           zl(2:nq)=xz(1:nq-1)
!!$           kks=1
!!$        else if (mq.eq.0) then   ! last FEM on Lobatto-only grid
!!$           k1=1
!!$           k2=nq-1
!!$           zl(1:nq-1)=xz((kel-1)*(nq-1):kel*(nq-1))
!!$           zl(nq)=az(nel)
!!$        else            ! last FEM (Radau) 
!!$           k1=1
!!$           k2=mq
!!$           zl(1:mq)=xz(nbas-mq+1:nbas)
!!$           kks=nbas-mq+1
!!$           ! Radau weight factor on ECS (Eq. 15 of Bill's notes)
!!$           fack=exp(-alphaRad*conjg(eit)*(z-R0))            
!!$        end if
!!$     end if
!!$
!!$     accum=0.d0
!!$     ii=iis
!!$     do i=i1, i2
!!$        !lsx=product(((x-xl)/(xl(i)-xl)), mask= xl.ne.xl(i))
!!$        lsx=1.d0
!!$        do iq=i1, i2
!!$           if(iq.eq.i) cycle
!!$           lsx=lsx*((x-xl(iq))/(xl(i)-xl(iq)))
!!$        end do
!!$        jj=jjs
!!$        do j=j1, j2
!!$           !lsy=product(((y-yl)/(yl(j)-yl)), mask= yl.ne.yl(j))
!!$           lsy=1.d0
!!$           do jq=j1, j2
!!$              if(jq.eq.j) cycle
!!$              lsy=lsy*((y-yl(jq))/(yl(j)-yl(jq)))
!!$           end do
!!$           kk=kks
!!$           do k=k1, k2
!!$              !lsz=product(((z-zl)/(zl(k)-zl)), mask= zl.ne.zl(k))
!!$              lsz=1.d0
!!$              do kq=k1, k2
!!$                 if(kq.eq.k) cycle
!!$                 lsz=lsz*((z-zl(kq))/(zl(k)-zl(kq)))
!!$              end do
!!$              accum=accum+lsx*lsy*lsz*vq(ii, jj, kk)/sqrt(wz(ii)*wz(jj)*wz(kk))
!!$              kk=kk+1
!!$           end do
!!$           jj=jj+1
!!$        end do
!!$        ii=ii+1
!!$     end do
!!$     
!!$     DVReval3dR=accum*faci*facj*fack
!!$
!!$  end if
!!$
!!$99 return
!!$end function DVReval3dR
!!$!+++++++++++++++++++++++++++++++++++!
!!$!-----------------------------------!
!!$subroutine lsf_z(x, m, n, xdvr, fls)
!!$
!!$  !f_m(x)
!!$  implicit none
!!$  integer m, n
!!$  complex*16 x, xdvr(n), fls, xm, xi
!!$  integer i
!!$
!!$  xm=xdvr(m)
!!$  fls=product(((x-xdvr)/(xm-xdvr)), mask= xdvr.ne.xm)
!!$
!!$end subroutine lsf_z
!!$!-----------------------------------!
!!$!-----------------------------------!
!!$subroutine lsfp_z(x, m, n, xdvr, flsp)
!!$
!!$  implicit none
!!$  integer m, n
!!$  complex*16 x, xdvr(n), flsp, pr
!!$  integer p, i
!!$
!!$
!!$  flsp=0.d0
!!$  do p=1, m-1
!!$     pr=1.d0
!!$     do i=1, p-1
!!$        pr=pr*((x-xdvr(i))/(xdvr(m)-xdvr(i)))
!!$     end do
!!$     do i=p+1, m-1
!!$        pr=pr*((x-xdvr(i))/(xdvr(m)-xdvr(i)))
!!$     end do
!!$     do i=m+1, n
!!$        pr=pr*((x-xdvr(i))/(xdvr(m)-xdvr(i)))
!!$     end do
!!$     flsp=flsp+pr/(xdvr(m)-xdvr(p))
!!$  end do
!!$  do p=m+1, n
!!$     pr=1.d0
!!$     do i=1, m-1
!!$        pr=pr*((x-xdvr(i))/(xdvr(m)-xdvr(i)))
!!$     end do
!!$     do i=m+1, p-1
!!$        pr=pr*((x-xdvr(i))/(xdvr(m)-xdvr(i)))
!!$     end do
!!$     do i=p+1, n
!!$        pr=pr*((x-xdvr(i))/(xdvr(m)-xdvr(i)))
!!$     end do
!!$     flsp=flsp+pr/(xdvr(m)-xdvr(p))
!!$  end do
!!$
!!$end subroutine lsfp_z
!!$!-----------------------------------!
