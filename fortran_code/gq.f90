!-----------------------------------!
subroutine gaussq(kind, n, alpha, beta, kpts, endpts, b, t, w)

  implicit none
  real*8 muzero, alpha, beta
  real*8 gam, t1, gbslve, h
  real*8 b(n), t(n), w(n), endpts(2)
  integer kind, n, kpts, ierr, nm1, i
     
  if(kind.eq.0) then
     if(2*(n/2).eq.n) then
        write(6,800) n
800     format(" n must be odd for simpson's rule ",i5)
        stop
     end if
     if(n.le.1) then
        t(1) = 0.
        w(1) = 2.0
        return
     end if
     h = 2.0/(n-1)
     t(1) = -1.0
     t(n) = 1.0
     w(1) = h/3.0
     w(n) = h/3.0
     nm1 = n-1
     do  i=2,nm1
        t(i) = t(i-1) + h
        w(i) = 4.0 - 2.0*(i-2*(i/2))
        w(i) = w(i)*h/3.0
     end do
     return
  end if

  call class (kind, n, alpha, beta, b, t, muzero)
  

!           the matrix of coefficients is assumed to be symmetric.
!           the array t contains the diagonal elements, the array
!           b the off-diagonal elements.
!           make appropriate changes in the lower right 2 by 2
!           submatrix.

  if (kpts.eq.0)  go to 100
  if (kpts.eq.2)  go to  50

!           if kpts=1, only t(n) must be changed

  t(n)=gbslve(endpts(1), n, t, b)*b(n-1)**2 + endpts(1)
  go to 100

!           if kpts=2, t(n) and b(n-1) must be recomputed

50 gam=gbslve(endpts(1), n, t, b)
  t1 = ((endpts(1) - endpts(2))/(gbslve(endpts(2), n, t, b) - gam))
  b(n-1) =  sqrt(t1)
  t(n) = endpts(1) + gam*t1
  
!           note that the indices of the elements of b run from 1 to n-1
!           and thus the value of b(n) is arbitrary.
!           now compute the eigenvalues of the symmetric tridiagonal
!           matrix, which has been modified as necessary.
!           the method used is a ql-type method with origin shifting

100 w(1) = 1.0d0
  do i = 2, n
     w(i) = 0.0d0
  end do
     
  call gbtql2 (n, t, b, w, ierr)
  do i = 1, n
     w(i) = muzero * w(i) * w(i)
  end do
  return
end subroutine gaussq
!-----------------------------------!
!+++++++++++++++++++++++++++++++++++!
function gbslve(shift, n, a, b)

  implicit none
  real*8 a(n),b(n)
  real*8 alpha, shift, gbslve
  integer i, nm1, n

  alpha=a(1)-shift
  nm1=n-1
  do i=2, nm1
     alpha=a(i)-shift-b(i-1)**2/alpha
  end do
  gbslve=1.0d0/alpha
  return

end function gbslve
!+++++++++++++++++++++++++++++++++++!
!-----------------------------------!
subroutine class(kind, n, alpha, beta, b, a, muzero)

!
!           this procedure supplies the coefficients a(j), b(j) of the
!        recurrence relation
!
!             b p (x) = (x - a ) p   (x) - b   p   (x)
!              j j            j   j-1       j-1 j-2
!
!        for the various classical (normalized) orthogonal polynomials,
!        and the zero-th moment
!
!             muzero = integral w(x) dx
!
!        of the given polynomial   weight function w(x).  since the
!        polynomials are orthonormalized, the tridiagonal matrix is
!        guaranteed to be symmetric.
!
!           the input parameter alpha is used only for laguerre and
!        jacobi polynomials, and the parameter beta is used only for
!        jacobi polynomials.  the laguerre and jacobi polynomials
!        require the gamma function.
!
!     ..................................................................
!

!      implicit real*8 (a-h,o-z)

  implicit none
  dimension  a(n),b(n)
  real*8 muzero, pi, ab, abi, a2b2, gamfun, a,  alpha
  real*8 b, beta
  integer n, nm1, i, kind
  data pi / 3.141592653589793d0  /


  nm1 = n - 1
  go to (10, 20, 30, 40, 50, 60), kind

!              kind = 1=  legendre polynomials p(x)
!              on (-1, +1), w(x) = 1.
!
10 muzero = 2.0d0
  do i = 1, nm1
     a(i) = 0.0d0
     abi = i
     b(i) = abi/ sqrt(4*abi*abi - 1.0d0  )
  end do
  a(n) = 0.0d0
  return

!
!              kind = 2=  chebyshev polynomials of the first kind t(x)
!              on (-1, +1), w(x) = 1 / sqrt(1 - x*x)
!
20 muzero = pi
  do i = 1, nm1
     a(i) = 0.0d0
     b(i) = 0.5d0
  end do
  b(1) =  sqrt(0.5d0  )
  a(n) = 0.0d0
  return

!
!              kind = 3=  chebyshev polynomials of the second kind u(x)
!              on (-1, +1), w(x) = sqrt(1 - x*x)
!
30 muzero = pi/2.0d0
  do i = 1, nm1
     a(i) = 0.0d0
     b(i) = 0.5d0
  end do
  a(n) = 0.0d0
  return
  
!
!              kind = 4=  hermite polynomials h(x) on (-infinity,
!              +infinity), w(x) = exp(-x**2)
!
40 muzero =  sqrt(pi)
  do i = 1, nm1
     a(i) = 0.0d0
     b(i) =  sqrt(i/2.0d0  )
  end do
  a(n) = 0.0d0
  return
!
!               kind = 5=  jacobi polynomials p(alpha, beta)(x) on
!              (-1, +1), w(x) = (1-x)**alpha + (1+x)**beta, alpha and
!              beta greater than -1
!
50 ab=alpha+beta
  abi=2.0d0+ab
  muzero=2.0d0**(ab+1.0d0)*gamfun(alpha+1.0d0)*gamfun(beta+1.0d0)/gamfun(abi)
  a(1)=(beta - alpha)/abi
  b(1)=sqrt(4.0d0  *(1.0d0  + alpha)*(1.0d0   + beta)/((abi + 1.0d0)*abi*abi))
  a2b2 = beta*beta - alpha*alpha
  do i = 2, nm1
     abi = 2.0d0  *i + ab
     a(i) = a2b2/((abi - 2.0d0  )*abi)
     b(i) =sqrt(4.0d0*i*(i+alpha)*(i+beta)*(i+ab)/((abi*abi-1)*abi*abi))
  end do
  abi = 2.0d0  *n + ab
  a(n) = a2b2/((abi - 2.0d0  )*abi)
  return

!
!              kind = 6=  laguerre polynomials l(alpha)(x) on
!              (0, +infinity), w(x) = exp(-x) * x**alpha, alpha greater
!              than -1.
!
60 muzero = gamfun(alpha + 1.0d0  )
  do i = 1, nm1
     a(i) = 2.0d0  *i - 1.0d0   + alpha
     b(i) =  sqrt(i*(i + alpha))
  end do
  a(n) = 2.0d0  *n - 1 + alpha
  return

end subroutine class
!-----------------------------------!
!-----------------------------------!
subroutine gbtql2(n, d, e, z, ierr)

  !     this subroutine is a translation of the algol procedure imtql2,
  !     num. math. 12, 377-383(1968) by martin and wilkinson,
  !     as modified in num. math. 15, 450(1970) by dubrulle.
  !     handbook for auto. comp., vol.ii-linear algebra, 241-248(1971).
  !
  !     this subroutine finds the eigenvalues and first components of the
  !     eigenvectors of a symmetric tridiagonal matrix by the implicit ql
  !     method, and is adapted from the eispak routine imtql2
  !
  !     on input=
  !
  !        n is the order of the matrix;
  !
  !        d contains the diagonal elements of the input matrix;
  !
  !        e contains the subdiagonal elements of the input matrix
  !          in its first n-1 positions.  e(n) is arbitrary;
  !
  !        z contains the first row of the identity matrix.
  !
  !      on output=
  !
  !        d contains the eigenvalues in ascending order.  if an
  !          error exit is made, the eigenvalues are correct but
  !          unordered for indices 1, 2, ..., ierr-1;
  !
  !        e has been destroyed;
  !
  !        z contains the first components of the orthonormal eigenvectors
  !          of the symmetric tridiagonal matrix.  if an error exit is
  !          made, z contains the eigenvectors associated with the stored
  !          eigenvalues;
  !
  !        ierr is set to
  !
  !        ierr is set to
  !          zero       for normal return,
  !          j          if the j-th eigenvalue has not been
  !                     determined after 30 iterations.
  !
  !     ------------------------------------------------------------------
  !
  !      implicit real*8 (a-h,o-z)

  implicit none
  integer i, j, k, l, m, n, ii, mml, ierr
  real*8  machep, c, f, b, s, p, g, r
  real*8 d(n),e(n),z(n)

  !     ========== machep is a machine dependent parameter specifying
  !                the relative precision of floating point arithmetic.
  !                machep = 16.0d0**(-13) for long form arithmetic
  !                on s360 ==========
  machep=1.0e-14
  !
  ierr = 0
  if (n .eq. 1) goto 1001
  !
  e(n) = 0.0d0
  do  l = 1, n
     j = 0
     !     ========== look for small sub-diagonal element ==========
105  do m = l, n
        if(m.eq.n) goto 120
        if(abs(e(m)).le.machep*(abs(d(m))+abs(d(m+1)))) goto 120
     end do
     !
120  p = d(l)
     if (m .eq. l) go to 240
     if (j .eq. 30) go to 1000
     j = j + 1
     !     ========== form shift ==========
     g = (d(l+1) - p) / (2.0d0   * e(l))
     r =  sqrt(g*g+1.0d0  )
     g = d(m) - p + e(l) / (g +  sign(r, g))
     s = 1.0d0
     c = 1.0d0
     p = 0.0d0
     mml = m - l
     !     ========== for i=m-1 step -1 until l do -- ==========
     do ii = 1, mml
        i = m - ii
        f = s * e(i)
        b = c * e(i)
        if ( abs(f) .lt.  abs(g)) go to 150
        c = g / f
        r =  sqrt(c*c+1.0d0  )
        e(i+1) = f * r
        s = 1.0d0   / r
        c = c * s
        go to 160
150     s = f / g
        r =  sqrt(s*s+1.0d0  )
        e(i+1) = g * r
        c = 1.0d0   / r
        s = s * c
160     g = d(i+1) - p
        r = (d(i) - g) * s + 2.0d0   * c * b
        p = s * r
        d(i+1) = g + p
        g = c * r - b
        !     ========== form first component of vector ==========
        f = z(i+1)
        z(i+1) = s * z(i) + c * f
        z(i) = c * z(i) - s * f
     end do
     !
     d(l) = d(l) - p
     e(l) = g
     e(m) = 0.0d0
     go to 105
240 end do
  !     ========== order eigenvalues and eigenvectors ==========
  do ii = 2, n
     i = ii - 1
     k = i
     p = d(i)
     !
     do j = ii, n
        if (d(j) .ge. p) go to 260
        k = j
        p = d(j)
260  end do
     !
     if (k .eq. i) go to 300
     d(k) = d(i)
     d(i) = p
     !
     p = z(i)
     z(i) = z(k)
     z(k) = p
     !
300 end do
  !
  go to 1001
  !     ========== set error -- no convergence to an
  !                eigenvalue after 30 iterations ==========
1000 ierr = l
1001 return
  !     ========== last card of gbtql2 ==========
end subroutine gbtql2
!-----------------------------------!
!+++++++++++++++++++++++++++++++++++!
function  gamfun(z)
  !  this is a procedure that evaluates gamma(z) for
  !     0 lt z le 3 to 16 significant figures
  !    it is based on a chebyshev-type polynomial
  !   approximation given in h. werner and r. collinge, math. comput.
  !    15 (1961), pp. 195-97.
  !   approximations to the gamma function, accurate up to 18 significant
  !   digits, may be found in the paper quoted above
  !
  !
  !
  !      implicit real*8 (a-h,o-z)
  implicit none
  real*8 a(18)
  real*8 p, gamfun, z, t
  integer k1, k

  a(1)=1.0d0
  a(2)=.4227843350984678d0
  a(3)=.4118403304263672d0
  a(4)=.0815769192502609d0
  a(5)=.0742490106800904d0
  a(6)=-.0002669810333484d0
  a(7)=.0111540360240344d0
  a(8)=-.0028525821446197d0
  a(9)=.0021036287024598d0
  a(10)=-.0009184843690991d0
  a(11)=.0004874227944768d0
  a(12)=-.0002347204018919d0
  a(13)=.0001115339519666d0
  a(14)=-.0000478747983834d0
  a(15)=.0000175102727179d0
  a(16)=-.0000049203750904d0
  a(17)=.0000009199156407d0
  a(18)=-.0000000839940496d0
  


  if(z.le.1.0d0  ) go to 10
  if(z.le.2.0d0  ) go to 20
  t=z-2.0d0
  go to 30
10 t=z
  go to 30
20 t=z-1.0d0
30 p=a(18)
  do  k1=1,17
     k=18-k1
     p=t*p+a(k)
40 end do
  
  if(z.gt.2.0d0  ) go to 50
  if(z.gt.1.0d0  ) go to 60
  gamfun=p/(z*(z+1.0d0  ))
  return
60 gamfun=p/z
  return
50 gamfun=p
  return
end function gamfun
!+++++++++++++++++++++++++++++++++++!





