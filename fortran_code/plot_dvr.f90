!----------------------------------!
subroutine plot_dvr(nev, v, nrs, rmax, nout)

  use dvrecs
  implicit none

  integer    i, j, irs, nrs, nev, nout
  real*8     drs, rmax
  complex*16 fval(nev), zrad, DVReval1d, v(nbas, nev)

  ! plot the first nplot radial orbitals in the DVR basis
  drs=(rmax)/dble(nrs)

  ! make the form look like 'orbs.dat' from MCHF code
  write(nout, *) '# ', nrs
  do irs=0, nrs
     zrad=dble(irs)*drs
     do i=1, nev
        fval(i)=DVReval1d(zrad, az, nq, nel, v(:, i), xz, wz)
     end do
     write(nout, *) dble(zrad), (dble(fval(i)), i=1, nev)
  end do

end subroutine plot_dvr
!----------------------------------!
