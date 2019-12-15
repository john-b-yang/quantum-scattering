!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
module dvrecs

  implicit none
  
  !integer, parameter :: dbl = Selected_Real_Kind(p=30,r=200)
  integer, parameter :: dbl = Selected_Real_Kind(p=30)
  integer,parameter::unit_stdout=6
  integer nq, mq, nel, rel, nbas, nbas2, rbas, ndvr
  integer lamax
  real*8 theta, R0, Zcharge
  complex*16 eit, eit_conj, eit2

  !element boundaries
  real*8, allocatable:: ar(:)
  complex*16, allocatable:: az(:)
  
  !one element DVR
  real*16, allocatable:: xq_l(:), wq_l(:)
  !real(kind=dbl), allocatable::  xq_r(:), wq_r(:)
  real*16, allocatable::  xq_r(:), wq_r(:)
  complex*16, allocatable:: xqz(:)
  complex*16, allocatable:: xs(:), ws(:)
  
  !full grid
  complex*16, allocatable:: xz(:), wz(:), siwz(:)
  
  !2 x kinetic energy and inverse
  complex*16, allocatable:: TXX(:, :, :), TIXX(:, :, :)


end module dvrecs
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
