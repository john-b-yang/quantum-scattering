!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
module dvrecs

  implicit none
  
  integer nq, nel, rel, nbas, nbas2, rbas
  integer lamax
  real*8 theta, R0, Zcharge

  !element boundaries
  real*8, allocatable:: ar(:)
  complex*16, allocatable:: az(:)
  
  !one element DVR
  real*8, allocatable:: xq_l(:), wq_l(:)
  complex*16, allocatable:: xqz(:)
  complex*16, allocatable:: xs(:), ws(:) 
  
  !full grid
  complex*16, allocatable:: xz(:), wz(:), siwz(:)
  
  !2 x kinetic energy and inverse
  complex*16, allocatable:: TXX(:, :, :), TIXX(:, :, :)


end module dvrecs
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
