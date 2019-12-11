program one_e_dvr

  use numbers
  use angular, only: lmax
  use dvrecs
  implicit none

  integer i, j, k, lval, nout
  real*8 rplot
  complex*16, allocatable :: AA(:, :), eig(:), ev(:, :)

  open(501, file='indvr.inp', form='formatted', status='old')
  read(501, *)
  read(501, *) lmax
  read(501, *)
  read(501, *) Zcharge

  ! set up radial DVR parameters and one-electron KE matrix
  call setupfls(90)
  call get_grid_param()
  call build_grids()
  lamax=2*lmax
  call build_KEmats()



  nout=2 ! number of bound orbitals to output
  allocate(AA(nbas, nbas), eig(nbas), ev(nbas, nout))
  AA=0.d0
  
  lval=0
! make the one-body hydrogenic Hamiltonian
  AA=0.5d0*TXX(:, :, lval)
  do i=1, nbas
     AA(i, i) = AA(i, i) - Zcharge/xz(i)
  end do

! diagonalize the one-body Ham
  call cdiag(nbas, AA, eig)
  write(6, *) 'writing eigenvectors for these eigenvalues '
  do i=1, nout
     write(6, *) i, dble(eig(i))
  end do

  do i=1, nbas
     write(555, *) dble(xz(i)), dble(AA(i, 1))
  end do

! output the orbitals along a radial portion up to some rmax
  ev=AA(:, 1:nout)
  rplot=25.d0
  call plot_dvr(nout, ev, 150, rplot, 705)
  



end program one_e_dvr

