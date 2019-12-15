!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
module numbers

  implicit none
  complex*16, parameter:: eye=(0.d0, 1.d0)
  complex*16, parameter:: zeroz=(0.d0, 0.d0)
  complex*16, parameter:: onez=(1.d0, 0.d0)  
  real*8, parameter:: pi=3.1415926535897932385d0
  real*8, parameter:: onethird=0.33333333333333333333333333

  real*8, parameter:: sol=137.0360285534d0
  real*8, parameter:: MBpersqa0=28.0028509d0
  real*8, parameter:: hartperev=1.d0/27.21140d0
  real*8, parameter:: eVperhart=27.211384d0

  integer nfls
  real*8, allocatable:: rfl(:)


end module numbers
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
