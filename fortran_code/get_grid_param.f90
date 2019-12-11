!----------------------------------!
subroutine get_grid_param()

  use numbers
  use dvrecs
  implicit none

  integer i

  read(501, *)
  read(501, *)
  read(501, *) nq, mq !Labatto order, Radau order
  read(501, *)
  read(501, *) nel, rel
  read(501, *)      
  allocate(ar(0:nel), az(0:nel))
  do i=0, nel
     read(501, *) ar(i)
  end do
  read(501, *)
  read(501, *) theta


  R0=ar(rel) 
  eit=exp(eye*(pi*theta/180.d0))
  az(0:rel)=ar(0:rel)
  az(rel+1:nel)=R0+eit*(ar(rel+1:nel)-R0)

  rbas=rel*(nq-1)-1
  if(rel==nel)then
    nbas=nel*(nq-1)-1
  else
    ndvr=(nel-1)*(nq-1)+mq !Generalized for Radua with different order mq
    nbas=ndvr-2 ! Only taking off basis function at the origin
  end if
  nbas2=nbas**2

end subroutine get_grid_param
!----------------------------------!
