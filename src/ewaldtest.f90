PROGRAM EWALDTEST

  USE STRUCTURE
  USE EWALD

  IMPLICIT NONE
  
  real :: t1(3), t2(3), t3(3)
  integer  :: nion

  type (coordinate), allocatable, dimension(:) :: coords
  integer, allocatable, dimension(:) :: valsel

  integer :: i
  real :: sumewald

  ! Read the structural information
  CALL READSTRUCTURE(5,6,t1,t2,t3,coords,nion)

  ! Prepare first valence set and calculate ewald sum
  ALLOCATE( valsel(nion) )
 
  DO i=1,nion
     valsel(i) = 1
  END DO

  call EWALDSUM(6,t1,t2,t3,coords,valsel,sumewald)

END
  
  
