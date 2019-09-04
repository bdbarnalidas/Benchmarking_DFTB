!!* The main dftb+ program
program wrapper
use main
implicit none
integer(kind=8) :: atoms=2
integer(kind=8) :: dynamics=0
double precision :: array1(3,2)
double precision :: array2(3,2)


double precision :: xyz(3,2)
double precision :: forces(3,2)
double precision :: energy
integer :: i,j

!array1=(reshape((/0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,37.794519d0/),shape(array1)))
!array2=(reshape((/0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,2.0d0/),shape(array1)))

!do i=1,3 !100

!dynamics=dynamics+1
!if(i.eq.1) xyz=array1
!if(i.gt.1) xyz=array2
!array1(3,2)=array1(3,2)-0.5
!xyz=array1
call dftb(xyz,atoms,dynamics,forces,energy)

!print*,"THE REPORTED ENERGY ", energy,i,array1(3,2)
!enddo

end program

