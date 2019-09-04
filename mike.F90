module main
    use iso_c_binding, only: c_double, c_int

!integer :: atoms
!integer :: dynamics
!double precision :: xyz(3,atoms)
!double precision :: forces(3,atoms)
contains

subroutine dftb(xyz,forces,atoms,dynamics) bind(c,name='dftb')
  implicit DOUBLE PRECISION (a-h,o-z)
  integer(c_int), intent(in), value :: atoms
  integer(c_int), intent(in), value :: dynamics
  real(c_double), intent(in) :: xyz(3,atoms)
  real(c_double), intent(out) :: forces(3,atoms)

  WRITE(6,*) "inital xyz"
  DO i=1,atoms
      WRITE(*,*) i,xyz(1:3,i)
  ENDDO

end subroutine

