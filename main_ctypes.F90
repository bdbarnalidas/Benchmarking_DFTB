!!* The main dftb+ program
module main_c
use iso_c_binding, only: c_double, c_int
use main

contains

subroutine dftb_call(xyz,forces,atoms,dynamics) bind(c)
integer(c_int), intent(in) :: atoms
logical,intent(in) :: dynamics
real(c_double),intent(in) :: xyz(3,atoms)
real(c_double),intent(out) :: forces(3,atoms)

call dftb(xyz,forces,atoms,dynamics)
end subroutine dftb_call

end module main_c
