program tevo
  use reader
  use mod_spectre
  use calculation

  implicit none

  integer :: lines, int
  type(spectre) :: test
  real*8, dimension(:), allocatable :: energy, redshifts, spez
  real*8 :: sigma
  type(spectre), dimension(:), allocatable :: spectre_qua

  energy = read_file("data/energy.dat")
  redshifts = read_file("data/redshifts.dat")
  !
  spectre_qua =  read_file_with_z("data/hm_qua.dat",redshifts) !массив спектров со смещением
  !
  sigma = integrate_gamma(energy, spectre_qua, ion_HI, redshifts(3))
  ! spez = get_spectre_by_z(redshifts(10),spectre_qua)
  !
   write(*,*) sigma
  ! call print_array(spez)

end program tevo
