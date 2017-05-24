module mod_spectre

  implicit none
  type spectre
    real*8 :: redshift
    real*8, dimension(:), allocatable :: j_spectre
  end type spectre

end module mod_spectre
