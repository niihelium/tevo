module constants

  implicit none

  integer, parameter :: ion_HI = 1
  integer, parameter :: ion_HeI = 11
  integer, parameter :: ion_HeII = 12


  real*8, parameter :: h_ev = 4.135667662e-15
  real*8, parameter :: h_erg = 6.6261e-27


  real*8, parameter :: HI_MIN = 13.5984
  real*8, parameter :: HeI_MIN = 24.5874
  real*8, parameter :: HeII_MIN = 54.417760

!cm^2
  real*8, parameter :: SIGMA_0_HI = 5.475e-14
  real*8, parameter :: SIGMA_0_HeI = 9.492e-16
  real*8, parameter :: SIGMA_0_HeII = 1.369e-14

!E_0(eV)
  real*8, parameter :: E_0_HI = 4.298e-1
  real*8, parameter :: E_0_HeI = 1.361e1
  real*8, parameter :: E_0_HeII = 1.720

!P
  real*8, parameter :: P_HI = 2.963
  real*8, parameter :: P_HeI = 3.188
  real*8, parameter :: P_HeII = 2.963

!y_a
  real*8, parameter :: y_a_HI = 32.88
  real*8, parameter :: y_a_HeI = 1.469
  real*8, parameter :: y_a_HeII = 32.88

!y w
  real*8, parameter :: y_w_HI = 0.
  real*8, parameter :: y_w_HeI = 2.039
  real*8, parameter :: y_w_HeII = 0.

!y 0
  real*8, parameter :: y_0_HI = 0.
  real*8, parameter :: y_0_HeI = 0.4434
  real*8, parameter :: y_0_HeII = 0.

!y 1
  real*8, parameter :: y_1_HI = 0.
  real*8, parameter :: y_1_HeI = 2.136
  real*8, parameter :: y_1_HeII = 0.

  real*8, parameter :: E_MAX = 9.6370342E+05

contains
  

  function erg_to_ev(erg) result(ev)
    implicit none
    real*8, intent(in) :: erg
    real*8 :: ev
    ev = erg * 6.242e+11
  end function erg_to_ev

  function ev_to_erg(ev) result(erg)
    implicit none
    real*8, intent(in) :: ev
    real*8 :: erg
    erg = ev * 1.6022e-12
  end function ev_to_erg

end module constants
