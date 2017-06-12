module constants

  implicit none

  integer, parameter :: ion_HI = 1
  integer, parameter :: ion_HII = 2
  integer, parameter :: ion_HeI = 3
  integer, parameter :: ion_HeII = 4
  integer, parameter :: ion_HeIII = 5

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

  !Cosmological constants
  real*8, parameter :: Omega_lambda = 0.7d0
  real*8, parameter :: Omega_m = 0.3d0
  real*8, parameter :: Omega_b = 0.044d0

  ! km/s/Mpc
  real*8, parameter :: H_0 = 67.80d0 !Used only in constants
  real*8, parameter :: H = H_0 * 1d5 / 3d24
  real*8, parameter :: H_c = H_0 / 100.0d0
  real*8, parameter :: H_pre = 2.169d-18

  real*8, parameter :: k_b = 1.38d-16
  real*8, parameter :: m_h = 1.67d-24

  real*8, parameter :: tiny = 1.0d-15

  ! T охлаждения
  real*8, parameter :: T_cooling = 1.0d-20
  real*8, parameter :: yr_to_s = 3.154e7


  ! gramms
  real*8, parameter :: m_p = 1.6726219d-24

  real*8, parameter :: n_0 = 1.88d-29 * (H_0/100.0d0)**2  * Omega_b / m_p

  ! cm/s
  real*8, parameter :: c = 2.998e+10


contains

  function get_concentration_by_z(z) result(n)
    implicit none
    real*8, intent(in) :: z
    real*8 :: n

    n = n_0 * (1+z)**3
  end function

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

  function get_E_min_by(kind) result(E_0)
    implicit none
    integer, intent(in) :: kind
    real*8 :: E_0

  select case (kind)
    case (ion_HI)
      E_0 = HI_MIN
    case (ion_HeI)
      E_0 = HeI_MIN
    case (ion_HeII)
      E_0 = HeII_MIN
    case default
  end select

  end function

  function get_T_in_ev(T) result(ev)
    implicit none
    real*8, intent(in) :: T
    real*8 :: ev

    ev = T * 8.6173303d-5
  end function get_T_in_ev

  ! function yr_to_s(yr) result(s)
  !   implicit none
  !   real*8, intent(in) :: yr
  !   real*8 :: s
  !
  !   s = yr *
  ! end function

end module constants
