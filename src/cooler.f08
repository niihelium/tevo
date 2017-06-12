module dq
  use calculation
  implicit none

  real*8 :: T, Tcmb, ne, nHI, nHII, nHeI, nHeII, nHeIII, z
  type(spectre), dimension(:), allocatable :: spectre_qua
  real*8, dimension(:), allocatable :: energy

contains

  !===========================EXCITATION

  function H_excitation() result(lambda)
    real*8 :: lambda

    lambda = (7.5d-19 * exp(-118348d0/T) * ne * nHI) / (1.0d0 + sqrt(T/1.0d5))
  end function H_excitation

  function He_excitation_11s() result(lambda)
    real*8 :: lambda

    lambda = 1.1d-19 *  T**0.082d0 * exp(-230d3/T) * ne * nHeI
  end function He_excitation_11s

  function He_excitation_23s() result(lambda)
    real*8 :: lambda

    lambda = 9.1d-27 *  T**0.1687d0 * exp(-13179d0/T) * ne**2 * nHeII / (1.0d0 + sqrt(T/1.0d5))
  end function He_excitation_23s

  function HeII_excitation() result(lambda)
    real*8 :: lambda

    lambda = 5.54d-17 *  T**(-0.397d0) * exp(-473638d0/T) * ne * nHeII / (1.0d0 + sqrt(T/1.0d5))
  end function


  !===========================COLLISIONAL IONIZATION
  function H_collisional_ionization() result(lambda)
    real*8 :: lambda

    lambda = 2.179d-11 * H_collisional_ionization_rate(get_T_in_ev(T)) * ne * nHI
  end function

  function He_collisional_ionization() result(lambda)
    real*8 :: lambda

    lambda = 3.94d-11 * He_collisional_ionization_rate(get_T_in_ev(T)) * ne * nHeI
  end function

  function HeII_collisional_ionization() result(lambda)
    real*8 :: lambda

    lambda = 4.95d-22 * sqrt(T) * (1.0d0 + sqrt(T*1d5))**(-1) * exp(-631515/T) * ne * nHeII
  end function

  !===========================RECOMBINATION

  function HII_recombination() result(lambda)
    real*8 :: lambda

    lambda = 1.38d-16 * T * HII_recombination_rate(T) * ne * nHII
  end function

  function HeII_recombination() result(lambda)
    real*8 :: lambda

    lambda = 1.38d-16 * T * HeII_recombination_rr_rate(T) * ne * nHeII
  end function

  function HeII_recombination_dielectronic() result(lambda)
    real*8 :: lambda

    lambda = 6.54d-11 * HeII_recombination_di_rate(T) * ne * nHeII
  end function

  function HeIII_recombination() result(lambda)
    real*8 :: lambda

    lambda = 3.48d-26 * sqrt(T) * ( ((T/1d3)**(-0.2)) /(1 + (T/1d6)**0.7)) * ne * nHeIII
  end function


  function compton_cooling() result(lambda)
    real*8 :: lambda

    lambda = 1.017d-37 * Tcmb**4 * (T - Tcmb) * ne
  end function

  function bremsstrahlung() result(lambda)
    real*8 :: lambda

    lambda = 1.426d-27 * sqrt(T) * &
    (nHII*gff(T, 1) + nHeII*gff(T, 1) + 4*nHeIII*gff(T, 2)) * ne
  end function

  function gff(T, ion) result(gff_result)
    integer, intent(in) :: ion
    real*8, intent(in) :: T
    real*8 ::TiZ2, gff_result

    TiZ2 = T/ion**2

    if (TiZ2 < 320000.0d0 ) then
      gff_result = 0.79464 + 0.1243 * log10(TiZ2)
    elseif ( TiZ2 > 320000.0d0 ) then
      gff_result = 2.13164 - 0.1240 * log10(TiZ2)
    end if
  end function

  !===========================HEATING
  function photoinization(kind) result(sum)
      implicit none
      integer, intent(in) :: kind

      real*8, dimension(:), allocatable :: spectre_z
      integer :: e_size, i
      real*8 :: sum, e_prev, de, e_curr, e1, e2, item_Z, E_min

      e_prev = 0.d0
      sum = 0.d0
      E_min = get_E_min_by(kind)

      do i = size(energy), 1, -1

          ! энергия в eV
          e_curr = energy(i)
          if (e_curr > E_min) then
            de = abs(e_prev - e_curr)

            ! берем спектры для нужного красного смещения
            spectre_z = get_spectre_by_z(z, spectre_qua)

            ! берем значение потока для текущей энергии
            item_Z = spectre_z(i)

            sum = sum + (item_Z* &
                    get_sigma_E(e_curr, kind)* &
                    (ev_to_erg(e_curr) - ev_to_erg(E_min))* &
                    (ev_to_erg(de)/ev_to_erg(e_curr)))



          end if
          e_prev =  e_curr
      end do
      sum = sum * (4.0d0 * 3.14d0)/6.6261e-27
  end function photoinization


  function get_cooling() result(dq)
    real*8 :: dq, type

    dq = H_excitation() + He_excitation_11s() + He_excitation_23s() + HeII_excitation() &
    + H_collisional_ionization() + He_collisional_ionization() + HeII_collisional_ionization() &
    + HII_recombination() + HeII_recombination() + HeII_recombination_dielectronic() + HeIII_recombination() &
    + compton_cooling() + bremsstrahlung()
  end function

  function get_heating() result(dq)
    real*8 :: dq

    dq = photoinization(ion_HI) * nHI &
       + photoinization(ion_HeI) * nHeI &
       + photoinization(ion_HeII)* nHeII
  end function
end module dq
