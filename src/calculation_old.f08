module calculation
    use constants
    use mod_spectre
    use reader

    implicit none

contains

    function get_sigma_E(E, ion_kind) result(sig_E)
        implicit none
        real*8, intent(in) :: E
        integer, intent(in):: ion_kind
        real*8 :: sig_E, sigma_E

        if (ion_kind.eq.ion_HI) then
            sig_E = sigma_Ee(E, SIGMA_0_HI, E_0_HI, y_0_HI, P_HI, y_a_HI, y_0_HI, y_1_HI)
        else if (ion_kind.eq.ion_HeI) then
            sig_E = sigma_Ee(E, SIGMA_0_HeI, E_0_HeI, y_0_HeI, P_HeI, y_a_HeI, y_0_HeI, y_1_HeI)
        else if (ion_kind.eq.ion_HeII) then
            sig_E = sigma_Ee(E, SIGMA_0_HeII, E_0_HeII, y_0_HeII, P_HeII, y_a_HeII, y_0_HeII, y_1_HeII)
        end if
    end function get_sigma_E

    function sigma_Ee(E, sigma_0, E_0, y_w, P, y_a, y_0, y_1) result(sig_E)
        implicit none
        real*8, intent(in):: E, sigma_0, E_0, y_w, P, y_a, y_0, y_1
        real*8 :: x, y
        real*8 :: sig_E

        x = E / E_0 - y_0
        y = sqrt(x**2. + y_1**2.)
        sig_E = sigma_0 * ((x - 1.)**2 + y_w**2.) * &
            (y**(0.5 * P - 5.5)) / &
            (1. + sqrt(y / y_a))**P
    end function sigma_Ee

    function gamma_i(e, J_e, sigma) result(gamma)
        implicit none
        real*8, intent(in):: e, J_e, sigma
        real*8 :: gamma

        gamma = J_e * sigma / e
    end function gamma_i

    function get_J_e(e, spector) result(result)
        implicit none
        real*8, intent(in) :: e
        real*8, dimension(:), allocatable :: spector
        real*8 :: d, result, diff
        integer :: ssize, i

        ssize = size(spector)
        result = -1.
        d = huge(0.0d0)

        do i=1,ssize
            diff = abs(e - spector(i))
            if ( diff < d ) then
                d = diff
                result = spector(i)
            end if
        end do
    end function get_J_e

    function get_spectre_by_z(z, spectres) result(spec)
        implicit none
        real*8, intent(in) ::z
        real*8  :: e, d, diff, zz
        integer :: length, i
        real*8, dimension(:), allocatable :: spec
        type(spectre), dimension(:), allocatable, intent(in) :: spectres
        type(spectre), dimension(:), allocatable :: spespecctres
        type(spectre) :: odin_spec

        length = size(spectres)
        d = huge(0.0d0)
        ! print *, "THIS IS length"
        ! print *, length
        !
        ! print *, "THIS IS Z"
        ! print *, z

        !  print *, "THIS IS spectres(i)%redshift"
        !  print *, spectres
        ! call print_spectres_array(spectres)

        do i = 1, length
            odin_spec = spectres(i)
            zz = odin_spec%redshift
            ! print *, "ZZZZZZZZZZZZZZZZZZZZZZZZZZZzz"
            ! print *, zz
            ! print *, i
            !
            ! read(*,*)

            diff = abs(z - spectres(i)%redshift)
            ! print *, diff
            if ( diff < d ) then
                d = diff

                spec = spectres(i)%j_spectre
            end if
        end do

        ! call print_array(spec)
    end function get_spectre_by_z

    function integrate_gamma(energy,spectres,kind,z) result(sum)
        implicit none
        real*8, dimension(:), allocatable, intent(in) :: energy
        type(spectre), dimension(:), allocatable, intent(in) :: spectres
        integer, intent(in) :: kind
        real*8, intent(in) :: z

        real*8, dimension(:), allocatable :: spectre_z
        integer :: e_size, i
        real*8 :: sum, e_prev, de, e_curr, e1, e2, item_Z, E_0

        e_prev = 0.d0
        sum = 0.d0
        E_0 = get_E_0_by(kind)

        do i = size(energy), 1, -1

            ! энергия в eV
            e_curr = energy(i)
            if (e_curr > E_0) then
              de = abs(e_prev - e_curr)

              ! берем спектры для нужного красного смещения
              spectre_z = get_spectre_by_z(z, spectres)
              ! print *, size(spectre_z)

              ! берем значение потока для текущей энергии
              item_Z = spectre_z(i)
              ! call print_array(spectre_z)
              ! считаем сечение взаимодействия для текущей энергии и иона
              ! print *, get_sigma_E(e_curr, kind)
              ! print *, e_curr
              ! print *, item_Z
              ! print *, de
              ! read (*,* )
  ! 13.6 11-18
  !
              sum = sum + (item_Z* &
                      get_sigma_E(e_curr, kind)* &
                      (ev_to_erg(de)/ev_to_erg(e_curr)))


            end if
            e_prev =  e_curr
        end do
        sum = sum * (4.0d0 * 3.14d0)/6.6261e-27
    end function integrate_gamma

    function sec_to_yr(sec) result(years)
      implicit none
      real*8, intent(in) :: sec
      real*8 :: years

      years = sec/60.0/60.0/24.0/365.0
    end function sec_to_yr

    function diff_time_z(z) result(time)
      implicit none
      real*8, intent(in) :: z
      real*8 :: time, z_plus_one

      z_plus_one = 1.0d0 + z
      time = 1d0/(z_plus_one*H_pre*sqrt(Omega_lambda + Omega_m*(z_plus_one)**3))
    end function diff_time_z

    function integrate_z(start, z, steps) result(result)
      implicit none
      real*8, intent(in) :: start, z, steps
      real*8 :: result
      real*8 :: n, h
      real*8 :: i

      result = 0d0

      i = 0.d0

      h = (z-start)/steps

      do
        result = result + diff_time_z(start + h * (i + 0.5d0))
        i = i + 1d0
        if (i > steps) then
          exit
        end if
      end do

      result = result * h

      result = sec_to_yr(result)/1d9
      ! print *, result
    end function integrate_z

    ! k11, k13, k24, k25rr, k25di,

    function H_collisional_ionization_rate(t_e) result(k11)
      implicit none
      real*8, intent(in) :: t_e
      real*8 :: k11, logT_e

      logT_e = log(t_e)

      k11 = exp(-3.271396786d1 + 1.35365560d1 * logT_e &
                -5.73932875d0 * logT_e**2 + 1.56315498d0 * logT_e**3 &
                -2.87705600d-1 * logT_e**4 + 3.48255977d-2 * logT_e**5 &
                -2.63197617d-3 * logT_e**6 + 1.11954395d-4 * logT_e**7 &
                -2.03914985d-6 * logT_e**8)
    end function H_collisional_ionization_rate

    function HII_recombination_rate(T) result(k13)
      implicit none
      real*8, intent(in) :: T
      real*8 :: k13

      k13 = 1.269d-13 * (315614.0d0/T)**1.503d0 * (1.0d0 + (604625.0d0/T)**0.470d0)**(-1.923d0)
    end

    function He_collisional_ionization_rate(T_e) result(k24)
      implicit none
      real*8, intent(in) :: T_e
      real*8 :: k24, logT_e

      logT_e = log(T_e)

      k24 = exp(-4.409864886d1 +  2.391596563d1 * logT_e &
                -1.07532302d1 * logT_e**2 +  3.05803875d0 * logT_e**3 &
                -5.6851189d-1 * logT_e**4 + 6.79539123d-2 * logT_e**5 &
                -5.0090561d-3 * logT_e**6 + 2.06723616d-4 * logT_e**7 &
                -3.64916141d-6 * logT_e**8)
    end function

    function HeII_collisional_ionization_rate(T) result(rate)
      implicit none
      real*8, intent(in) :: T
      real*8 :: rate

      rate = 5.86d-12 * sqrt(T) * (1.0d0 + sqrt(T/1d5))**(-1) * exp(-631515/T)
    end function

    function HeII_recombination_rr_rate(T) result(k25rr)
      implicit none
      real*8, intent(in) :: T
      real*8 :: k25rr, logT
      logT = log10(T)

      k25rr = 1.0d-11* T**(-0.5d0) * (12.72d0 - 1.615d0 * logT - 0.3162*(logT**2) + 0.0493*(logT**3))
    end function

    function HeII_recombination_di_rate(T) result(k25di)
      implicit none
      real*8, intent(in) :: T
      real*8 :: k25di

      k25di = 1.9d-3 * T**(-1.5d0) * exp(-473421.0d0/T) * (1.0d0 + 0.3d0*exp(-94684/T))
    end function

    function HeIII_recombination_rate(T) result(rate)
      implicit none
      real*8, intent(in) :: T
      real*8 :: rate

      rate = 3.36d-10 * T**(-0.5d0) * (T/1d3)**(-0.2d0) / (1.0d0 + (T/1d6)**0.7)
    end function

    function n_b(z) result(n_b_result)
      implicit none
      real*8, intent(in) :: z
      real*8 :: n_b_result

      n_b_result = Omega_b * 1.88e-29 * H_c**2 * (1.0d0 + z)**3/(0.63*m_h)

    end function n_b

    function get_Tcmb(z) result(Tcmb)
      implicit none
      real*8, intent(in) :: z
      real*8 :: Tcmb

      Tcmb = 2.7d0*(1 + z)
    end function get_Tcmb

    function get_H(z) result(Hz)
      implicit none
      real*8, intent(in) :: z
      real*8 :: Hz

      Hz = H * sqrt(Omega_lambda + Omega_m*(1 + z)**3)
    end function get_H

end module calculation
