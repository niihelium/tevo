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

        x = E / E_0
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
        print *, "THIS IS length"
        print *, length

        print *, "THIS IS Z"
        print *, z

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
        real*8, dimension(:), allocatable :: spectre_z
        integer, intent(in) :: kind
        integer :: e_size, i
        real*8, intent(in) :: z
        real*8 :: sum, e_prev, de, e_curr, e1, e2, item_Z

        e_prev = 0.d0
        sum = 0.d0

        e_size = size(energy)



        do i=1,e_size
            e_curr = energy(i)
            de = abs(e_prev - e_curr)
            spectre_z = get_spectre_by_z(z, spectres)
            ! print *, size(spectre_z)
            item_Z = spectre_z(i)
            ! call print_array(spectre_z)
            sum = sum + (item_Z* &
                    get_sigma_E(e_curr, kind)* &
                    (erg_to_ev(de)/erg_to_ev(e_curr)))
        end do
    end function integrate_gamma
end module calculation
