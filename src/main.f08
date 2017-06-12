module equations
  use dq
  use constants
  use mod_spectre
  use calculation


  implicit none

  real*8 :: nH_tot, nHe_tot
  real*8 :: z_step
  real*8 :: F_ion_HI, F_ion_HII, F_ion_HeI, F_ion_HeII, F_ion_HeIII
  real*8 :: D_ion_HI, D_ion_HII, D_ion_HeI, D_ion_HeII, D_ion_HeIII
  real*8 :: heating, cooling


  contains

    function D(ion_kind) result(result)
      integer, intent(in) :: ion_kind
      real*8 :: result

      select case (ion_kind)
        case (ion_HI)
          result = H_collisional_ionization_rate(get_T_in_ev(T))*ne &
            + integrate_gamma(energy, spectre_qua, ion_HI, z)
        case (ion_HII)
          result = HII_recombination_rate(T)*ne
        case (ion_HeI)
          result = He_collisional_ionization_rate(get_T_in_ev(T))*ne &
            + integrate_gamma(energy, spectre_qua, ion_HeI, z)
        case (ion_HeII)
          result = HeII_collisional_ionization_rate(T)*ne &
            + integrate_gamma(energy, spectre_qua, ion_HeII, z) &
            + HeII_recombination_rr_rate(T)*ne &
            + HeII_recombination_di_rate(T)*ne
        case (ion_HeIII)
          result = HeIII_recombination_rate(T)*ne
        case default
      end select
    end function

    function F(ion_kind) result(result)
      integer, intent(in) :: ion_kind
      real*8 :: result

      select case (ion_kind)
        case (ion_HI)
          result = HII_recombination_rate(T)*ne*nHII
        case (ion_HII)
          result = H_collisional_ionization_rate(get_T_in_ev(T))*ne*nHI &
            + integrate_gamma(energy, spectre_qua, ion_HI, z)*nHI
        case (ion_HeI)
          result = HeII_recombination_rr_rate(T)*ne*nHeII &
            + HeII_recombination_di_rate(T)*ne*nHeII
        case (ion_HeII)
          result = He_collisional_ionization_rate(get_T_in_ev(T))*ne*nHeI &
            + integrate_gamma(energy, spectre_qua, ion_HeI, z)*nHeI &
            + HeIII_recombination_rate(T)*nHeIII*ne
        case (ion_HeIII)
          result = HeII_collisional_ionization_rate(T)*ne*nHeII &
             + integrate_gamma(energy, spectre_qua, ion_HeII, z)*nHeII
        case default
      end select
    end function

    SUBROUTINE FEX(NEQ, Time, Y, YDOT)
      IMPLICIT NONE
      INTEGER NEQ !6
      DOUBLE PRECISION Time, Y, YDOT
      DIMENSION Y(NEQ), YDOT(NEQ)

      intent(in) NEQ, Time, Y
      intent(out) YDOT

      YDOT(1) = F_ion_HI - D_ion_HI * Y(1) !HI
      YDOT(2) = F_ion_HII - D_ion_HII * Y(2) !HII
      YDOT(3) = F_ion_HeI - D_ion_HeI * Y(3) !HeI
      YDOT(4) = F_ion_HeII - D_ion_HeII * Y(4) !HeII
      YDOT(5) = F_ion_HeIII - D_ion_HeIII * Y(5) !HeIII
      YDOT(6) = (2.0d0/(3*k_b*(nH_tot + nHe_tot + ne)))*(heating - cooling) - 2*get_H(Z)*T!T
      ! write(*,*) "YDOT"
      ! write(*,*) YDOT(1), YDOT(2), YDOT(3), YDOT(4), YDOT(5)
      ! write(*,*) YDOT(6), heating, cooling
      ! read(*,*)
    RETURN
    END SUBROUTINE FEX

    function get_Tau_HII() result(tau)
      implicit none
      real*8 :: tau

      tau = ne*get_sigma_E(HI_MIN,ion_HI)*z_step/ (get_H(z)*(1+z))

    end function get_Tau_HII

end module equations



module utils
  use calculation
  use dq
  implicit none
  contains
    ! subroutine write_gamma
    !   implicit none
    !   integer :: i
    ! open(unit=11, file="gamma_H.out", status="new", action="write")
    !
    ! do i = 1, 49
    !   write(11, "(E12.3E2 E12.3E2)") redshifts(i), integrate_gamma(energy, spectre_qua, ion_HI, redshifts(i))
    ! end do
    !
    ! close (11)
    !
    ! open(unit=12, file="gamma_HeI.out", status="new", action="write")
    !
    ! do i = 1, 49
    !   write(12, "(E12.3E2 E12.3E2)") redshifts(i), integrate_gamma(energy, spectre_qua, ion_HeI, redshifts(i))
    ! end do
    !
    ! close (12)
    !
    ! open(unit=13, file="gamma_HeII.out", status="new", action="write")
    !
    ! do i = 1, 49
    !   write(13, "(E12.3E2 E12.3E2)") redshifts(i), integrate_gamma(energy, spectre_qua, ion_HeII, redshifts(i))
    ! end do
    !
    ! close (13)
    ! end subroutine write_gamma
    !
    subroutine write_sigma(energy, z)
      implicit none
      real*8, dimension(:), allocatable, intent(in) :: energy
      real*8, intent(in) :: z
      integer :: i

      open(unit=12, file="sigma.out", status="new", action="write")
      write(12, "(E12.3E2)") z
      do i = size(energy), 1, -1
        write(12, "(E12.3E2 E12.3E2)") energy(i), get_sigma_E(energy(i), ion_HeII)
      end do

      close (12)

    end subroutine write_sigma

    subroutine write_coll_ionization
      use constants
      implicit none
      real*8 :: T

      T = 0.0d0

      open(unit=12, file="coll_H.out", status="new", action="write")
      open(unit=13, file="coll_HeI.out", status="new", action="write")
      open(unit=14, file="coll_HeII.out", status="new", action="write")

      !10^3 - 10^6
      do while (T < 1.0d6)
        write(12, "(E12.3E2 E12.3E2)") T, H_collisional_ionization_rate(get_T_in_ev(T))
        write(13, "(E12.3E2 E12.3E2)") T, He_collisional_ionization_rate(get_T_in_ev(T))
        write(14, "(E12.3E2 E12.3E2)") T, HeII_collisional_ionization_rate(T)

        T = T+100.0
      end do

      close (12)
      close (13)
      close (14)

    end subroutine write_coll_ionization

    subroutine write_recombination
      use constants
      implicit none
      real*8 :: T

      T = 0.0d0

      open(unit=12, file="rec_HII.out", status="new", action="write")
      open(unit=13, file="rec_HeII.out", status="new", action="write")
      open(unit=14, file="rec_HeIII.out", status="new", action="write")

      !10^3 - 10^6
      do while (T < 1.0d6)
        write(12, "(E12.3E2 E12.3E2)") T, HII_recombination_rate(T)
        write(13, "(E12.3E2 E12.3E2)") T, HeII_recombination_rr_rate(T) + HeII_recombination_di_rate(T)
        write(14, "(E12.3E2 E12.3E2)") T, HeIII_recombination_rate(T)

        T = T+100.0
      end do

      close (12)
      close (13)
      close (14)

    end subroutine write_recombination

    subroutine write_cooling
      use constants
      implicit none

      T = 1.0d0

      open(unit=11, file="cooling.out", status="new", action="write")

      !10^3 - 10^6
      do while (T < 1.0d7)
        write(11, "(14(E12.3,2x))") T, H_excitation(), &
                                         He_excitation_11s(), &
                                         He_excitation_23s(), &
                                         HeII_excitation(), &
                                         H_collisional_ionization(), &
                                         He_collisional_ionization(), &
                                         HeII_collisional_ionization(), &
                                         HII_recombination(), &
                                         HeII_recombination(), &
                                         HeII_recombination_dielectronic(), &
                                         HeIII_recombination(), &
                                         compton_cooling(), &
                                         bremsstrahlung()

        T = T+100.0

      end do

      close (11)

    end subroutine write_cooling
end module utils



program tevo
  use reader
  use mod_spectre
  use calculation
  use constants
  use utils
  use DVODE_F90_M
  use equations

  implicit none
    !=========INPUT
    real*8, dimension(:), allocatable :: redshifts, spectre_test
    !=========VARIABLES
    real*8 :: z_min, Time_start, Time_step, universe_age, z_start, n_coeff
    real*8 :: Tau_HII
    !=========DVODE_F90
    INTEGER, parameter :: NEQ = 6
    DOUBLE PRECISION ATOL, RTOL, Time, TOUT, Y, RSTATS
    INTEGER ITASK, ISTATE, ISTATS, IOUT, IERROR
    INTEGER i, j
    DIMENSION Y(NEQ), ATOL(NEQ), RSTATS(22), ISTATS(31)
    logical :: write_data, calculating

    TYPE(VODE_OPTS) :: OPTIONS

    write_data = .false.
    calculating = .true.

    !=========READ INPUT DATA
    energy = read_file("data/energy.dat")
    redshifts = read_file("data/redshifts.dat")
    spectre_qua =  read_file_with_z("data/hm_qua.dat", redshifts) !массив спектров со смещением

    Tau_HII = 0.0d0

    ! call print_array(spectre_qua(5)%j_spectre)
    ! read(*,*)

    !initial conditions
    ! z_min = redshifts(size(redshifts))

    ! n_tot = get_concentration_by_z(z)

    ! ============= INITIAL CONDITIONS =============
    z_start = 8.0d0
    ! z_start = .5d0
    z_step = 1d-3

    z = z_start

    nH_tot = n_0 * (1+z)**3
    nHI = 1.0d-4 * nH_tot
    nHII =  nH_tot - nHI

    nHe_tot = 0.06d0 * nH_tot
    nHeI = 1.0d-4 * nHe_tot
    nHeII = 1.0d-4 * nHe_tot
    nHeIII = nHe_tot - nHeI - nHeII

    nE = nHII + nHeII + 2*nHeIII

    T = 1d4

    universe_age = integrate_z(0.0d0, 10000.0d0, 1.0d6)


    Tau_HII = Tau_HII + get_Tau_HII()
    ! do i = 1, size(redshifts)
    !   do j = 1, size(energy)
    !     spectre_qua(i)%j_spectre(j) = 1d-18 * (energy(j)/HI_MIN)**(-1)
    !   end do
    ! end do


  ! ============= ============= ============= =============

  ! spectre_test =  get_spectre_by_z(z, spectre_qua)
  ! call print_array(get_spectre_by_z(z, spectre_qua))

  ! write(*,*) integrate_z(1.0d0, 0.0d0, 1.0d6)
  !
  ! write(*,*) get_heating(), get_cooling()

    if (write_data) then
      call write_coll_ionization()
      call write_recombination()
      call write_cooling()

      open(unit=11, file="concentration.out", status="new", action="write")

      do i = 1, size(redshifts)
          write(11, "(E12.3E2 E12.3E2)") redshifts(i), get_concentration_by_z(redshifts(i))
      end do

      close(11)
    end if

    !Initial conditions
    !evolution loop

    open(unit=123, file="output.out", status="unknown", action="write")
    open(unit=444, file="lambda.out", status="unknown", action="write")

    ! =============================== MAIN LOOP ===============================
    do while (z_start > 1.0d-2 .and. calculating)

      Y(1) = nHI
      Y(2) = nHII
      Y(3) = nHeI
      Y(4) = nHeII
      Y(5) = nHeIII
      Y(6) = T

      RTOL = 1.D-4
      ATOL(:) = 1.d-4 * Y(:)
      ITASK = 1
      ISTATE = 1
      IERROR = 0
      OPTIONS = SET_OPTS(ABSERR_VECTOR=ATOL, &
        RELERR=RTOL, MXSTEP=1000000)

      Time_start = universe_age - integrate_z(0.0d0, z, 1d6)
      TOUT = universe_age - integrate_z(0.0d0, z - z_step, 1d6)



      F_ion_HI = F(ion_HI)
      F_ion_HII = F(ion_HII)
      F_ion_HeI = F(ion_HeI)
      F_ion_HeII = F(ion_HeII)
      F_ion_HeIII = F(ion_HeIII)

      D_ion_HI = D(ion_HI)
      D_ion_HII = D(ion_HII)
      D_ion_HeI = D(ion_HeI)
      D_ion_HeII = D(ion_HeII)
      D_ion_HeIII = D(ion_HeIII)

      Tcmb = 2.7d0 * (1 + z)

      heating = get_heating()
      cooling = get_cooling()

      ! write(*,*) heating, cooling
      ! read(*,*)

      CALL DVODE_F90(FEX,NEQ,Y,Time_start,TOUT,ITASK,ISTATE,OPTIONS)

      do i=1, NEQ-1
        if (Y(i) < tiny) then
          Y(i) = tiny
        end if
      end do

      n_coeff = (1 + z - z_step)**3/(1+z)**3

      nHI = Y(1) * n_coeff
      nHII = Y(2) * n_coeff
      nHeI = Y(3) * n_coeff
      nHeII = Y(4) * n_coeff
      nHeIII = Y(5) * n_coeff
      T = Y(6)

      nH_tot = nHI + nHII
      nHe_tot = nHeI + nHeII + nHeIII
      nE = nHII + nHeII + 2*nHeIII

      ! z_step =  2d-3 * k_b * T/get_cooling()


      Tau_HII = Tau_HII + get_Tau_HII()


      WRITE(123,"(10(E15.6,2x))") z, nHI, nHII, nH_tot, nHeI, nHeII, nHeIII, nHe_tot, T, &
                        n_0 * (1 + z)**3, Tau_HII
      write(444, "(16(E12.3,2x))") T, H_excitation(), &
                                      He_excitation_11s(), &
                                      He_excitation_23s(), &
                                      HeII_excitation(), &
                                      H_collisional_ionization(), &
                                      He_collisional_ionization(), &
                                      HeII_collisional_ionization(), &
                                      HII_recombination(), &
                                      HeII_recombination(), &
                                      HeII_recombination_dielectronic(), &
                                      HeIII_recombination(), &
                                      compton_cooling(), &
                                      bremsstrahlung(), &
                                      cooling, &
                                      heating

      z = z - z_step

      if (z < 0) then
        exit
      end if

    END DO

    close(123)
    close(444)

end program tevo
