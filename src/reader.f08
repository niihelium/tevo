module reader
  use mod_spectre

  implicit none

contains
  function count_lines(filename) result(lines)
    implicit none
    character(len=*) :: filename
    real*8 :: dummy
    integer :: lines, reason
    !===========================================

    lines = 0
    open(unit=1, file=filename, status="old", action="read")
    do
      read(1,*,IOSTAT=reason)  dummy
      if (reason > 0)  then
        ! write(*,*) "Something wrong"
        exit
      else if (reason < 0) then
        ! write(*,*) "End of file"
        exit
      else
        lines = lines + 1
      end if
    end do
    close(1)
  end function count_lines

  function read_file(filename) result(data)
    implicit none
    character(len=*) :: filename
    real*8, dimension(:), allocatable :: data
    integer lines, reason
    !===========================================


    lines = count_lines(filename)
    allocate(data(lines))
    ! write(*,*) "Lines:",lines


    lines = 1
    open(unit=12, file=filename, status="old", action="read")

    do
      read(12,*,IOSTAT=reason) data(lines)
      if (reason < 0) then
        ! write(*,*) "End of file"
        exit
      end if
      lines = lines + 1
    end do
    close(12)

  end function read_file

  function read_file_with_z(filename, redshifts) result(data)
    implicit none
    character(len=*) :: filename
    type(spectre), dimension(:), allocatable :: data
    type(spectre) :: temp_spectre
    real*8, dimension(:), allocatable :: redshifts, temp, temp_e
    real*8, dimension(:,:), allocatable :: energy
    ! real*8 :: zz
    integer lines, reason, i, redshifts_size

    !===========================================

    lines = count_lines(filename)
    ! write(*,*) "lines", lines
    ! write(*,*) "size(redshifts)", size(redshifts)
    redshifts_size = size(redshifts)

    ! write(*,*) "redshifts_size", redshifts_size

    allocate(energy(redshifts_size+1, lines))
    allocate(data(redshifts_size))


    open(unit=13, file=filename, status="old", action="read")
    read(13,*,IOSTAT=reason) energy

    energy = energy(2:,:)




    do i = 1, redshifts_size
      temp_e = energy(i,:)

    ! call print_array(temp_e)
    ! write (*,"(ES15.3E2)") zz
    ! read (*,*)
      temp_spectre = spectre(redshifts(i), temp_e)
      ! call print_array(temp_spectre%j_spectre)
      ! write(*,"(I3 E12.3E2)") 3333333, temp_spectre%redshift
      ! read (*,*)
      data(i) = temp_spectre
    end do
    ! data(lines) = spectre(redshifts(lines), energy)
    close(13)
  end function read_file_with_z

  subroutine print_array(array)
    implicit none
    real*8, dimension(:), allocatable :: array
    integer length, i
    !===========================================

    length = size(array)

    do i=1,length
      write(*,"(I3 ES15.3E2)") i, array(i)
    end do

  end subroutine print_array

  subroutine print_spectres_array(array)
    implicit none
    type(spectre), dimension(:), allocatable :: array
    integer length, i
    !===========================================

    length = size(array)

    do i=1,length
      write(*,"(I3 E12.3E2)") i, array(i)%redshift
      call print_array(array(i)%j_spectre)
    end do

  end subroutine print_spectres_array

end module reader
