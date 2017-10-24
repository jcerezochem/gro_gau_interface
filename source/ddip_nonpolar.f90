program ddip_nonpolar
    
    implicit none

    real(8),dimension(1000) :: q
    real(8) :: qtmp
    integer :: nat
    ! Counters
    integer :: i
    ! Read stuff
    integer :: IOstatus

    ! Read unknown number of lines
    i=0
    do 
        i=i+1
        read(5,*,iostat=IOstatus) qtmp
        if (IOstatus /= 0) exit
        q(i) = qtmp
    enddo
    nat = i-1

    ! Now write ddip
    do i=1,nat
        ! wrt x_i
        write(6,'(3D20.12)') q(i), 0.d0, 0.d0
        ! wrt y_i
        write(6,'(3D20.12)') 0.d0, q(i), 0.d0 
        ! wrt z_i
        write(6,'(3D20.12)') 0.d0, 0.d0, q(i)
    enddo

    stop

end program ddip_nonpolar
