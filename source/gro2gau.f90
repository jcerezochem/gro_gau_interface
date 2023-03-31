program export_gmx


    ! ======================================================
    ! This program is part of GAUSSIAN -- GROMACS interface
    ! truough External="gmxMM4gauss.sh parameters.dat" keyword
    ! within Gaussian input file
    !=======================================================
    !
    ! Description
    ! -----------
    ! Get energies, dipole moment and forces calculated by GROMACS
    ! and update the corresponding output file to feed Gaussian
    !
    ! Compilation instructions
    ! -------------------------
    !make$ gfortran export_gmx_v5.f90 -o export_gmx_v5.exe

    IMPLICIT NONE

    integer,parameter::MAX_ATOMS=10000 !Below Gaussian09 limits (it is 250 000, they say)
    real(8),parameter:: &
                           BOHRtoANGS= 5.2917720859D-1, &
                           UMAtoKG   = 1.66053873d-27,  &
                           UMAtoAU   = 1.82288839d3,    &
                           AUtoKG    = 9.10938291d-31,  &
                           BOHRtoM   = 5.291772083d-11, &
                           BOHRtoNM  = 5.291772083d-2,  &
                           AMStoM    = 1.d-10,          &
                           ANGStoM   = 1.d-10,          &
                           HARTtoJ   = 4.3597482d-18,   &
                           HtoKCALM  = 627.5095d0,      &
                           CALtoJ    = 4.184,           &
                           HtoeV     = 27.2114,         &
                           autown    = 2.1947463068d5    !From FCclasses Freq from AU to cm-1
    real(8),parameter:: &
                       DEBYE_TO_eBOHR=0.39343     , &
                       MASS_C=12.0110             , &
                       MASS_H=1.00800

    ! Derived conversion factors
    real(8) :: NM_TO_BOHR, KJMOL_TO_HARTREES


    real(8)::energy, t, force, force_max, mass_tot
    integer::natoms, natoms_gmx_vis, natoms_ua, natoms_aa, nghost, at_max, N, aa2ua
    !Contadores
    integer::i, ii, j

!   File names
    character(len=150)::gmx_energy_mu,gmx_forces,gmx_Ddip,gmx_hessian, & 
                        gmx_pol,gauss_out_file,gauss_msg_file,aa2ua_list

!   OUTPUT DATA MATRICES
    real(8),dimension(1:3)::dipole
    real(8),dimension(1:3,1:MAX_ATOMS)::grad, grad_gauss
    real(8),dimension(1:6)::polarizability
    real(8),dimension(1:9,1:MAX_ATOMS)::Ddipole
    real(8),dimension(1:MAX_ATOMS*MAX_ATOMS)::hessian 
    real(8),dimension(MAX_ATOMS,MAX_ATOMS) :: hess_pp

    ! Units
    integer:: I_ener=10,   & !Input with the energies and dipoles
              I_forc=11,   & !Input with the gradients (forces)
              I_pol=12,    & !Input with the polarizabilities and dipole derivatives
              I_ddip=15,   & !Input with the polarizabilities and dipole derivatives
              I_hess=13,   & !Input with the hessian
              O_gout=21 ,  & !Output for gaussian
              O_msg=22,    & !Output with the messages
              S_hess=30,   &
              I_ua=14        !Input list for ua<->aa conversion
    ! Stuff to read from command line
    integer :: iargc

    ! Input data from argument list
    call getarg(1,gmx_energy_mu)
    call getarg(2,gmx_forces)
    call getarg(3,gmx_pol)
    call getarg(4,gmx_Ddip)
    call getarg(5,gmx_hessian)
    call getarg(6,gauss_out_file)
    call getarg(7,gauss_msg_file)
    call getarg(8, aa2ua_list) 

    ! Open files
    open(I_ener,file=adjustl(gmx_energy_mu),status='old')
    if ( adjustl(gmx_forces) /= 'NO' ) then
        open(I_forc,file=adjustl(gmx_forces),status='old')
    endif
    if ( gmx_pol /= 'ZERO' ) then
        open(I_pol,file=adjustl(gmx_pol),status='old')
    endif
    if ( gmx_Ddip /= 'ZERO' ) then
        open(I_ddip,file=adjustl(gmx_Ddip),status='old')
    endif
    if ( gmx_hessian /= 'NO' ) then
        open(I_hess,file=adjustl(gmx_hessian),status='old')
    endif
    open(O_gout,file=adjustl(gauss_out_file),status='new')
    open(O_msg,file=adjustl(gauss_msg_file),status='old',position='append')
    if ( aa2ua_list /= 'NO' ) then
        open(I_ua,file=adjustl(aa2ua_list),status='old')
    endif

    ! READ DATA
    !1. Energy and dipole moment
    read(I_ener,*) t, energy, dipole(1:3)
    read(I_ener,*) natoms         !This is g09 natoms
    read(I_ener,*) natoms_gmx_vis !This is gmx natoms 
                                  ! as given by update_trr_vX

!    !2. Forces
!    if ( adjustl(gmx_forces) /= 'NO' ) then
!        read(I_forc,*) t, grad(1:3,1:natoms)
!        close(I_forc)
!    endif
    !2. Forces
    if ( adjustl(gmx_forces) /= 'NO' ) then
        do i=1,natoms_gmx_vis
            read(I_forc,'(14X,3(2X,G12.5))') grad(1:3,i)
        enddo
        close(I_forc)
    endif


    if ( aa2ua_list /= 'NO' ) then
!==============================================
!   Transform UA to AA
    read(I_ua,*) natoms_aa, natoms_ua, nghost
    ii=1
    do i=1,natoms_ua+nghost
        read(I_ua,*) aa2ua
        if ( aa2ua /= 0 ) then
            mass_tot = MASS_C + (aa2ua-1)*MASS_H
            grad_gauss(1:3,ii) = grad(1:3,i)*MASS_C/mass_tot
            do j=1,aa2ua-1
                grad_gauss(1:3,ii+j) = grad(1:3,i)*MASS_H/mass_tot
            enddo
            ii = ii + aa2ua
        else
            grad_gauss(1:3,ii) = 0.
            ii = ii + 1
        endif
    enddo
    close(I_ua)

    do i=natoms_ua+1,natoms_gmx_vis
        grad_gauss(1:3,ii) = grad(1:3,i)
        ii = ii + 1
    enddo
!    natoms = natoms + natoms_aa - natoms_ua + ighost
!===============================================
    else
        grad_gauss=grad
    endif


    !3. Polizabilities and dipole derivatives
    if ( adjustl(gmx_pol) == 'ZERO' ) then
        polarizability(1:6)=0.d0
!   else
!        read(I_pol,*) !Don't know the file/format yet
!        close(I_pol)
    endif
    
    if ( adjustl(gmx_Ddip) == 'ZERO' ) then
        Ddipole(1:9,1:natoms)=0.d0
    else
        do i=1,natoms
            read(I_ddip,'(3D20.12)') Ddipole(1:3,i)
            read(I_ddip,'(3D20.12)') Ddipole(4:6,i)
            read(I_ddip,'(3D20.12)') Ddipole(7:9,i)
        enddo
        close(I_ddip)
    endif

    !4. Hessian 
    if ( gmx_hessian /= 'NO' ) then
        read(I_hess,*) N
!! First preprocess
!do i=1,N
!    read(I_hess,*) hess_pp(i,1:N)
!    do j=1,N
!        if (i==13 .or. i==14 .or. i==12 .or. i==9 .or. &
!            j==13 .or. j==14 .or. j==12 .or. j==9 ) hess_pp(i,j) = 0.d0
!    enddo
!enddo
!! Now rewrite
!open(S_hess,file="hess_pp.dat",status='replace')
!do i=1,N
!write(S_hess,*) hess_pp(i,1:N)
!enddo
!rewind(S_hess)
        j=1
        do i=1,N
            read(I_hess,*) hessian(j:j+i-1)
            j=j+i
        enddo
!    close(S_hess)
    close(I_hess)
    endif


    ! CONVERSION FACTORS
    KJMOL_TO_HARTREES = 1.d0/CALtoJ/HtoKCALM
    NM_TO_BOHR = 1.d0/BOHRtoNM
    ! UNIT CONVERSION (and sign for grad)           ! GROMACS       --> Atomic Units
    energy=energy*KJMOL_TO_HARTREES                 ! KJ/mol        --> Hartree
    dipole=dipole*DEBYE_TO_eBOHR                    ! Debye         --> e * bohr
    grad_gauss=-grad_gauss*KJMOL_TO_HARTREES/NM_TO_BOHR         ! KJ/mol * nm-1 --> Hartree * bohr-1
!   polarizbility                                   ! NOT AVAILABLE
!   Ddipole=Ddipole                                 ! e             --> e
    hessian=hessian*KJMOL_TO_HARTREES/NM_TO_BOHR**2 ! KJ/mol * nm-2 --> Hartree * bohr-2


    ! WRITE DATA
    !1. Energy and dipole moment
    write(O_gout,'(4D20.12)') energy, dipole(1:3)

    force_max=0.d0
    !2. Forces
    if ( adjustl(gmx_forces) /= 'NO' ) then
        do i=1,natoms
!if (i==13 .or. i==14 .or. i==12 .or. i==9) grad_gauss(1:3,i)=0.d0
            write(O_gout,'(3D20.12)') grad_gauss(1:3,i)
            force=dsqrt(grad_gauss(1,i)**2+grad_gauss(2,i)**2+grad_gauss(3,i)**2)
!write(O_msg,'(A,1X,I3, F10.5)') "gmx Force on atom", i, force
            force_max=max(force,force_max)
            if ( force_max == force ) at_max=i
        enddo
        write(O_msg,'(A16,1X,F10.5,A,G10.5,A,1X,I5,/)')   &
              'gmxMM  Max force', force_max, " au (", &
              force_max*NM_TO_BOHR/KJMOL_TO_HARTREES, &
              'KJ/mol * nm-1) on atom', at_max
    endif

    !3. Polizabilities and dipole derivatives
    write(O_gout,'(3D20.12)') polarizability(1:3)
    write(O_gout,'(3D20.12)') polarizability(4:6)
    do i=1,natoms
        write(O_gout,'(3D20.12)') Ddipole(1:3,i)
        write(O_gout,'(3D20.12)') Ddipole(4:6,i)
        write(O_gout,'(3D20.12)') Ddipole(7:9,i)
    enddo        

    !4. Hessian
    if ( gmx_hessian /= 'NO' ) then
        write(O_msg,'(/,A)') 'gmxMM  Writting the Hessian'
        write(O_msg,'(A,2(G10.3),/)') 'gmxMM  First and second elements: ', hessian(1), hessian(2)
!    *****************************************************************************
!    GAUSSIAN HELP (http://www.gaussian.com/g_tech/g_ur/k_external.htm)
!    Pseudocode: FFX(I), I=1,(3*NAtoms*(3*NAtoms+1))/2
!    Line format: 3D20.12
!    Note: The Hessian is given in lower triangular form: Alpha_ij, i=1 to N, j=1 to i
!    *****************************************************************************
        do i=1,3*N*(3*N+1)/2,3
            write(O_gout,'(3D20.12)') hessian(i:i+2)
        enddo
    endif

    close(O_gout)


    !WRITE MESSAGE
    write(O_msg,'(A,/)') 'gmxMM  The output file for Gaussian has been updated'

    close(O_msg)

    stop

end program export_gmx


