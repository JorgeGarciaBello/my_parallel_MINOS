!################################################
!
!   totalHamiltonianMassBase: is a subroutine
!       that return the the total Hamiltonian 
!       which is the sum of hamiltonian in 
!       vasuum and the potential matrix in
!       mass base
!
!################################################
subroutine totalHamiltonianMassBase(t_Hm,t12,t23,t13,delta,sm,aM,P,nu,Ne)
    implicit none
    double complex :: t_Hm(3,3)   ! t_Hm is the sum of Hamiltonian in vacumm (mass base)  and the potental in the base of mass
    real(8) :: t12,t23,t13,delta  ! Are the three mixing angles and the CP-violation phase of the mixing matrix
    real(8) :: sm,aM                ! sm,aM are the squared mass difference m=m_21 y M=m_32
    real(8) :: P                  ! P es el momento del neutrino
    integer :: nu                 ! nu is 1 for neutrinos an 2 for antineutrino    
    real(8) :: Ne                 ! Ne is the electron density

    double complex :: Hv(3,3)            ! Hv is the hamiltonian for the propagation of the neutrinos in vacuum
    double complex :: Vm(3,3)            ! Vm is the potential matrix in mass base

    call hamiltonianInVacuum(Hv,sm,aM,P)
    !print*, ''
    !print*,'Hv',Hv(1,:)
    !print*,'Hv',Hv(2,:)
    !print*,'Hv',Hv(3,:)
    !print*, ''
    call potentialMatrixMassBase(Vm,t12,t23,t13,delta,nu,Ne)    
    !print*, ''
    !print*,'Vm',Vm(1,:)
    !print*,'Vm',Vm(2,:)
    !print*,'Vm',Vm(3,:)
    !print*, ''
    t_Hm=Hv+Vm
    return
end subroutine totalHamiltonianMassBase