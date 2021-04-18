!################################################
!
!   traceHm: is a function that return the
!       trace of the Hamiltonian-mass matrix
!
!################################################
double complex function traceHm(t12,t23,t13,delta,sm,aM,P,nu,Ne)
    implicit none
    real(8) :: t12,t23,t13,delta  ! Are the three mixing angles and the CP-violation phase of the mixing matrix
    real(8) :: sm,aM                ! sm,aM are the squared mass difference m=m_21 y M=m_32
    real(8) :: P                  ! P es el momento del neutrino
    integer :: nu                 ! nu is 1 for neutrinos an 2 for antineutrino    
    real(8) :: Ne                 ! Ne is the electron density

    double complex :: t_Hm(3,3)          ! t_Hm is the sum of Hamiltonian in vacumm (mass base)  and the potental in the base of mass
    
    traceHm=(0.0D0,0.0d0)

    call totalHamiltonianMassBase(t_Hm,t12,t23,t13,delta,sm,aM,P,nu,Ne)
    traceHm=t_Hm(1,1)+t_Hm(2,2)+t_Hm(3,3)
    return
end function traceHm