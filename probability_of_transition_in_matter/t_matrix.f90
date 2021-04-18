!################################################
!
!   tMatrix: is a subroutine that load
!       the T-matrix (3,3)
!
!################################################
subroutine tMatrix(T,t12,t23,t13,delta,sm,aM,P,nu,Ne)
    implicit none
    double complex :: T(3,3)      ! T is a traceless matrix from the model
    real(8) :: t12,t23,t13,delta  ! Are the three mixing angles and the CP-violation phase of the mixing matrix
    real(8) :: sm,aM                ! sm,aM are the squared mass difference m=m_21 y M=m_32
    real(8) :: P                  ! P es el momento del neutrino
    integer :: nu                 ! nu is 1 for neutrinos an 2 for antineutrino    
    real(8) :: Ne                 ! Ne is the electron density

    double complex :: t_Hm(3,3)   ! t_Hm is the Hamiltonian in the mass base
    double complex :: I(3,3)             ! I is the identity matrix
    double complex :: traceHm,trHm! traceHm is the trace of the Hamiltonian-mass

    
    call totalHamiltonianMassBase(t_Hm,t12,t23,t13,delta,sm,aM,P,nu,Ne)
    
    call identityMatrix(I)    

    trHm=traceHm(t12,t23,t13,delta,sm,aM,P,nu,Ne)

    T=t_Hm-((trHm/3.0d0)*I)
    return
end subroutine tMatrix