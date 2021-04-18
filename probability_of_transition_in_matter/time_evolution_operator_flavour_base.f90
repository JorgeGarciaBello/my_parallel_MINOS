!################################################
!
!   timeEvolutionOperatorFlavourBase: is a subroutine 
!       that load the time evolution operator
!       in the flavour base as function of length L
!
!################################################
subroutine timeEvolutionOperatorFlavourBase(UfL,L,t12,t23,t13,delta,sm,aM,P,nu,Ne)
    implicit none
    double complex :: UfL(3,3)    ! UfL is the time evolution operator matrix in the flavour base
    real(8) :: L                  ! L is the length between the source of neutrinos an the position
    
    real(8) :: t12,t23,t13,delta  ! Are the three mixing angles and the CP-violation phase of the mixing matrix
    real(8) :: sm,aM                ! sm,aM are the squared mass difference m=m_21 y M=m_32
    real(8) :: P                  ! P es el momento del neutrino
    integer :: nu                 ! nu is 1 for neutrinos an 2 for antineutrino    
    real(8) :: Ne                 ! Ne is the electron density

    double complex :: UmL(3,3)           ! UmL is the time evolution operator matrix in the mass base
    double complex :: U(3,3)      ! U is the mixing matrix of the oscillation model
    double complex :: U_1(3,3)    ! U is the mixing matrix of the oscillation model


   call timeEvolutionOperatorMassBase(UmL,L,t12,t23,t13,delta,sm,aM,P,nu,Ne)
   call mixingMatrix(U,t12,t23,t13,delta)
   call inverseMixingMatrix(U_1,t12,t23,t13,delta)
    UfL(:,:) = 0.0d0
   
    UfL = matmul(U,matmul(UmL,U_1))
    return
end subroutine timeEvolutionOperatorFlavourBase



