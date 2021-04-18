subroutine iterativeTimeEvolutionOperator(num_layer,iterUfL,P,l_a,ro_a,t12,t23,t13,delta,sm,aM,nu,Z,A)
    use types
    use constants
    implicit none
    integer ::  num_layer
    double complex :: iterUfL(3,3)! iterUfL is the time evolution operator matrix in the flavour base iterative for diferents lengths and electron densities        
    real(dp) :: P                  ! P es el momento del neutrino, total energy
    real(dp) :: l_a(num_layer)           ! ls an array with lengths 
    real(dp) :: ro_a(num_layer)          ! is an array with matter density values
    real(dp) :: t12,t23,t13,delta  ! Are the three mixing angles and the CP-violation phase of the mixing matrix
    real(dp) :: sm,aM              ! sm,aM are the squared mass difference m=m_21 y M=m_32    
    integer  :: nu                 ! nu is 1 for neutrinos an 2 for antineutrino
    integer  :: Z                   ! Z
    integer  :: A                   ! A
    real(dp) :: Ne
    double complex :: UfL(3,3)    ! UfL is the time evolution operator matrix in the flavour base
    double complex :: I(3,3)      ! I is the identity matrix
    real(dp) :: ld    
    integer  :: k

    call identityMatrix(I)      
    iterUfL=I
    do k=1,num_layer        
        ld=l_a(k)
        Ne=N_A*ro_a(k)*REAL(Z)/REAL(A)
        call timeEvolutionOperatorFlavourBase(UfL,ld,t12,t23,t13,delta,sm,aM,P,nu,Ne)
        iterUfL=matmul(UfL,iterUfL)
    enddo    
    return
end subroutine iterativeTimeEvolutionOperator