real(8) function probability_of_transition_in_matter_a_b(flvr1,flvr2,L,t12,t23,t13,delta,sm,aM,P,nu,ro,Z,A)
    implicit none
    integer :: flvr1              ! flvr1 is the flavour with which the neutrino is generated
    integer :: flvr2              ! flvr2 is the flavour that is transited 
    real(8) :: L                  ! L is the length between the source of neutrinos an the position [Km]
    real(8) :: t12,t23,t13,delta  ! Are the three mixing angles and the CP-violation phase of the mixing matrix
    real(8) :: sm,aM              ! sm,aM are the squared mass difference m=m_21 y M=m_32
    real(8) :: P                  ! P es el momento del neutrino [GeV]
    integer :: nu                 ! nu is 1 for neutrinos an 2 for antineutrino    
    real(8) :: Ne                 ! Ne is the electron density
    real(8) :: ro                 ! ro [g/cm^3]
    integer :: Z                  ! Z = 1
    integer :: A                  ! A = 2
    real(8) :: probabilityOfTransitionAB
    real(8),parameter :: N_A=6.0221415D23 ! N_A is the Avogadro's constatnt [mo^(-1)]    
    Ne=N_A*ro*REAL(Z)/REAL(A)
    probability_of_transition_in_matter_a_b=probabilityOfTransitionAB(flvr1,flvr2,L,t12,t23,t13,delta,sm,aM,P,nu,Ne)

    return
end function probability_of_transition_in_matter_a_b