function iterativeProbabilityOfTransitionAB(dim,flvr1,flvr2,l_a,ro_a,t12,t23,t13,delta,sm,aM,P,nu,Z,A)
    use types
    implicit none
    real(dp) :: iterativeProbabilityOfTransitionAB
    integer  :: dim
    integer  :: flvr1              ! flvr1 is the flavour with which the neutrino is generated
    integer  :: flvr2              ! flvr2 is the flavour that is transited 
    real(dp) :: l_a(dim)            ! L is the length between the source of neutrinos an the position
    real(dp) :: ro_a(dim)           ! is an array with matter density values
    real(dp) :: t12,t23,t13,delta  ! Are the three mixing angles and the CP-violation phase of the mixing matrix
    real(dp) :: sm,aM              ! sm,aM are the squared mass difference m=m_21 y M=m_32
    real(dp) :: P                  ! P es el momento del neutrino
    integer  :: nu                 ! nu is 1 for neutrinos an 2 for antineutrino    
    real(dp) :: Ne                 ! Ne is the electron density
    integer  :: Z                  ! Z
    integer  :: A                  ! A

    double complex :: iterativeProbabilityAmplitude

    iterativeProbabilityOfTransitionAB = 0.0d0
    iterativeProbabilityOfTransitionAB = iterativeProbabilityAmplitude(dim,flvr1,flvr2,l_a,ro_a,t12,t23,t13,delta,sm,aM,P,nu,Z,A)*&
                                CONJG(iterativeProbabilityAmplitude(dim,flvr1,flvr2,l_a,ro_a,t12,t23,t13,delta,sm,aM,P,nu,Z,A))
    return
end function iterativeProbabilityOfTransitionAB