real(8) function probabilityOfTransitionAB(flvr1,flvr2,L,t12,t23,t13,delta,sm,aM,P,nu,Ne)
    implicit none
    integer :: flvr1              ! flvr1 is the flavour with which the neutrino is generated
    integer :: flvr2              ! flvr2 is the flavour that is transited 
    real(8) :: L                  ! L is the length between the source of neutrinos an the position
    real(8) :: t12,t23,t13,delta  ! Are the three mixing angles and the CP-violation phase of the mixing matrix
    real(8) :: sm,aM              ! sm,aM are the squared mass difference m=m_21 y M=m_32
    real(8) :: P                  ! P es el momento del neutrino
    integer :: nu                 ! nu is 1 for neutrinos an 2 for antineutrino    
    real(8) :: Ne                 ! Ne is the electron density

    double complex :: probabilityAmplitude

    probabilityOfTransitionAB = 0.0d0
    probabilityOfTransitionAB = probabilityAmplitude(flvr1,flvr2,L,t12,t23,t13,delta,sm,aM,P,nu,Ne)*        &
                                CONJG(probabilityAmplitude(flvr1,flvr2,L,t12,t23,t13,delta,sm,aM,P,nu,Ne))
    return
end function probabilityOfTransitionAB