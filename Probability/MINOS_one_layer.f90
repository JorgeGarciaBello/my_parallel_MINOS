subroutine MINOS_one_layer(flvr1,flvr2,t12,t23,t13,delta,sm,aM,P,nu,result)
    use types
    use constants
    implicit none  
    real(dp) :: t12,t23,t13,delta  ! angle mixing and phase
    real(dp) :: sm,aM              ! sm,aM are the squared mass difference m=m_21 y M=m_32
    real(dp) :: P                  ! Energy
    integer  :: nu                 ! nu is 1 for neutrinos an 2 for antineutrino
    integer  :: flvr1              ! flavor intial 1 2 3 = e, m, t
    integer  :: flvr2              ! flavor final 1 2 3 = e, m, t
    real(dp) :: L    ! L is the length of slabs in serie. unitis [km]
    real(dp) :: Ne   ! Ne is the electron density. unitis [mol/cm^{3}]

    real(dp) :: probabilityOfTransitionAB    ! fuction
    integer  :: k                  ! counter
    real(dp) :: result

    L= 735.0_dp       ! [km] far detector minos for de source
    ! electron density
    Ne= 1.36d0  ! rho1 Z1/A1 [mol/cm^{3}] acelerator minos
    result=0.0d0
    result=probabilityOfTransitionAB(flvr1,flvr2,L,t12,t23,t13,delta,sm,aM,P,nu,Ne)
    return
end subroutine MINOS_one_layer
