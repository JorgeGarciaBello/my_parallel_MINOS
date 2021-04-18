subroutine minos_prop_one_slab(flvr1,flvr2,t12,t23,t13,delta,sm,aM,P,nu,result)
    implicit none    
    integer :: flvr1              ! flavor intial 1 2 3 = e, m, t
    integer :: flvr2              ! flavor final 1 2 3 = e, m, t
    real(8) :: t12,t23,t13,delta  ! angle mixing and phase
    real(8) :: sm,aM              ! sm,aM are the squared mass difference m=m_21 y M=m_32
    real(8) :: P                  ! Energy
    integer :: nu                 ! nu is 1 for neutrinos an 2 for antineutrino
    
    real(8) :: L                  ! L is the length of slabs in serie. unitis [km]
    real(8) :: Ne                 ! Ne is the electron density. unitis [mol/cm^{3}]    
    real(8) :: probabilityOfTransitionAB    ! fuction
    real(8) :: result
    
    ! slab longitude
    L= 735           ! [km] far detector minos for de source
    
    ! electron density
    Ne= 1.36d0  ! rho1 Z1/A1 [mol/cm^{3}] acelerator minos
    
    result=0.0d0
    result=probabilityOfTransitionAB(flvr1,flvr2,L,t12,t23,t13,delta,sm,aM,P,nu,Ne)
 return
end subroutine minos_prop_one_slab