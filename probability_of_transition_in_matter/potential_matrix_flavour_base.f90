!################################################
!
!   potentialMatrixFlavourBase: is a subroutine 
!       that load the potential Vf matrix 
!
!################################################
subroutine potentialMatrixFlavourBase(Vf,nu,Ne)
    implicit none
    double complex :: Vf(3,3)     ! Vf is the hamiltonian for the propagation of the neutrinos in vacuum
    
    integer :: nu                 ! nu is 1 for neutrinos an 2 for antineutrino    
    real(8) :: Ne                 ! Ne is the electron density    
    real(8) :: A                  ! A is the matter density
    real(8) :: matterDensityUnits ! is a functions that return the matter density value

    A=matterDensityUnits(nu,Ne)   ! Without units

    Vf(1,1)=cmplx(A,0.0d0); Vf(1,2)=0.0d0; Vf(1,3)=0.0d0
    Vf(2,1)=0.0d0;          Vf(2,2)=0.0d0; Vf(2,3)=0.0d0
    Vf(3,1)=0.0d0;          Vf(3,2)=0.0d0; Vf(3,3)=0.0d0

    return
end subroutine potentialMatrixFlavourBase