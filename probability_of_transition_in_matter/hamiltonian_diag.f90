!################################################
!
!   hamiltonianDiag: is a subroutine that 
!       load the hamiltonian-diag-matrix 
!
!################################################
subroutine hamiltonianDiag(Hdiag)
    implicit none
    double complex :: Hdiag(3,3)  ! Hdiag is the hamiltonian differential 
    real(8) :: E2                 ! E2 is the energy of the mass eigenstate two
    E2=-1.0d0
    
    Hdiag(1,1)=cmplx(E2,0.0d0);     Hdiag(1,2)=0.0d0;             Hdiag(1,3)=0.0d0
    Hdiag(2,1)=0.0d0;               Hdiag(2,2)=cmplx(E2,0.0d0);   Hdiag(2,3)=0.0d0
    Hdiag(3,1)=0.0d0;               Hdiag(3,2)=0.0d0;             Hdiag(3,3)=cmplx(E2,0.0d0)
    return
end subroutine hamiltonianDiag