!################################################
!
!   inverseMixingMatrix: is a subroutine that load
!       the inverse mixing matrix of the model of
!       neutrino oscillation
!
!################################################
subroutine inverseMixingMatrix(U_1,t12,t23,t13,delta)
    implicit none
    double complex :: U_1(3,3)    ! U is the mixing matrix of the oscillation model    
    real(8) :: t12,t23,t13        ! Are the mixing angles of hte oscillation model
    real(8) :: delta              ! delta is the phase factor for CP violations

    real(8) :: s12,s23,s13        ! Are the sin's os the mixing angles of the oscillation model
    real(8) :: c12,c23,c13        ! Are the cos's os the mixing angles of the oscillation model
    double complex :: arg,exp_delta   

    double complex :: U(3,3)      ! U is the mixing matrix of the oscillation model
    
    call mixingMatrix(U,t12,t23,t13,delta) 
    call inverseMatrix(U,U_1)
    
    return
end subroutine inverseMixingMatrix