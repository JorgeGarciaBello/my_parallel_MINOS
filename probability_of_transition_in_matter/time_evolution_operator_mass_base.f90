!################################################
!
!   timeEvolutionOperatorMassBase: is a subroutine 
!       that load the time evolution operator
!       in the mass base as function of length L
!
!################################################
subroutine timeEvolutionOperatorMassBase(UmL,L,t12,t23,t13,delta,sm,aM,P,nu,Ne)
    implicit none
    double complex :: UmL(3,3)    ! UmL is the time evolution operator matrix in the mass base
    real(8) :: L                  ! L is the length between the source of neutrinos an the position

    !double complex :: complexPhaseFactor ! is the function that return the value of the complex fase factor phi of the model    
    !double complex :: phi         ! is the complex phase factor value        
    real(8) :: t12,t23,t13,delta  ! Are the three mixing angles and the CP-violation phase of the mixing matrix
    real(8) :: sm,aM              ! sm,aM are the squared mass difference m=m_21 y M=m_32
    real(8) :: P                  ! P es el momento del neutrino
    integer :: nu                 ! nu is 1 for neutrinos an 2 for antineutrino    
    real(8) :: Ne                 ! Ne is the electron density

    double complex :: I(3,3)      ! I is the identity matrix
    double complex :: T(3,3)      ! T is the T-matrix (3,3)
    double complex :: T2(3,3)     ! T2 is the T-matrix squared (3,3)
    double complex :: vectA(3)    ! vectA is the vector of the model

    double complex :: Ls(3)       ! Ls is an array with the values of coefficients lambda

    !print*,'timeEvolutionOperatorMassBase'

    call identityMatrix(I)    
    call tMatrix(T,t12,t23,t13,delta,sm,aM,P,nu,Ne)
    !print*,'T',T(1,:)
    !print*,'t12',t12
    !print*,'t23',t23
    !print*,'t13',t13
    !print*,'delta',delta
    !print*,'sm',sm
    !print*,'aM',aM
    !print*,'P',P
    !print*,'nu',nu
    !print*,'Ne',Ne
    call tMatrix2(T2,t12,t23,t13,delta,sm,aM,P,nu,Ne)    
    call vectorA(vectA,L,t12,t23,t13,delta,sm,aM,P,nu,Ne)
    UmL(:,:) =0.0d0
    UmL = vectA(1)*I - complex(0.0d0,L)*vectA(2)*T - (L**2)*matmul(T,T)*vectA(3)
       
    UmL=UmL!*phi
    !print*,'UmL',UmL(1,:)
    !print*,'UmL',UmL(2,:)
    !print*,'UmL',UmL(3,:)
    return
end subroutine timeEvolutionOperatorMassBase