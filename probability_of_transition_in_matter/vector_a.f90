subroutine vectorA(vectA,L,t12,t23,t13,delta,sm,aM,P,nu,Ne)
    implicit none
    double complex :: vectA(3)    ! vectA is the vector of the model
    real(8) :: L                  ! L is the length between the source of neutrinos an the position
    
    real(8) :: t12,t23,t13,delta  ! Are the three mixing angles and the CP-violation phase of the mixing matrix
    real(8) :: sm,aM              ! sm,aM are the squared mass difference m=m_21 y M=m_32
    real(8) :: P                  ! P es el momento del neutrino
    integer :: nu                 ! nu is 1 for neutrinos an 2 for antineutrino    
    real(8) :: Ne                 ! Ne is the electron density

    double complex :: lambdaMtrx(3,3)! lambdaMtrx is the matrix of the model
    double complex :: vectE(3)    ! vectE is the vector of the model
    double complex :: inverseLambdaMtrx(3,3)! Inverse matrix

    call lambdaMatrix(lambdaMtrx,L,t12,t23,t13,delta,sm,aM,P,nu,Ne) 
    call vectorE(vectE,L,t12,t23,t13,delta,sm,aM,P,nu,Ne)
    call inverseMatrix(lambdaMtrx,inverseLambdaMtrx)

    vectA=matmul(inverseLambdaMtrx,vectE)
 
    return
end subroutine vectorA