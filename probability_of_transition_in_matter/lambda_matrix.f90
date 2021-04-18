subroutine lambdaMatrix(lambdaMtrx,L,t12,t23,t13,delta,sm,aM,P,nu,Ne)
    implicit none
    double complex :: lambdaMtrx(3,3)! lambdaMtrx is the matrix of the model
    real(8) :: L                  ! L is the length between the source of neutrinos an the position

    double complex :: Ls(3)       ! Ls is an array with the values of coefficients lambda    
    real(8) :: t12,t23,t13,delta  ! Are the three mixing angles and the CP-violation phase of the mixing matrix
    real(8) :: sm,aM              ! sm,aM are the squared mass difference m=m_21 y M=m_32
    real(8) :: P                  ! P es el momento del neutrino
    integer :: nu                 ! nu is 1 for neutrinos an 2 for antineutrino    
    real(8) :: Ne                 ! Ne is the electron density
    
    call lambdaFromEISPACK(Ls,t12,t23,t13,delta,sm,aM,P,nu,Ne)

    lambdaMtrx(1,1)=1.0d0;  lambdaMtrx(1,2)=complex(0.d0,-L)*Ls(1);  lambdaMtrx(1,3)=-(L**2)*(Ls(1)**2);
    lambdaMtrx(2,1)=1.0d0;  lambdaMtrx(2,2)=complex(0.d0,-L)*Ls(2);  lambdaMtrx(2,3)=-(L**2)*(Ls(2)**2);
    lambdaMtrx(3,1)=1.0d0;  lambdaMtrx(3,2)=complex(0.d0,-L)*Ls(3);  lambdaMtrx(3,3)=-(L**2)*(Ls(3)**2);

    return
end subroutine lambdaMatrix