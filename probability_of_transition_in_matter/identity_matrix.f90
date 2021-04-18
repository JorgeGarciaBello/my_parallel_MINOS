!################################################
!
!   identityMatrix: is a subroutine that load
!       the Identity matrix I(3,3)
!
!################################################
subroutine identityMatrix(I)
    implicit none
    double complex :: I(3,3)             ! I is the identity matrix

    I(1,1)=1.0d0; I(1,2)=0.0d0; I(1,3)=0.0d0;
    I(2,1)=0.0d0; I(2,2)=1.0d0; I(2,3)=0.0d0;
    I(3,1)=0.0d0; I(3,2)=0.0d0; I(3,3)=1.0d0;

    return
end subroutine identityMatrix