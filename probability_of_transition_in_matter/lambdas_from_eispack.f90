subroutine lambdaFromEISPACK(Ls,t12,t23,t13,delta,sm,aM,P,nu,Ne)
    implicit none
    double complex :: Ls(3)       ! Ls is an array with the values of coefficients lambda    
    real(8) :: t12,t23,t13,delta  ! Are the three mixing angles and the CP-violation phase of the mixing matrix
    real(8) :: sm,aM              ! sm,aM are the squared mass difference m=m_21 y M=m_32
    real(8) :: P                  ! P es el momento del neutrino
    integer :: nu                 ! nu is 1 for neutrinos an 2 for antineutrino    
    real(8) :: Ne                 ! Ne is the electron density

    double complex :: T(3,3)      ! T is a traceless matrix from the model

    integer,parameter :: n=3
    real(8) :: ar(n,n), ai(n,n)

    real (8) :: wi(n)
    real (8) :: wr(n)
    logical matz
    real (8) :: zr(n,n)
    real (8) :: zi(n,n)
    integer ( kind = 4 ) ierr
    matz=.FALSE.

    call tMatrix(T,t12,t23,t13,delta,sm,aM,P,nu,Ne)

    ar(1,1)=REAL(T(1,1));   ar(1,2)=REAL(T(1,2));  ar(1,3)=REAL(T(1,3));
    ar(2,1)=REAL(T(2,1));   ar(2,2)=REAL(T(2,2));  ar(2,3)=REAL(T(2,3));
    ar(3,1)=REAL(T(3,1));   ar(3,2)=REAL(T(3,2));  ar(3,3)=REAL(T(3,3));
    
    ai(1,1)=IMAG(T(1,1));   ai(1,2)=IMAG(T(1,2));  ai(1,3)=IMAG(T(1,3));
    ai(2,1)=IMAG(T(2,1));   ai(2,2)=IMAG(T(2,2));  ai(2,3)=IMAG(T(2,3));
    ai(3,1)=IMAG(T(3,1));   ai(3,2)=IMAG(T(3,2));  ai(3,3)=IMAG(T(3,3));
    
    call cg_lr ( n, ar, ai, wr, wi, matz, zr, zi, ierr )   ! gets eigenvalues and eigenvectors of a complex general matrix.
    Ls(1)=cmplx(wr(1), wi(1))
    Ls(2)=cmplx(wr(2), wi(2))
    Ls(3)=cmplx(wr(3), wi(3))

    return
end subroutine lambdaFromEISPACK