real(8) function pm1slab(flvr1,flvr2,L,t12,t23,t13,delta,sm,aM,P,nu,rho,za)
 implicit none
 integer :: flvr1,flvr2                   ! flavor initial, flavor final
 real(8) :: L                             ! L is the length of slab [km] [m]
 real(8) :: t12,t23,t13, delta            ! angles of mixing and CP fase
 real(8) :: sm,aM                         ! sm,aM are the squared mass difference m=m_21 y M=m_32 [eV^2]
 real(8) :: P                             ! Energy [GeV] [MeV]
 integer :: nu                            ! nu is 1 for neutrinos an 2 for antineutrino
 real(8) :: rho                           ! density [g/cm^3]
 real(8) :: za                            ! 0<=Z/A<=1 [mol/g] numero de protones entre numero masico

 double complex :: probabilityAmplitude   ! function
 real(8) :: Ne                            ! semy electron density
    
 if ( (za .LT. 0) .OR. (za .GT. 1) ) then
  write (*,*) 'Z/A no puede ser menor a cero ni mayor a uno', za
  stop
 end if

 if ( (rho .LT. 0.0d0)  ) then
  write (*,*) 'no hay densidad negativa', rho
  stop
 endif


 Ne = rho * za

 pm1slab = 0.0d0

 pm1slab =  probabilityAmplitude(flvr1,flvr2,L,t12,t23,t13,delta,sm,aM,P,nu,Ne) *             &
     CONJG( probabilityAmplitude(flvr1,flvr2,L,t12,t23,t13,delta,sm,aM,P,nu,Ne) )


 return
end function pm1slab



double complex function probabilityAmplitude(flvr1,flvr2,L,t12,t23,t13,delta,sm,aM,P,nu,Ne)
    implicit none
    integer :: flvr1              ! flvr1 is the flavour with which the neutrino is generated
    integer :: flvr2              ! flvr2 is the flavour that is transited 
    real(8) :: L                  ! L is the length between the source of neutrinos an the position
    real(8) :: t12,t23,t13,delta  ! Are the three mixing angles and the CP-violation phase of the mixing matrix
    real(8) :: sm,aM                ! sm,aM are the squared mass difference m=m_21 y M=m_32
    real(8) :: P                  ! P es el momento del neutrino
    integer :: nu                 ! nu is 1 for neutrinos an 2 for antineutrino    
    real(8) :: Ne                 ! Ne is the electron density

    double complex :: UfL(3,3)    ! UfL is the time evolution operator matrix in the flavour base
    real(8) :: fState1(3)         ! is the initial flavour eigenstate 
    real(8) :: fState2(3)         ! ! is the final flavour eigenstate

    probabilityAmplitude=0.0d0

    select case (flvr1)
        case(1)
            fState1=(/1.0d0,0.0d0,0.0d0/)
        case(2)
            fState1=(/0.0d0,1.0d0,0.0d0/)
        case(3)
            fState1=(/0.0d0,0.0d0,1.0d0/)
        case DEFAULT
            fState1=(/0.0d0,0.0d0,0.0d0/)
            print*, 'No seleccionaste un eigenestado de sabor inicila definido'
    end select

    select case (flvr2)
        case(1)
            fState2=(/1.0d0,0.0d0,0.0d0/)
        case(2)
            fState2=(/0.0d0,1.0d0,0.0d0/)
        case(3)
            fState2=(/0.0d0,0.0d0,1.0d0/)
        case DEFAULT
            
            print*, 'No seleccionaste un eigenestado de sabor final definido'
    end select

    call timeEvolutionOperatorFlavourBase(UfL,L,t12,t23,t13,delta,sm,aM,P,nu,Ne)

    
    probabilityAmplitude = dot_product(fState1,matmul(UfL,fState2))

    return
end function probabilityAmplitude



!################################################
!
!   timeEvolutionOperatorFlavourBase: is a subroutine 
!       that load the time evolution operator
!       in the flavour base as function of length L
!
!################################################
subroutine timeEvolutionOperatorFlavourBase(UfL,L,t12,t23,t13,delta,sm,aM,P,nu,Ne)
 implicit none
 double complex :: UfL(3,3)    ! UfL is the time evolution operator matrix in the flavour base
 real(8) :: L                  ! L is the length between the source of neutrinos an the position
    
 real(8) :: t12,t23,t13,delta  ! Are the three mixing angles and the CP-violation phase of the mixing matrix
 real(8) :: sm,aM              ! sm,aM are the squared mass difference m=m_21 y M=m_32
 real(8) :: P                  ! P es el momento del neutrino
 integer :: nu                 ! nu is 1 for neutrinos an 2 for antineutrino    
 real(8) :: Ne                 ! Ne is the electron density

 double complex :: UmL(3,3)    ! UmL is the time evolution operator matrix in the mass base
 double complex :: U(3,3)      ! U is the mixing matrix of the oscillation model
 double complex :: U_1(3,3)    ! U is the mixing matrix of the oscillation model

 call timeEvolutionOperatorMassBase(UmL,L,t12,t23,t13,delta,sm,aM,P,nu,Ne)
 call mixingMatrix(U,t12,t23,t13,delta)
 call inverseMixingMatrix(U_1,t12,t23,t13,delta)

 UfL(:,:) = 0.0d0
   
 UfL = matmul(U,matmul(UmL,U_1))
 return
end subroutine timeEvolutionOperatorFlavourBase



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
    !double complex :: phi        ! is the complex phase factor value        
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
    UmL = vectA(1)*I - cmplx(0.0d0,L)*vectA(2)*T - (L**2)*matmul(T,T)*vectA(3)
       
    UmL=UmL!*phi
    !print*,'UmL',UmL(1,:)
    !print*,'UmL',UmL(2,:)
    !print*,'UmL',UmL(3,:)
    return
end subroutine timeEvolutionOperatorMassBase


!################################################
!
!   mixingMatrix: is a subroutine that load
!       the mixing matrix of the model of
!       neutrino oscillation
!
!################################################
subroutine mixingMatrix(U,t12,t23,t13,delta)
    implicit none
    double complex :: U(3,3)      ! U is the mixing matrix of the oscillation model
    real(8) :: t12,t23,t13        ! Are the mixing angles of hte oscillation model
    real(8) :: delta              ! delta is the phase factor for CP violations

    real(8) :: s12,s23,s13        ! Are the sin's os the mixing angles of the oscillation model
    real(8) :: c12,c23,c13        ! Are the cos's os the mixing angles of the oscillation model
    double complex :: arg,exp_delta   

    s12=sin(t12); s23=sin(t23); s13=sin(t13)
    c12=cos(t12); c23=cos(t23); c13=cos(t13)
    
    exp_delta=exp(arg)

    U(1,1)=cmplx(c13*c12,0.0d0);  U(1,2)=cmplx(c13*s12,0.0d0); U(1,3)=s13*exp(cmplx(0.0d0,-delta))    
    U(2,1)=-s12*c23-c12*s23*s13*exp(cmplx(0.0d0,delta)); U(2,2)=c12*c23-s12*s23*s13*exp(cmplx(0.0d0,delta));  U(2,3)=s23*c13
    U(3,1)=s12*s23-c12*c23*s13*exp(cmplx(0.0d0,delta));  U(3,2)=-c12*s23-s12*c23*s13*exp(cmplx(0.0d0,delta)); U(3,3)=c23*c13
        
    return
end subroutine mixingMatrix


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



!################################################
!
!   inverseMatrix: is a subroutine that load
!
!################################################



subroutine inverseMatrix(A,B)
    !! Performs a direct calculation of the inverse of a 3×3 matrix.
    double complex :: A(3,3)    ! Matrix to invert
    double complex :: B(3,3)    ! Inverse matrix
    double complex :: detinv    ! inverse determinant
    double complex :: numerator ! Numerator of the inverse determinant

!    print*, 'A = ', A
    ! Calculate the inverse determinant of the matrix
    numerator=(A(1,1)*A(2,2)*A(3,3) - A(1,1)*A(2,3)*A(3,2) &
             - A(1,2)*A(2,1)*A(3,3) + A(1,2)*A(2,3)*A(3,1) &
             + A(1,3)*A(2,1)*A(3,2) - A(1,3)*A(2,2)*A(3,1))

!    print*, 'numerator = ',  numerator
 
    if (numerator.eq.0.0d0) then 
        print*, 'EL  det|A| = 0 , A no tiene inversa'
        stop
    end if

    detinv = 1.0d0/numerator

    ! Calculate the inverse of the matrix
    B(1,1) = +detinv * (A(2,2)*A(3,3) - A(2,3)*A(3,2))
    B(2,1) = -detinv * (A(2,1)*A(3,3) - A(2,3)*A(3,1))
    B(3,1) = +detinv * (A(2,1)*A(3,2) - A(2,2)*A(3,1))
    B(1,2) = -detinv * (A(1,2)*A(3,3) - A(1,3)*A(3,2))
    B(2,2) = +detinv * (A(1,1)*A(3,3) - A(1,3)*A(3,1))
    B(3,2) = -detinv * (A(1,1)*A(3,2) - A(1,2)*A(3,1))
    B(1,3) = +detinv * (A(1,2)*A(2,3) - A(1,3)*A(2,2))
    B(2,3) = -detinv * (A(1,1)*A(2,3) - A(1,3)*A(2,1))
    B(3,3) = +detinv * (A(1,1)*A(2,2) - A(1,2)*A(2,1))

    ! write(*,*) B  !bug
    return
end subroutine inverseMatrix


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


!################################################
!
!   tMatrix: is a subroutine that load
!       the T-matrix (3,3)
!
!################################################
subroutine tMatrix(T,t12,t23,t13,delta,sm,aM,P,nu,Ne)
    implicit none
    double complex :: T(3,3)      ! T is a traceless matrix from the model
    real(8) :: t12,t23,t13,delta  ! Are the three mixing angles and the CP-violation phase of the mixing matrix
    real(8) :: sm,aM                ! sm,aM are the squared mass difference m=m_21 y M=m_32
    real(8) :: P                  ! P es el momento del neutrino
    integer :: nu                 ! nu is 1 for neutrinos an 2 for antineutrino    
    real(8) :: Ne                 ! Ne is the electron density

    double complex :: t_Hm(3,3)   ! t_Hm is the Hamiltonian in the mass base
    double complex :: I(3,3)             ! I is the identity matrix
    double complex :: traceHm,trHm! traceHm is the trace of the Hamiltonian-mass

    
    call totalHamiltonianMassBase(t_Hm,t12,t23,t13,delta,sm,aM,P,nu,Ne)
    
    call identityMatrix(I)    

    trHm=traceHm(t12,t23,t13,delta,sm,aM,P,nu,Ne)

    T=t_Hm-((trHm/3.0d0)*I)
    return
end subroutine tMatrix



!################################################
!
!   tMatrix2: is a subroutine that load
!       the T-matrix squared (3,3)
!
!################################################
subroutine tMatrix2(T2,t12,t23,t13,delta,sm,aM,P,nu,Ne)
    implicit none
    double complex :: T2(3,3)            ! T2 is the T-matrix squared
    real(8) :: t12,t23,t13,delta  ! Are the three mixing angles and the CP-violation phase of the mixing matrix
    real(8) :: sm,aM                ! sm,aM are the squared mass difference m=m_21 y M=m_32
    real(8) :: P                  ! P es el momento del neutrino
    integer :: nu                 ! nu is 1 for neutrinos an 2 for antineutrino    
    real(8) :: Ne                 ! Ne is the electron density

    double complex :: T(3,3)             ! T is the T-matrix of the model

    call tMatrix(T,t12,t23,t13,delta,sm,aM,P,nu,Ne)
    T2 = matmul(T,T)
    return
end subroutine tMatrix2



!################################################
!
!   vectorA: is a subroutine that load
!
!################################################

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


!################################################
!
!   matterDensityUnits: is a functions that return
!       the matter density value for neutrinos 
!       (nu=1) or antineutrinos (nu=2) giving 
!       an electron density Ne
!
!################################################
real(8) function matterDensityUnits(nu,Ne)
    implicit none
    integer :: nu                 ! nu is 1 for neutrinos an 2 for antineutrino
    real(8) :: Ne                 ! Ne is the electron density [N_A/cm^{3}]
    
    real(8) :: matterDensity      ! is a functions that return the matter density value [eV]
    

    matterDensityUnits=0.0d0
    matterDensityUnits=5.067731d9*matterDensity(nu,Ne)   ! Without units
    return
end function matterDensityUnits


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


!################################################
!
!   hamiltonianDiff: is a subroutine that 
!       load the hamiltonian-diff-matrix 
!
!################################################
subroutine hamiltonianDiff(Hdiff,sm,aM,P)
    implicit none
    double complex :: Hdiff(3,3)  ! Hdiff is the hamiltonian differential 
    real(8) :: sm,aM              ! sm,aM are the squared mass difference m=m_21 y M=m_32
    real(8) :: P                  ! P es el momento del neutrino o la enerǵia total
    real(8),parameter :: unityFactor=2.5399811853d0

    Hdiff(1,1)=cmplx((-unityFactor*sm)/P,0.0d0);  Hdiff(1,2)=0.0d0;  Hdiff(1,3)=0.0d0
    Hdiff(2,1)=0.0d0;                             Hdiff(2,2)=0.0d0;  Hdiff(2,3)=0.0d0
    Hdiff(3,1)=0.0d0;                             Hdiff(3,2)=0.0d0;  Hdiff(3,3)=cmplx((unityFactor*aM)/P,0.0d0)
    return
end subroutine hamiltonianDiff





!################################################
!
!   hamiltonianInVacuum: is a subroutine that 
!       load the hamiltonian-vacuum-matrix 
!
!################################################
subroutine hamiltonianInVacuum(Hv,sm,aM,P)
    implicit none
    double complex :: Hv(3,3)     ! Hv is the hamiltonian for the propagation of the neutrinos in vacuum
    real(8) :: sm,aM              ! sm,aM are the squared mass difference m=m_21 y M=m_32
    real(8) :: P                  ! P es el momento del neutrino o la enerǵia total

    double complex :: Hdiff(3,3)  ! Hdiff is the hamiltonian differential 
    double complex :: Hdiag(3,3)  ! Hdiag is the hamiltonian differential 
    

    call hamiltonianDiff(Hdiff,sm,aM,P)
    !print*, 'Hdiff',Hdiff(1,:)
    !print*, 'Hdiff',Hdiff(2,:)
    !print*, 'Hdiff',Hdiff(3,:)

    call hamiltonianDiag(Hdiag)
    !print*, 'Hdiag', Hdiag(1,:)
    !print*, 'Hdiag', Hdiag(2,:)
    !print*, 'Hdiag', Hdiag(3,:)
    Hv = Hdiff + Hdiag
    return
end subroutine hamiltonianInVacuum


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

 lambdaMtrx(1,1)=1.0d0;  lambdaMtrx(1,2)=cmplx(0.d0,-L)*Ls(1);  lambdaMtrx(1,3)=-(L**2)*(Ls(1)**2);
 lambdaMtrx(2,1)=1.0d0;  lambdaMtrx(2,2)=cmplx(0.d0,-L)*Ls(2);  lambdaMtrx(2,3)=-(L**2)*(Ls(2)**2);
 lambdaMtrx(3,1)=1.0d0;  lambdaMtrx(3,2)=cmplx(0.d0,-L)*Ls(3);  lambdaMtrx(3,3)=-(L**2)*(Ls(3)**2);

 return
end subroutine lambdaMatrix

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



!################################################
!
!   matterDensity: is a functions that return
!       the matter density value for neutrinos 
!       (nu=1) or antineutrinos (nu=2) giving 
!       an electron density Ne
!
!################################################
real(8) function matterDensity(nu,Ne) ! Units [eV]
    implicit none
    integer :: nu                 ! nu is 1 for neutrinos an 2 for antineutrino
    real(8) :: Ne                 ! Ne is the electron density [N_A/cm^{3}]
    real(8), parameter :: sqrt2GFNa=7.6324D-14! Fermi constant [eV cm^{3}/N_A]
    real(8), parameter :: N_A=6.0221415D23  ! N_A is the Avogadro's number [1/mol]

    matterDensity=0.0d0
    !print*, 'Ne in matterDensity: ', Ne
    matterDensity=(sqrt2GFNa)*Ne ! Utilizar para cuando Ne = [N_A/cm^{3}], en el caso la densidad en el Sol dada por John Bahcall
                               ! los datos vienen dados en  [cm^{-3}/N_A] por lo que hay que multiplicarlos por N_A**2
    select case(nu)
        case(1)
            matterDensity=matterDensity            
        case(2)
            matterDensity=-matterDensity            
        case default
           !print*, nu            
            print*, 'Error: no existe la opción'
            print*,'nu',nu
            print*,'Ne',Ne
            stop
    end select
end function matterDensity





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




!################################################
!
!   potentialMatrixMassBase: is a subroutine 
!       that load the potential Vf matrix 
!
!################################################
subroutine potentialMatrixMassBase(Vm,t12,t23,t13,delta,nu,Ne)
    implicit none
    double complex :: Vm(3,3)     ! Vm is the potential matrix in mass base
    real(8) :: t12,t23,t13        ! Are the mixing angles of hte oscillation model
    real(8) :: delta              ! delta is the phase factor for CP violations
    integer :: nu                 ! nu is 1 for neutrinos an 2 for antineutrino    
    real(8) :: Ne                 ! Ne is the electron density

    double complex :: U(3,3)      ! U is the mixing matrix of the oscillation model    
    double complex :: U_1(3,3)    ! U_1 is the inverse matrix of U U*U_1=I
    double complex :: Vf(3,3)     ! Vf is the hamiltonian for the propagation of the neutrinos in vacuum
    
    call inverseMixingMatrix(U_1,t12,t23,t13,delta)
    !print*, ''
    !print*,'U_1',U_1(1,:)
    !print*,'U_1',U_1(2,:)
    !print*,'U_1',U_1(3,:)
    !print*, ''
    call mixingMatrix(U,t12,t23,t13,delta)
    !print*, ''
    !print*,'U',U(1,:)
    !print*,'U',U(2,:)
    !print*,'U',U(3,:)
    !print*, ''
    call potentialMatrixFlavourBase(Vf,nu,Ne)
    !print*, ''
    !print*,'Vf',Vf(1,:)
    !print*,'Vf',Vf(2,:)
    !print*,'Vf',Vf(3,:)
    !print*, ''
    

    Vm=matmul( matmul(U_1,Vf),U )

    return
end subroutine potentialMatrixMassBase




!################################################
!
!   totalHamiltonianMassBase: is a subroutine
!       that return the the total Hamiltonian 
!       which is the sum of hamiltonian in 
!       vasuum and the potential matrix in
!       mass base
!
!################################################
subroutine totalHamiltonianMassBase(t_Hm,t12,t23,t13,delta,sm,aM,P,nu,Ne)
    implicit none
    double complex :: t_Hm(3,3)   ! t_Hm is the sum of Hamiltonian in vacumm (mass base)  and the potental in the base of mass
    real(8) :: t12,t23,t13,delta  ! Are the three mixing angles and the CP-violation phase of the mixing matrix
    real(8) :: sm,aM                ! sm,aM are the squared mass difference m=m_21 y M=m_32
    real(8) :: P                  ! P es el momento del neutrino
    integer :: nu                 ! nu is 1 for neutrinos an 2 for antineutrino    
    real(8) :: Ne                 ! Ne is the electron density

    double complex :: Hv(3,3)            ! Hv is the hamiltonian for the propagation of the neutrinos in vacuum
    double complex :: Vm(3,3)            ! Vm is the potential matrix in mass base

    call hamiltonianInVacuum(Hv,sm,aM,P)
    !print*, ''
    !print*,'Hv',Hv(1,:)
    !print*,'Hv',Hv(2,:)
    !print*,'Hv',Hv(3,:)
    !print*, ''
    call potentialMatrixMassBase(Vm,t12,t23,t13,delta,nu,Ne)    
    !print*, ''
    !print*,'Vm',Vm(1,:)
    !print*,'Vm',Vm(2,:)
    !print*,'Vm',Vm(3,:)
    !print*, ''
    t_Hm=Hv+Vm
    return
end subroutine totalHamiltonianMassBase



!################################################
!
!   traceHm: is a function that return the
!       trace of the Hamiltonian-mass matrix
!
!################################################
double complex function traceHm(t12,t23,t13,delta,sm,aM,P,nu,Ne)
    implicit none
    real(8) :: t12,t23,t13,delta  ! Are the three mixing angles and the CP-violation phase of the mixing matrix
    real(8) :: sm,aM                ! sm,aM are the squared mass difference m=m_21 y M=m_32
    real(8) :: P                  ! P es el momento del neutrino
    integer :: nu                 ! nu is 1 for neutrinos an 2 for antineutrino    
    real(8) :: Ne                 ! Ne is the electron density

    double complex :: t_Hm(3,3)          ! t_Hm is the sum of Hamiltonian in vacumm (mass base)  and the potental in the base of mass
    
    traceHm=(0.0D0,0.0d0)

    call totalHamiltonianMassBase(t_Hm,t12,t23,t13,delta,sm,aM,P,nu,Ne)
    traceHm=t_Hm(1,1)+t_Hm(2,2)+t_Hm(3,3)
    return
end function traceHm



subroutine vectorE(vectE,L,t12,t23,t13,delta,sm,aM,P,nu,Ne)
    implicit none
    double complex :: vectE(3)    ! vectE is the vector of the model
    real(8) :: L                  ! L is the length between the source of neutrinos an the position

    double complex :: Ls(3)       ! Ls is an array with the values of coefficients lambda    
    real(8) :: t12,t23,t13,delta  ! Are the three mixing angles and the CP-violation phase of the mixing matrix
    real(8) :: sm,aM              ! sm,aM are the squared mass difference m=m_21 y M=m_32
    real(8) :: P                  ! P es el momento del neutrino
    integer :: nu                 ! nu is 1 for neutrinos an 2 for antineutrino    
    real(8) :: Ne                 ! Ne is the electron density
    
    call lambdaFromEISPACK(Ls,t12,t23,t13,delta,sm,aM,P,nu,Ne)

    vectE(1)=exp(cmplx(0.0d0,-L)*Ls(1))
    vectE(2)=exp(cmplx(0.0d0,-L)*Ls(2))
    vectE(3)=exp(cmplx(0.0d0,-L)*Ls(3))

    return
end subroutine vectorE