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