!##################################################################################################################
!
!       iterativeProbabilityAmplitude: is the amplitud of probability for a neutrino that
!           was born is the zone 'zone' of the Sun with density a_density and energy P
!
!##################################################################################################################
double complex function iterativeProbabilityAmplitude(dim,flvr1,flvr2,l_a,ro_a,t12,t23,t13,delta,sm,aM,P,nu,Z,A)
    use types
    implicit none
    integer  :: dim
    integer  :: flvr1               ! flvr1 is the flavour with which the neutrino is generated
    integer  :: flvr2               ! flvr2 is the flavour that is transited
    real(dp) :: l_a(dim)            ! L is the length between the source of neutrinos an the position
    real(dp) :: ro_a(dim)           ! is an array with matter density values
    real(dp) :: t12,t23,t13,delta   ! Are the three mixing angles and the CP-violation phase of the mixing matrix
    real(dp) :: sm,aM               ! sm,aM are the squared mass difference m=m_21 y M=m_32
    real(dp) :: P                   ! P es el momento del neutrino
    integer  :: nu                  ! nu is 1 for neutrinos an 2 for antineutrino
    integer  :: Z                   ! Z
    integer  :: A                   ! A
    double complex :: iterUfL(3,3)! iterUfL is the time evolution operator matrix in the flavour base iterative for diferents lengths and electron densities
    real(dp) :: fState1(3)         ! is the initial flavour eigenstate
    real(dp) :: fState2(3)         ! ! is the final flavour eigenstate
    iterativeProbabilityAmplitude=0.0d0
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
    call iterativeTimeEvolutionOperator(dim,iterUfL,P,l_a,ro_a,t12,t23,t13,delta,sm,aM,nu,Z,A)    
    iterativeProbabilityAmplitude = dot_product(fState1,matmul(iterUfL,fState2))
    return
end function iterativeProbabilityAmplitude