!################################################
!
!   matterDensity: is a functions that return
!       the matter density value for neutrinos 
!       (nu=1) or antineutrinos (nu=2) giving 
!       an electron density Ne
!
!################################################
!real(8) function matterDensity(nu,Ne,ro,Z,A) ! Units [eV]
real(8) function matterDensity(nu,Ne)
    implicit none
    integer :: nu                 ! nu is 1 for neutrinos an 2 for antineutrino
    real(8) :: Ne                 ! Ne is the electron density [N_A/cm^{3}]
    real(8) :: ro                 ! ro is the matter density
    integer :: Z,A                ! Z is the number of electrons of the nucleo and A is the molar mass [g/mol]    
    real(8), parameter :: N_A=6.0221415D23  ! N_A is the Avogadro's number [1/mol]}
    real(8), parameter :: GF=8.96180870D-38  ! GF is the Fermi constatnt [ eV cm^{3} ]
    real(8), parameter :: sqrt2=sqrt(2.d0)  ! GF is the Fermi constatnt [ eV cm^{3} ]

    matterDensity=0.0d0    
    select case(nu)
        case(1)            
            matterDensity = sqrt2*GF*Ne
        case(2)
            matterDensity = -sqrt2*GF*Ne
        case default
            print*, nu            
            print*, 'Error: no existe la opcion-matterDensity=0.0'
            print*,'nu',nu
            print*,'Ne',Ne
            stop
    end select
end function matterDensity