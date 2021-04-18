module constants
    use types
    implicit none
    !integer,  parameter :: dim    = 5                            ! Number of Earth density layers [3,5]
    real(dp), parameter :: PI     = 4._dp*DATAN(1._dp)
    real(dp), parameter :: N_A    = 6.0221415D23                 ! N_A is the Avogadro's number [1/mol]
    real(dp), parameter :: GF     = 8.96180870D-38               ! GF is the Fermi constatnt [ eV cm^{3} ]    
end module constants