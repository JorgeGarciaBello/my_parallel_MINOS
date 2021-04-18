program minos_main
    use types
    use minos_data
    implicit none
    integer :: bin
    integer :: i,j
    real(dp) :: t12,t23,t13,delta  ! angle mixing and phase
    real(dp) :: sm,aM              ! sm,aM are the squared mass difference m=m_21 y M=m_32    
    integer  :: nu                 ! nu is 1 for neutrinos an 2 for antineutrino            
    real(dp) :: average_probability_per_bin
    real(dp) :: probability_of_transition_in_matter_a_b
    real(dp) :: results,x
    call MINOS_read_data()

    t12=0.5872523687_dp
    t23=0.8304591359_dp
    t13=0.1481900178_dp
    delta=0.0_dp
    sm=7.53d-5
    aM=2.5d-3
    nu=1
    !bin=1

    !do bb=1,39
    call MINOS_gaussian_analysis()
        !call MINOS_chi_squared_analysis()
    !enddo
    !do bin=1,NBIN        
        !print*, bin, bins(bin,:), average_probability_per_bin(bin,t12,t23,t13,delta,sm,aM,nu)
    !enddo
    !x=0.00001
    !do bin=1,1000
        !x=bin*(50.0_dp)/1000
        !call logarithmic_partition(1000,bin,LOG10(0.0001_dp),LOG10(50.0_dp),x)
        !print*, x,probability_of_transition_in_matter_a_b(2,2,735.0_dp,t12,t23,t13,delta,sm,aM,x,nu,3.3_dp,1,2)
    !enddo
end program