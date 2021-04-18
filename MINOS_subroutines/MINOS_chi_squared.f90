function MINOS_chi_squared(Y,P)
    use types
    use MINOS_data, only: NBIN,NUM_U,best_fit,bins,bottom_sigma_limits,data_points,no_oscillation, &
                          sys_uncertainties, sigma,bb
    implicit none
    real(dp) :: MINOS_chi_squared
    real(dp) :: average_probability_per_bin
    real(dp) :: mod
    real(dp) :: obs,exp
    real(dp) :: Y(*)
    real(dp) :: P(NUM_U)
    real(dp) :: chi
    integer  ::  i,j,k,u

    real(dp) :: t12,t23,t13,delta          ! angle mixing and phase
    real(dp) :: sm,aM                      ! sm,aM are the squared mass difference m=m_21 y M=m_32
    real(dp) :: x                          ! Energy [GeV]
    integer  :: nu                         ! nu is 1 for neutrinos an 2 for antineutrino
    integer  :: flvr1                      ! flavor intial 1 2 3 = e, m, t
    integer  :: flvr2                      ! flavor final 1 2 3 = e, m, t
    real(dp) :: L                          ! L is the length of slabs in serie. unitis [km]    
    real(dp) :: ro                         ! ro [g/cm^3]
    integer  :: tested_bin
    
    t12=Y(1)
    t23=Y(2)
    t13=Y(3)
    delta=Y(4)
    sm=Y(5)
    aM=Y(6)
    nu=1
    chi=0.0_dp    
        
    !tested_bin=bb
    do i=1,NBIN
        if(i==1) then
            obs= abs(data_points(i,2) + P(1))
        elseif(i==2) then
            obs= abs(data_points(i,2) + P(2))
        elseif(i==3) then
            obs= abs(data_points(i,2) + P(3))
        elseif(i==4) then
            obs= abs(data_points(i,2) + P(4))
        elseif(i==5) then
            obs= abs(data_points(i,2) + P(5))
        elseif(i==6) then
            obs= abs(data_points(i,2) + P(6))
        elseif(i==7) then
            obs= abs(data_points(i,2) + P(7))
        elseif(i==8) then
            obs= abs(data_points(i,2) + P(8))
        else
           obs=data_points(i,2)
        endif

        exp=no_oscillation(i,2)*average_probability_per_bin(i,t12,t23,t13,delta,sm,aM,nu)*abs( 1.0_dp+P(9) )
        !exp=no_oscillation(i,2)*average_probability_per_bin(i,t12,t23,t13,delta,sm,aM,nu)*abs( 1.0_dp+P(1)+P(2) )
        !print*,'chi_i',i, data_points(i,2),P(1),exp,obs/exp
        !chi = chi + ( exp - obs + obs*LOG(obs/exp) ) + (P(i)/sigma(i,2))**2
        !print*, 'chi_i', i, ( exp - obs + obs*LOG(obs/exp) )
        chi = chi + ( exp - obs + obs*LOG(obs/exp) )
    enddo

    do i=1,8
    !chi = chi + 0.5*( P(i)/sigma(tested_bin,2) )**2
    chi = chi + 0.5*( P(i)/sigma(i,2) )**2
    enddo

    !do i=1,5
    !chi = chi + 0.5_dp*( P(1)/(sys_uncertainties(1)/100.0_dp) )**2
    chi = chi + 0.5_dp*( P(9)/(sys_uncertainties(4)/100.0_dp) )**2
    !enddo
    print*, 'chi', chi
    MINOS_chi_squared=chi
    return
end function MINOS_chi_squared