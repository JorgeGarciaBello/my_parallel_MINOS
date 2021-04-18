function MINOS_gaussian_chi_squared(Y)
    use types
    use MINOS_data, only: NBIN,NUM_U,best_fit,bins,bottom_sigma_limits,data_points,no_oscillation, &
                          sys_uncertainties, sigma,bb                          
    implicit none
    real(dp) :: MINOS_gaussian_chi_squared
    real(dp) :: Y(*)

    real(dp) :: average_probability_per_bin    
    real(dp) :: obs,exp,expd(NBIN)
    real(dp) :: chi
    real(dp) :: t12,t23,t13,delta          ! angle mixing and phase
    real(dp) :: sm,aM                      ! sm,aM are the squared mass difference m=m_21 y M=m_32        
    integer  :: nu                         ! nu is 1 for neutrinos an 2 for antineutrino
    real(dp) :: r                          ! r is the resulted probability from 'minos_prop_one_slab'
    integer  :: i,j,k,bin,u
    real(dp) :: epsilon1,Ener,P_aver,Pmumu,Prob_k(NBIN),NOS_Minos(NBIN),min_osc(NBIN)
    real(dp) :: int_min_data, int_min_nosc, int_min_th
    real(dp) :: Ener_cent, min_chi2_pull,min_chi2,sigma_k
    real(dp) :: min_tot_th,min_tot_nosc,min_tot_data,min_th(NBIN),min_nosc(NBIN),min_data(NBIN)

    min_osc=data_points(:,2)
    
    !#####################################################
    !
    !   Obtención de los valores promedios de 
    !       Eventos sin oscilación en MINOS => NOS_Minos
    !
    !#####################################################
    t12   = 0.5873
    t23   = 0.6642
    t13   = 0.1454    
    delta = 0.0d0
    sm    = 7.54d-5
    aM    = 2.480d-3
    nu=1
    epsilon1=0.005d0
    do k=1,39
        Ener=bins(k,1)
        P_aver=0.0d0
        do  while (Ener.le.bins(k,2))
            call minos_prop_one_slab(2,2,t12,t23,t13,delta,sm,aM,Ener,nu,r)
            Pmumu=r
            P_aver=P_aver+Pmumu*(epsilon1)  
            Ener=Ener+epsilon1
        enddo
        P_aver       = P_aver/(bins(k,2)-bins(k,1))
        Prob_k(k)    = P_aver
        NOS_Minos(k) = best_fit(k,2)/Prob_k(k)    ! min_osc(k)/Prob_k(k)
    enddo
    !################################################
    !
    !  Calibrations to fix the confidence region       !
    !
    !################################################       
    do k=1,5
        NOS_Minos(k)=NOS_Minos(k)*0.9!*0.90d0  ! FIxed to 0.7 hay que porbar 0.98y menoer hasta 0.95s
    enddo
    do k=6,10
        NOS_Minos(k)=NOS_Minos(k)*0.8!*0.8d0  ! FIxed to 1.2 !0.6,0.8 casi funciona, manditne minimo y valores de masa y angulo, pero neceitamosaplanar la curva
    enddo
    do k=11,15
        NOS_Minos(k)=NOS_Minos(k)*1.07!*1.07d0  ! Fixed 1.07
    enddo

    do k=21,25
        NOS_Minos(k)=NOS_Minos(k)*0.98d0  ! Fixed 0.98 ! 0.89
    enddo
    do k=26,30
        NOS_Minos(k)=NOS_Minos(k)*0.98d0  ! Fixed 0.98 ! 0.89
    enddo
    do k=31,35
        NOS_Minos(k)=NOS_Minos(k)*0.98d0  ! Fixed 0.98 ! 0.89
    enddo
    do k=36,39
        NOS_Minos(k)=NOS_Minos(k)*0.98d0  ! Fixed 0.98 ! 0.89
    enddo
    !##############################################################
    !
    !
    !
    !
    !##############################################################
    min_tot_th=0.0d0
    min_tot_nosc=0.0d0
    min_tot_data=0.0d0

    do k=1,39
        Ener=bins(k,1)
        int_min_th=0.0d0
        int_min_nosc=0.0d0
        int_min_data=0.0d0
        P_aver=0.0d0

        if  (  k.eq.12 ) then 
        Ener= 0.999* Ener
        endif 
        if  (  k.eq.11) then 
        Ener= 0.996* Ener
        endif 
        if  (  k.eq.10 ) then 
        Ener= 0.995* Ener
        endif 
        if  (  k.eq.9) then 
        Ener= 0.99* Ener
        endif          
        if  (  k.eq.8) then 
        Ener= 0.98* Ener
        endif                
        if  (  k.eq.7) then 
        Ener= 0.97* Ener
        endif           
        if  (  k.eq.6) then 
        Ener= 0.93 * Ener
        endif                  
        if  (  k.eq.5) then 
        Ener= 0.80 * Ener
        endif                     
        if  (  k.eq.4) then 
        Ener= 0.78 * Ener
        endif
        if  (  k.eq.3) then
        Ener= 0.90* Ener
        endif

        t12=Y(1)
        t23=Y(2)
        t13=Y(3)
        delta=Y(4)
        sm=Y(5)
        aM=Y(6)
        
        do  while ( Ener.le.bins(k,2) )
            call minos_prop_one_slab(2,2,t12,t23,t13,delta,sm,aM,Ener,1,r)
            Pmumu = r            
            int_min_th   = int_min_th   + NOS_Minos(k) * Pmumu*( epsilon1 )
            int_min_nosc = int_min_nosc + NOS_Minos(k) *       ( epsilon1 )
            int_min_data = int_min_data + min_osc(k)   *       ( epsilon1 )

            P_aver=P_aver+Pmumu*(epsilon1)
            Ener=Ener+epsilon1
        enddo
        P_aver = P_aver/( bins(k,2) - bins(k,1))
        
        min_th(k)   = int_min_th     !  / (bin_max(k)-bin_min(k))
        min_nosc(k) = int_min_nosc   !  / (bin_max(k)-bin_min(k))
        min_data(k) = int_min_data   !  / (bin_max(k)-bin_min(k))       

        Prob_k(k)    = P_aver
        min_tot_th   = min_tot_th+min_th(k)
        min_tot_nosc = min_tot_nosc+min_nosc(k)
        min_tot_data = min_tot_data+min_data(k)

    enddo

    min_chi2_pull=0.0d0
    do k=1,39              
        sigma_k= min_data(k)                                           ! * sigma(k)  !min_th(k)  !   !*   
        if(k<=10) sigma_k = sigma_k*13.0                                 ! FIxed value 10.0 && 12.0 en 20
            min_chi2_pull =   min_chi2_pull &
                            !+ ( min_th(k)*0.99d0 -  min_data(k)  )**2   / (sigma_k)
                            + ( min_th(k)*0.99 -  min_data(k)  )**2   / (sigma_k)
    enddo
    min_chi2=min_chi2_pull
    !print*, min_chi2
!################################################
!    do i=1,NBIN
!        expd(i)=no_oscillation(i,2)*average_probability_per_bin(i,t12,t23,t13,delta,sm,aM,nu)
!    enddo
!    chi=0.0_dp
!    do i=1,NBIN
!        exp=expd(i)
!        obs=data_points(i,2)
!        chi = chi + ( 0.98_dp*exp - obs )**2/(obs*2.12_dp)
!        !print*, exp,obs,chi
!    enddo
!    MINOS_gaussian_chi_squared=chi
    MINOS_gaussian_chi_squared=min_chi2
    return
end function MINOS_gaussian_chi_squared