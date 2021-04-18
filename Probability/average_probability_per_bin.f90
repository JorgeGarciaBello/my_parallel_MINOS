function average_probability_per_bin(bin,t12,t23,t13,delta,sm,aM,nu)
    use types
    use MINOS_data, only: NBIN, bins
    implicit none
    integer,parameter  :: n=50
    real(dp) :: average_probability_per_bin
    integer  :: bin    
    real(dp) :: t12,t23,t13,delta          ! angle mixing and phase
    real(dp) :: sm,aM                      ! sm,aM are the squared mass difference m=m_21 y M=m_32    
    integer  :: nu                         ! nu is 1 for neutrinos an 2 for antineutrino

    real(dp) :: probability_of_transition_in_matter_a_b
    integer  :: flvr1                      ! flavor intial 1 2 3 = e, m, t
    integer  :: flvr2                      ! flavor final 1 2 3 = e, m, t
    real(dp) :: L                          ! L is the length of slabs in serie. unitis [km]        
    real(dp) :: ro                         ! ro [g/cm^3]
    integer  :: Z                          ! Z = 1
    integer  :: A                          ! A = 2    
    real(dp) :: result,result1
    
    real(dp) :: x                          ! Energy [GeV]
    real(dp) :: a_bin,b_bin,r,h
    integer  :: k
    
    flvr1=2; flvr2=2
    L=735.0_dp             ! [km] far detector minos for de source
    ro=3.3_dp              ! ro [g/cm^3]
    Z=1; A=2;    

    a_bin=bins(bin,1);    b_bin=bins(bin,2)
    h=(b_bin-a_bin)/real(n)    
    r=0.0_dp
    do k=1,n
        x=a_bin + h*real(k-1)
        if (bin==1) then
            call logarithmic_partition(n,bin,LOG10(a_bin),LOG10(b_bin),x)
        end if
        r = r + h*( probability_of_transition_in_matter_a_b(flvr1,flvr2,L,t12,t23,t13,delta,sm,aM,x,nu,ro,Z,A) + &
                    probability_of_transition_in_matter_a_b(flvr1,flvr2,L,t12,t23,t13,delta,sm,aM,x+h,nu,ro,Z,A) )/real(2)
    enddo
    average_probability_per_bin = r/(b_bin-a_bin)
    return
end function average_probability_per_bin