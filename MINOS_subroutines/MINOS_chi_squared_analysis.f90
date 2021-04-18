subroutine MINOS_chi_squared_analysis()
    use types
    use constants
    use MINOS_data
    implicit none
    integer,parameter :: m=15
    integer,parameter :: n=50    
    integer,parameter :: l=50
    real(dp) :: MINOS_chi_squared
    real(dp) :: X(NUM_U)   ! X => Pulls
    real(dp) :: results(m,n,l)
    real(dp) :: res(m)
    real(dp) :: t12,t23,t13,delta          ! angle mixing and phase
    real(dp) :: sm,aM                      ! sm,aM are the squared mass difference m=m_21 y M=m_32
    integer :: nu                         ! nu=1 neutrinos - nu=2 antineutrinos
    real(dp) :: Y(12)
    real(dp) :: P(NUM_U)
    real(dp) :: chi
    integer  :: i,j,k,u
    real(dp) :: t23_i,t23_f,dm32_i,dm32_f,t13_i,t13_f
    real(dp) :: t23_p(m),dm32_p(n),t13_p(l)
    real(dp) :: average_probability_per_bin,chis

    nu=1
    Y(1)=asin(sqrt(0.307_dp)) ! 0.5872523687_dp   ! t12=0.5872523687_dp
    !Y(2)=t23    !t23=0.8304591359_dp
    Y(3)=asin(sqrt(0.021_dp)) ! 0.1481900178_dp   ! t13=0.1481900178_dp
    Y(4)=0.0_dp            ! delta=0.0_dp
    Y(5)=7.54d-5           ! sm=7.53d-5
    Y(6)=2.5d-3            ! aM=2.5d-3

    t23_i=0.5_dp
    t23_f=1.1_dp
    do i=1,m
        t23_p(i)=t23_i + real(i-1)*(t23_f-t23_i)/real(m)
    enddo    
    
    !!$omp parallel do private(chi,Y)
    do i=1,m
         Y(1)=asin(sqrt(0.307_dp)) ! 0.5872523687_dp   ! t12=0.5872523687_dp
        !Y(2)=t23    !t23=0.8304591359_dp
        Y(3)=asin(sqrt(0.021_dp)) ! 0.1481900178_dp   ! t13=0.1481900178_dp
        Y(4)=0.0_dp            ! delta=0.0_dp
        Y(5)=7.54d-5           ! sm=7.53d-5
        Y(6)=2.5d-3            ! aM=2.5d-3
        !Y(6)=2.48d-3            ! aM=2.48d-3
        print*, i        
            Y(2) = t23_p(i) !0.28215877686164392_dp !t23_p(i)
            !Y(2) = asin(sqrt(0.38_dp)) ! t23
            call MINOS_chi_minimization(Y,chi)
            results(i,1,1) = chi            
            !do i=1,NBIN
    !print*,data_points(i,2), &
           !no_oscillation(i,2)*average_probability_per_bin(i,Y(1),Y(2),Y(3),Y(4),Y(5),Y(6),nu)*abs( 1.0_dp+X(1)+X(2) )
            !enddo
    enddo
    !!$omp end parallel do

    open(newunit=u,file='MINOS_data/results.dat',position='append',status='old')
    do i=1,m
        Y(2) = t23_p(i)
        write(u,*) sin(Y(2))**2, results(i,1,1) 
    enddo
    write(u,*) ' '
    write(u,*) ' '
    close(u)

    res = results(:,1,1)
    print*, minval(res)
    return
end subroutine MINOS_chi_squared_analysis