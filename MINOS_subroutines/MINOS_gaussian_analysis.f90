subroutine MINOS_gaussian_analysis()
    use types
    use constants
    use MINOS_data
    implicit none
    integer,parameter :: m=30
    integer,parameter :: n=30
    integer,parameter :: l=30
    real(dp) :: MINOS_gaussian_chi_squared    
    real(dp) :: results(m,n,l)
    real(dp) :: results_sk(44800)
    real(dp) :: t12,t23,t13,delta          ! angle mixing and phase
    real(dp) :: Y(12)
    real(dp) :: chi    
    real(dp) :: dm32_i,dm32_f
    real(dp) :: t13_i,t13_f
    real(dp) :: t23_i,t23_f
    real(dp) :: t23_p(m),dm32_p(n),t13_p(l)
    integer  :: i,j,k,u

    real(dp) :: var_the23min, var_the23max
    real(dp) :: var_the13min, var_the13max
    real(dp) :: var_dm23min,  var_dm23max
    !###############################################################
    !
    !   Setting regions of the parameters
    !
    !###############################################################
    var_the23min = asin(sqrt(0.35)) ! 0.5
    var_the23max = asin(sqrt(0.65)) ! 1.0 

    var_the13min=-0.45
    var_the13max=0.45
 
    var_dm23min=2.20d-3
    var_dm23max=2.70d-3
    

    !t23_i=asin(sqrt(0.28_dp))
    !t23_f=asin(sqrt(0.75_dp))
    t23_i=var_the23min
    t23_f=var_the23max
    do i=1,m
        t23_p(i)=t23_i + real(i-1)*(t23_f-t23_i)/real(m)
    enddo
    !dm32_i=0.0021_dp
    !dm32_f=0.0028_dp
    dm32_i=var_dm23min
    dm32_f=var_dm23max
    do i=1,n
        dm32_p(i)=dm32_i + real(i-1)*(dm32_f-dm32_i)/real(m)
    enddo
    !t13_i=asin(sqrt(0.0000_dp))
    !t13_f=asin(sqrt(0.0775_dp))
    t13_i=var_the13min
    t13_f=var_the13max
    do i=1,l
        t13_p(i)=t13_i + real(i-1)*(t13_f-t13_i)/real(m)
    enddo
    print*, 'hola'   
!###############################################################    
    do i=1,m
        print*, i
        !$omp parallel do private(Y)
        do j=1,n
            Y(1) = 0.5873            ! t12
            !Y(2) = 0.6642           ! t23
            Y(3) = 0.0d0 ! 0.1454    ! t13
            Y(4) = 0.0d0             ! delta
            Y(5) = 7.540d-5          ! sm = dm21
            !Y(6) = 2.480d-3         ! aM = dm32 
            Y(7) = 0.0d0
            Y(8) = 0.0d0
            Y(9) = 0.0d0
            Y(10)= 0.0

            !Y(1)=asin(sqrt(0.307_dp)) ! 0.5872523687_dp   ! t12=0.5872523687_dp    
            Y(2) = t23_p(i)      ! t23_p(i)
            !Y(3)=asin(sqrt(0.021_dp)) ! 0.1481900178_dp   ! t13=0.1481900178_dp
            !Y(4)=0.0_dp            ! delta=0.0_dp
            !Y(5)=7.54d-5           ! sm=7.53d-5
            Y(6) = dm32_p(j) ! dm32_p(j)
            results(i,j,1) = MINOS_gaussian_chi_squared(Y)
        enddo
        !$omp end parallel do
    enddo
    open(newunit=u,file='MINOS_data/pp.MINOS.gauss.chi.dat')
    do i=1,m
        do j=1,n
            write(u,*) sin(t23_p(i))**2, dm32_p(j), results(i,j,1) 
        enddo
        write(u,*) ' '
    enddo
    close(u)
    print*, MINVAL(results(:,:,1))

!###############################################################
!!$omp parallel do private(Y)
!    do i=1,44800
!        print*, i
!        Y(1) = 0.5873            ! t12
!        Y(2) = asin(sqrt(grid_sk(i,3)))      ! t23
!        Y(3) = asin(sqrt(grid_sk(i,2)))      ! t13
!        Y(4) = 0.0d0             ! delta
!        Y(5) = 7.540d-5          ! sm = dm21
!        Y(6) = grid_sk(i,1)      ! aM = dm32 
!        results_sk(i) = MINOS_gaussian_chi_squared(Y)
!    enddo
!!$omp end parallel do
!open(newunit=u,file='MINOS_data/3D.MINOS.gauss.chi.dat')
!do i=1,44800    
!    write(u,*) grid_sk(i,1), grid_sk(i,2), grid_sk(i,3), results_sk(i)
!enddo
!close(u)
!print*, MINVAL(results_sk)
!###############################################################    
end subroutine MINOS_gaussian_analysis