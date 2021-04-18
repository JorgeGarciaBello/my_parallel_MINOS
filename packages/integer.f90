subroutine integer(f,n,a,b,r)
    use types
    implicit none
    real(dp), external :: f
    integer  :: n
    real(dp) :: a,b,r,h    
    integer  :: i
    real(dp) :: x 

    h=(b-a)/real(n)
    r=0.0_dp
    do i=1,n
        x=a + h*real(i-1)
        r = r + h*( f(x) + f(x+h) )/real(2)
    enddo
    return
end subroutine integer