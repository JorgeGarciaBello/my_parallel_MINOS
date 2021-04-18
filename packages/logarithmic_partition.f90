!##################################################################################
!
!   logarithmic_partition: es una subroutina que
!       genera puntos distribuidos logaritmicamente 
!       en un ragno determinado.
!
!   input: 
!       n  (int)  => es el tamaño de la partición
!       i  (int)  => es la posición del punto a generar dentro de la partición
!       v1 (real) => es el límite inferior del intervalo a partir,
!                    Debe ser dado como el exponenete de la base 10:
!                    10^{8} => v1 = 8.
!       v2 (real) => es el límite superior del intervalo a partir,
!                    Debe ser dado como el exponenete de la base 10:
!                    10^{15} => v1 = 15.
!   output:
!       x (real)  => es el resultado de a partición en la posición i
!
!
!##################################################################################
subroutine logarithmic_partition(n,i,v1,v2,x)
    use types
    implicit none
    integer  :: n,i
    real(dp) :: v1,v2,x
    real(dp) :: arg

    arg=v1 + i*((v2-v1)/real(n))
    x = exp(log(10._dp)*(arg))
    return
end subroutine logarithmic_partition