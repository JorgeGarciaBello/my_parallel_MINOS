      SUBROUTINE CALFUN_M (N,X,F,XX)
      use types
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION X(*),XX(*)      
      real(dp) :: MINOS_chi_squared
      ! XX = Parámetros de oscilación
      ! X => Pulls
      F=MINOS_chi_squared(XX,X)      
      !F=X(1)**2 + X(2)**2+XX(3)      
      RETURN
      END