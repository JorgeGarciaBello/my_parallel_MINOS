subroutine MINOS_chi_minimization(XX,chi)
    use types
    use MINOS_data, only: NUM_U
    implicit none
    real(dp) :: t13,dm,chi
    real(dp) :: YX
    real(dp) :: X(NUM_U)   ! X => Pulls
    real(dp) :: XX(*)     ! XX = Parámetros de oscilación
    real(dp) :: W(10000),RHOEND,RHOBEG
    INTEGER  :: IPRINT,MAXFUN,N,NPT,I
    IPRINT=1
    MAXFUN=5000
    RHOEND=1.0D-8
    N=NUM_U
    NPT=2*N+1

      DO 10 I=1,N
   10 X(I)=DFLOAT(I)/DFLOAT(N+1)
      RHOBEG=0.2D0*X(1)      
      CALL NEWUOA_M (N,NPT,X,RHOBEG,RHOEND,IPRINT,MAXFUN,W,XX,YX)
      chi=YX
    return
end subroutine MINOS_chi_minimization