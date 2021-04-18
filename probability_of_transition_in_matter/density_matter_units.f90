!################################################
!
!   matterDensityUnits: is a functions that return
!       the matter density value for neutrinos 
!       (nu=1) or antineutrinos (nu=2) giving 
!       an electron density Ne
!
!################################################
real(8) function matterDensityUnits(nu,Ne)
    implicit none
    integer :: nu                 ! nu is 1 for neutrinos an 2 for antineutrino
    real(8) :: Ne                 ! Ne is the electron density [N_A/cm^{3}]
    
    real(8) :: matterDensity      ! is a functions that return the matter density value [eV]
    

    matterDensityUnits=0.0d0
    matterDensityUnits=5.067731d9*matterDensity(nu,Ne)   ! Without units
    return
end function matterDensityUnits