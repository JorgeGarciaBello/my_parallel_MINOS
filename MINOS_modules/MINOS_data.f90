module MINOS_data
    use types    
    implicit none
    integer, parameter :: NUM_U=9                           ! Amout of uncertainties
    integer, parameter :: NBIN=39                         ! Amout of energy bins
    integer, parameter :: NBIN_beam=24                    ! Amout of energy bins
    integer, parameter :: NBIN_nonf_beam=7                ! Amout of energy bins for non fiducial beam
    integer, parameter :: NBIN_cv_beam=11                 ! Amout of energy bins for contained vertex beam
    integer, parameter :: NBIN_enh_beam=13                ! Amout of energy bins for enhanced beam
    integer :: bb
    
    real(dp) :: grid_sk(44800,3)
    real(dp) :: best_fit(NBIN,2)                ! Amout of energy bins
    real(dp) :: bins(NBIN,2)                    ! Amout of energy bins
    real(dp) :: bottom_sigma_limits(NBIN,2)
    real(dp) :: data_points(NBIN,2)
    real(dp) :: no_oscillation(NBIN,2)
    real(dp) :: POT_numu_minos_mas(NBIN,2)
    real(dp) :: POT_numu_numubar(NBIN,2)
    real(dp) :: upper_sigma_limits(NBIN,2)
    real(dp) :: sigma(NBIN,2)
    

    real(dp) :: numu_beam_best_fit(NBIN_beam,2)
    real(dp) :: numu_beam_bins(NBIN_beam,2)
    real(dp) :: numu_beam_bottom_sigma(NBIN_beam,2)
    real(dp) :: numu_beam_data_points(NBIN_beam,2)
    real(dp) :: numu_beam_NC_background(NBIN_beam,2)
    real(dp) :: numu_beam_no_oscillations(NBIN_beam,2)
    real(dp) :: numu_beam_upper_sigma(NBIN_beam,2)

    real(dp) :: beam_nonf_mu_best_fit(NBIN_nonf_beam,2)
    real(dp) :: beam_nonf_mu_bins(NBIN_nonf_beam,2)
    real(dp) :: beam_nonf_mu_bottom_sigma(NBIN_nonf_beam,2)
    real(dp) :: beam_nonf_mu_data_points(NBIN_nonf_beam,2)
    real(dp) :: beam_nonf_mu_NC_background(NBIN_nonf_beam,2)
    real(dp) :: beam_nonf_mu_no_oscillations(NBIN_nonf_beam,2)
    real(dp) :: beam_nonf_mu_upper_sigma(NBIN_nonf_beam,2)

    real(dp) :: beam_cv_best_fit(NBIN_cv_beam,2)
    real(dp) :: beam_cv_bins(NBIN_cv_beam,2)
    real(dp) :: beam_cv_bottom_sigma(NBIN_cv_beam,2)
    real(dp) :: beam_cv_data_points(NBIN_cv_beam,2)
    real(dp) :: beam_cv_NC_background(NBIN_cv_beam,2)
    real(dp) :: beam_cv_no_oscillations(NBIN_cv_beam,2)
    real(dp) :: beam_cv_upper_sigma(NBIN_cv_beam,2)

    real(dp) :: beam_enh_best_fit(NBIN_enh_beam,2)
    real(dp) :: beam_enh_bins(NBIN_enh_beam,2)
    real(dp) :: beam_enh_bottom_sigma(NBIN_enh_beam ,2)
    real(dp) :: beam_enh_data_points(NBIN_enh_beam,2)    
    real(dp) :: beam_enh_no_oscillations(NBIN_enh_beam,2)
    real(dp) :: beam_enh_upper_sigma(NBIN_enh_beam,2)
    
    real(dp) :: sys_uncertainties(5)

end module MINOS_data