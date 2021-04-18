subroutine MINOS_read_data()
    use MINOS_data
    implicit none    
    integer :: i,j,u

    open(newunit=u,file='MINOS_data/sk-nh.dat',status='old')
        read(u,*) ((grid_sk(i,j), j=1,3), i=1,44800)
    close(u)

    open(newunit=u,file='MINOS2/MINOS-data/MINOS-best-fit.dat',status='old')
        read(u,*) ((best_fit(i,j), j=1,2), i=1,NBIN)
    close(u)
    open(newunit=u,file='MINOS2/MINOS-data/MINOS-bines.dat',status='old')
        read(u,*) ((bins(i,j), j=1,2), i=1,NBIN)
    close(u)
    open(newunit=u,file='MINOS2/MINOS-data/MINOS-bottom_sigma_limits.dat',status='old')
        read(u,*) ((bottom_sigma_limits(i,j), j=1,2), i=1,NBIN)
    close(u)
    open(newunit=u,file='MINOS2/MINOS-data/MINOS-data-points.dat',status='old')
        read(u,*) ((data_points(i,j), j=1,2), i=1,NBIN)
    close(u)
    open(newunit=u,file='MINOS2/MINOS-data/MINOS-no-oscillation.dat',status='old')
        read(u,*) ((no_oscillation(i,j), j=1,2), i=1,NBIN)
    close(u)
    open(newunit=u,file='MINOS2/MINOS-data/MINOS-POT-numu-minos-mas.dat',status='old')
        read(u,*) ((POT_numu_minos_mas(i,j), j=1,2), i=1,NBIN)
    close(u)
    open(newunit=u,file='MINOS2/MINOS-data/MINOS-POT-numu-numubar.dat',status='old')
        read(u,*) ((POT_numu_numubar(i,j), j=1,2), i=1,NBIN)
    close(u)
    open(newunit=u,file='MINOS2/MINOS-data/MINOS-upper-sigma-limits.dat',status='old')
        read(u,*) ((upper_sigma_limits(i,j), j=1,2), i=1,NBIN)
    close(u)
    open(newunit=u,file='MINOS2/MINOS-data/MINOS-sigam-stat.dat',status='old')
        read(u,*) ((sigma(i,j), j=1,2), i=1,NBIN)
    close(u)

    


    open(newunit=u,file='MINOS2/MINOS-numu-dominated-beam/MINOS-numu-dominated-beam-best-fit-oscillation.dat',status='old')
        read(u,*) ((numu_beam_best_fit(i,j), j=1,2), i=1,NBIN_beam)
    close(u)
    open(newunit=u,file='MINOS2/MINOS-numu-dominated-beam/MINOS-numu-dominated-beam-bins.dat',status='old')
        read(u,*) ((numu_beam_bins(i,j), j=1,2), i=1,NBIN_beam)
    close(u)
    open(newunit=u,file='MINOS2/MINOS-numu-dominated-beam/MINOS-numu-dominated-beam-bottom-sigma.dat',status='old')
        read(u,*) ((numu_beam_bottom_sigma(i,j), j=1,2), i=1,NBIN_beam)
    close(u)
    open(newunit=u,file='MINOS2/MINOS-numu-dominated-beam/MINOS-numu-dominated-beam-data-points.dat',status='old')
        read(u,*) ((numu_beam_data_points(i,j), j=1,2), i=1,NBIN_beam)
    close(u)
    open(newunit=u,file='MINOS2/MINOS-numu-dominated-beam/MINOS-numu-dominated-beam-NC-background.dat',status='old')
        read(u,*) ((numu_beam_NC_background(i,j), j=1,2), i=1,NBIN_beam)
    close(u)
    open(newunit=u,file='MINOS2/MINOS-numu-dominated-beam/MINOS-numu-dominated-beam-no-oscillations.dat',status='old')
        read(u,*) ((numu_beam_no_oscillations(i,j), j=1,2), i=1,NBIN_beam)
    close(u)
    open(newunit=u,file='MINOS2/MINOS-numu-dominated-beam/MINOS-numu-dominated-beam-upper-sigma.dat',status='old')
        read(u,*) ((numu_beam_upper_sigma(i,j), j=1,2), i=1,NBIN_beam)
    close(u)


    open(newunit=u,file='MINOS2/MINOS-numu-dominated-non-fiducial-mu/'//&
                      'MINOS-numu-dom-beam-non-fiducial-mu-best-fit-oscillations.dat',status='old')
        read(u,*) ((beam_nonf_mu_best_fit(i,j), j=1,2), i=1,NBIN_nonf_beam)
    close(u)
    open(newunit=u,file='MINOS2/MINOS-numu-dominated-non-fiducial-mu/'//&
                        'MINOS-numu-dom-beam-non-fiducial-mu-bins.dat',status='old')
        read(u,*) ((beam_nonf_mu_bins(i,j), j=1,2), i=1,NBIN_nonf_beam)
    close(u)
    open(newunit=u,file='MINOS2/MINOS-numu-dominated-non-fiducial-mu/'//&
                        'MINOS-numu-dom-beam-non-fiducial-mu-bottom-sigma.dat',status='old')
        read(u,*) ((beam_nonf_mu_bottom_sigma(i,j), j=1,2), i=1,NBIN_nonf_beam)
    close(u)
    open(newunit=u,file='MINOS2/MINOS-numu-dominated-non-fiducial-mu/'//&
                        'MINOS-numu-dom-beam-non-fiducial-mu-data-points.dat',status='old')
        read(u,*) ((beam_nonf_mu_data_points(i,j), j=1,2), i=1,NBIN_nonf_beam)
    close(u)
    open(newunit=u,file='MINOS2/MINOS-numu-dominated-non-fiducial-mu/'//&
                        'MINOS-numu-dom-beam-non-fiducial-mu-NC-background.dat',status='old')
        read(u,*) ((beam_nonf_mu_NC_background(i,j), j=1,2), i=1,NBIN_nonf_beam)
    close(u)
    open(newunit=u,file='MINOS2/MINOS-numu-dominated-non-fiducial-mu/'//&
                        'MINOS-numu-dom-beam-non-fiducial-mu-no-oscillations.dat',status='old')
        read(u,*) ((beam_nonf_mu_no_oscillations(i,j), j=1,2), i=1,NBIN_nonf_beam)
    close(u)
    open(newunit=u,file='MINOS2/MINOS-numu-dominated-non-fiducial-mu/'//&
                        'MINOS-numu-dom-beam-non-fiducial-mu-upper-sigma.dat',status='old')
        read(u,*) ((beam_nonf_mu_upper_sigma(i,j), j=1,2), i=1,NBIN_nonf_beam)
    close(u)


    open(newunit=u,file='MINOS2/MINOS-numubar-dominated-beam-contained-vertex/'//&
                      'MINOS-numubar-dom-beam-cont-ver-best-fit-oscillations.dat',status='old')
        read(u,*) ((beam_cv_best_fit(i,j), j=1,2), i=1,NBIN_cv_beam)
    close(u)    
    open(newunit=u,file='MINOS2/MINOS-numubar-dominated-beam-contained-vertex/'//&
                        'MINOS-numubar-dom-beam-cont-ver-bins.dat',status='old')
        read(u,*) ((beam_cv_bins(i,j), j=1,2), i=1,NBIN_cv_beam)
    close(u)
    open(newunit=u,file='MINOS2/MINOS-numubar-dominated-beam-contained-vertex/'//&
                        'MINOS-numubar-dom-beam-cont-ver-bottom-sigma.dat',status='old')
        read(u,*) ((beam_cv_bottom_sigma(i,j), j=1,2), i=1,NBIN_cv_beam)
    close(u)
    open(newunit=u,file='MINOS2/MINOS-numubar-dominated-beam-contained-vertex/'//&
                        'MINOS-numubar-dom-beam-cont-ver-data-points.dat',status='old')
        read(u,*) ((beam_cv_data_points(i,j), j=1,2), i=1,NBIN_cv_beam)
    close(u)
    open(newunit=u,file='MINOS2/MINOS-numubar-dominated-beam-contained-vertex/'//&
                        'MINOS-numubar-dom-beam-cont-ver-NC-background.dat',status='old')
        read(u,*) ((beam_cv_NC_background(i,j), j=1,2), i=1,NBIN_cv_beam)
    close(u)
    open(newunit=u,file='MINOS2/MINOS-numubar-dominated-beam-contained-vertex/'//&
                        'MINOS-numubar-dom-beam-cont-ver-no-oscillations.dat',status='old')
        read(u,*) ((beam_cv_no_oscillations(i,j), j=1,2), i=1,NBIN_cv_beam)
    close(u)
    open(newunit=u,file='MINOS2/MINOS-numubar-dominated-beam-contained-vertex/'//&
                        'MINOS-numubar-dom-beam-cont-ver-upper-sigma.dat',status='old')
        read(u,*) ((beam_cv_upper_sigma(i,j), j=1,2), i=1,NBIN_cv_beam)
    close(u)


    open(newunit=u,file='MINOS2/MINOS-numubar-enhanced-beam/'//&
                      'MINOS-numubar-enhanced-beam-best-fit-oscillations.dat',status='old')
        read(u,*) ((beam_enh_best_fit(i,j), j=1,2), i=1,NBIN_enh_beam)
    close(u)    
    open(newunit=u,file='MINOS2/MINOS-numubar-enhanced-beam/'//&
                        'MINOS-numubar-enhanced-beam-bins.dat',status='old')
        read(u,*) ((beam_enh_bins(i,j), j=1,2), i=1,NBIN_enh_beam)
    close(u)
    open(newunit=u,file='MINOS2/MINOS-numubar-enhanced-beam/'//&
                        'MINOS-numubar-enhanced-beam-bottom-sigma.dat',status='old')
        read(u,*) ((beam_enh_bottom_sigma(i,j), j=1,2), i=1,NBIN_enh_beam)
    close(u)
    open(newunit=u,file='MINOS2/MINOS-numubar-enhanced-beam/'//&
                        'MINOS-numubar-enhanced-beam-data-points.dat',status='old')
        read(u,*) ((beam_enh_data_points(i,j), j=1,2), i=1,NBIN_enh_beam)
    close(u)
    open(newunit=u,file='MINOS2/MINOS-numubar-enhanced-beam/'//&
                        'MINOS-numubar-enhanced-beam-no-oscillations.dat',status='old')
        read(u,*) ((beam_enh_no_oscillations(i,j), j=1,2), i=1,NBIN_enh_beam)
    close(u)
    open(newunit=u,file='MINOS2/MINOS-numubar-enhanced-beam/'//&
                        'MINOS-numubar-enhanced-beam-upper-sigma.dat',status='old')
        read(u,*) ((beam_enh_upper_sigma(i,j), j=1,2), i=1,NBIN_enh_beam)
    close(u)

    open(newunit=u,file='MINOS2/MINOS_uncertainties/MINOS_systematic_uncertainties.dat',status='old')
        read(u,*) sys_uncertainties
    close(u)    

end subroutine