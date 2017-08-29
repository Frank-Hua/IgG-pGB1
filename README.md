# IgG-pGB1
MatLab codes for the IgG-pGB1 project

Code structure:

single_flow_analysis.m(1)
    if do drift correction
        drift_correction_control.m(2)
    end

    if count spot number per frame
        spot_number_per_frame.m(3)
    end

    if correct baseline
        baseline_site.m(4)
    end

    if site list does not exist
        finding_site_radius.m(5)
        redundant_binding_sites.m(6)
    end

    while not at the end of site list
        analyze_site.m(7)
    end

drift_correction_control.m(2)
    drift_correction.m(8)

baseline_site.m(4)
    smm_intensity.m(9)
    baseline_correction.m(10)

smm_intensity.m(9)
    flip_drift_correction.m(11)
    STORM_xynm2conventional_xypixel.m(12)

finding_site_radius.m(5)
    circle_generator(13)

analyze_site.m(7)
    smm_intensity.m(9)
    line_generator.m(14)
    mtr_modifier.m(15)
    display_movie.m(16)
    tr_modifier.m(17)

display_movie.m(16)
    flip_drift_correction.m(11)
    STORM_xynm2conventional_xypixel.m(12)

../utilities
    command_input.m(18)
    plot_formatter.m(19)