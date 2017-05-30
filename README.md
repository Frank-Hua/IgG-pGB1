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
        baseline_localization.m(4)
    end

    if localization list does not exist
        finding_localization.m(5) or finding_localization_radius.m(6)
        redundant_binding_sites.m(7)
    end

    while not at the end of localization list
        analyze_localization.m(8)
    end

drift_correction_control.m(2)
    drift_correction.m(9)

baseline_localization.m(4)
    smm_intensity.m(10)
    baseline_correction.m(11)

smm_intensity.m(10)
    flip_drift_correction.m(12)
    STORM_xynm2conventional_xypixel.m(13)

finding_localization_radius.m(6)
    circle_generator(14)

analyze_localization.m(8)
    smm_intensity.m(10)
    line_generator.m(15)
    mtr_modifier.m(16)
    display_movie.m(17)
    analyze_localization2vbFRET(18)
    tr_modifier.m(19)

display_movie.m(17)
    flip_drift_correction.m(12)
    STORM_xynm2conventional_xypixel.m(13)

../utilities
    command_input.m
    plot_formatter.m