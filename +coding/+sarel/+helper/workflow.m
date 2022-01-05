function out = workflow(spikes, fields, Binning, Opt)
    % This workflow repeats, but marginalized in different ways

    tic
    % For every combination of fields, histcount the major variables
    out = coding.sarel.tuning(spikes, fields, Binning, Opt);
    out = coding.sarel.occupancyNormalize(out, Opt.occupancy);

    % occ norm based copmutations
    out.rayleigh = coding.sarel.computeDirectionalityIndex(out.occNorm, Binning.angleCenters, -1,  "Angle");
    out.vm = coding.sarel.computeVonMises(out.occNorm, Binning.angleCenters, -1,  "Angle");

    % coding indices
    out.maxmean_indices = coding.sarel.computeMaxMeanIndices(out,    Binning);
    out.goalplaceindex  = coding.sarel.computerGoalPlaceIndex(spikes, fields, Binning);

    fprintf("\nProcessing %s required %2.0f seconds\n", join(fields, " "), toc);
