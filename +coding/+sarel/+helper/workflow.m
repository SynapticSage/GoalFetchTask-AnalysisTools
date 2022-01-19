function out = workflow(spikes, fields, Binning, Opt)
    % This workflow repeats, but marginalized in different ways
    %

tic
%% ------------- %%
%% Tuning Curves %%
%% ------------- %%
% For every combination of fields, histcount the major variables
out = coding.sarel.tuning(spikes, fields, Binning, Opt);
out.Dimensions = [fields, "bins"];
out = coding.sarel.occupancyNormalize(out, Opt.occupancy);


%% ------- %%
%% METRICS %%
%% ------- %%

% occ norm based copmutations
out.rayleigh = coding.sarel.metric.directionalityIndex(out.occNorm, Binning.angleCenters, -1,  "Angle");
out.vm = coding.sarel.metric.fitVonMises(out.occNorm, Binning.angleCenters, -1,  "Angle");
out.maxmean_indices = coding.sarel.metric.maxMeanIndices(out,    Binning);
%out.goalplaceindex  = coding.sarel.metric.GoalPlaceIndex(spikes, fields, Binning);

fprintf("\nProcessing %s required %2.0f seconds\n", join(fields, " "), toc);
