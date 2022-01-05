function comptueGoalPlaceIndex(spikes, beh,...
         marginal_name, varargin)
%
%  

ip = inputParser;
ip.parse(varargin{:})
Opt = ip.Results;

assert(isfield(spikes.sarel, marginal_name));
msarel = spikes.sarel.(marginal_name);
Binning = spikes.sarel.Binning;

scaff = Binning.Scaffold;

% Add information about the spikes bins
[spikes, Ssummary, Opt] = coding.field.scaffold2binning(spikes, scaff, Opt);


% Classify GOAL cell
% ------------------
% MEASUREMENT RayleighVector
% (1) dir index (rayleigh vector) significant based on shuffles: for eeach
% neuron, concatonatee all in spike sequences for trials and time-shift spikes
% by random interval, circularly. directonal tuned if greater than 95%
% percentile of shifts
signal  = msarel.rayleigh.currentAngle.directionality_index;
shuffle = spikes.shuffle;
% (2)) dir index > 0.25
% (3) more than 50 spikes
% (4) pearson corr > 0.5 between tuningg for even and odd minutes
% (5) stronger goal than place index

% Classify PLACE cell
% ------------------
% (1) sig spatial information on shuffling above

