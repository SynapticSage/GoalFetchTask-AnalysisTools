function [spikeOut, behOut, scaff] = calc(spikes, varargin)
% function field = calcfield(beh, props)
%
% Generalizes calculating a field from behavior properties props
%
% Inputs
% ------
%
% Returns
% -------
% Fields of a behavior over spiking per neuron and overall behaviooutr

ip = inputParser;
ip.KeepUnmatched = true;
ip.addParameter('useGPU', false);
ip.addParameter('useGPU_binning', true); % if useGPU && useGPU_binning
ip.addParameter('beh', []);
ip.addParameter('grid', 25);
ip.addParameter('scaff', []);
ip.addParameter('props', []);
ip.addParameter('blocksize', -1);
ip.addParameter('samplingPeriod', []);
ip.addParameter('suppressProgress', false);
ip.parse(varargin{:})
Opt = ip.Results;
Opt.props = string(Opt.props);
%% -------------------------
%% CREATE SHORTCUT VARIABLES
%% -------------------------
if isempty(Opt.beh)
    B = spikes.beh;
else
    B = Opt.beh;
end
if isstruct(spikes)
    S = spikes.beh;
elseif istable(spikes)
    S = spikes;
else
    error("unimplemented");
end
Opt.samplingPeriod = median(diff(B.time));
if iscell(S)
    %uNeuron = 1:numel(S);
    error("Rest of the code for cell mode is unimplemented")
else
    if ismember(fieldnames(S), 'uNeuron')
        uNeuron = spikes.uNeuron;
    else
        uNeuron = unique(S.neuron);
    end
end
clear spikes
Opt.beh = [];

% -------------
% Bring to GPU?
% -------------
if Opt.useGPU
    B = util.table.table2GPUtable(B);
    S = util.table.table2GPUtable(S);
    kws = {'gpuArray'};
    uNeuron = gpuArray(uNeuron);
else
    kws = {};
end

% ----------
% Scaffold
% ----------
scaff = coding.field.scaffold(B, Opt);
% Store overall center counts : The field grid!
nCenters     = arrayfun(@(x) numel(x.center), scaff);
nCentersCell = num2cell(nCenters);

% ----------
% Bin spikes
% ----------
[S, Ssummary, Opt] = coding.field.scaffold2binning(S, scaff, Opt);

% ------------
% Bin behavior
% ------------
[B, Bsummary, Opt] = coding.field.scaffold2binning(B, scaff, Opt);
assert(~all(Bsummary.visit_time == 0,'all'), 'Behavior times for all bins are zero')


% Prepare to absorb results
spikeOut = struct(...
    'spikeCount', nan(max(uNeuron), nCentersCell{:}, kws{:}),...
    'spikeTime',  nan(max(uNeuron), nCentersCell{:}, kws{:})...
);

if Opt.suppressProgress
    N = uNeuron';
else
    N = progress(uNeuron', 'Title' , 'Calculating fields' );
end
neurons = gpuArray(S.neuron);
overall_bin  = util.type.castefficient(S.overall_bin);
%overall_bin  = S.overall_bin;
uOverallBins = unique(overall_bin);
uOverallBins = uOverallBins(~isnan(uOverallBins) & uOverallBins > 0);
for neuron = N 

    %Subset the bin table and the spike table - they are 1 to 1
    % Which entries
    nfilt  = find(neurons == neuron); % somehow for gpu(table), find is 1/3 the time
    subset = overall_bin(nfilt);

    % Get the distribution of times the neuron falls into all possible explicit
    % bin combinations
    [N, E] = histcounts(subset, uOverallBins);
    E = double(E);
    E = E(1:end-1);
    clean = E>0;
    N = N(clean);
    E = E(clean);
    assert(util.num.fractionQueryMatch(N) >  0);
    spikeOut.spikeCount(neuron, E) = N;
    spikeOut.spikeTime(neuron, E) =  single(N) * Opt.samplingPeriod; %TODO need to subtract times in the same bin
    
    % One way to do the TODO above, is this
    % Bin specific calculations
    %for bin = unique(S.overall_bin)'
    %    filt = subset.S_xy_bin == bin;
    %    uTimes = unique(spikes.beh.time(filt));
    %    spikeOut.spikeTime(neuron, bin) = sum(numel(uTimes)) * samplingPeriod; % Sum of sampling periods that have spikes
    %end
    
end

% The only TRUE zero times should be those where the animal visited
% the location, but no spikes were emitted (this paragraph)
behtime_nonzero = repmat(shiftdim(Bsummary.visit_time,-1), [size(spikeOut.spikeTime,1),1,1,1]);
behtime_nonzero = behtime_nonzero ~= 0;
filt = isnan(spikeOut.spikeTime) & behtime_nonzero;
spikeOut.spikeTime(filt) = 0;


% Occupancy normalize, dividing time spent visiting by the spikeCount
spikeOut.FR_occNorm = bsxfun(@rdivide, ...
    spikeOut.spikeCount, shiftdim(Bsummary.visit_time, -1));
spikeOut.multiunit = Ssummary;
behOut             = Bsummary;

if Opt.useGPU
    spikeOut = util.struct.GPUstruct2struct(spikeOut);
    behOut   = util.struct.GPUstruct2struct(behOut);
end
%------------------- HELPER FUNCTIONS -------------------------------------
