function FP = shuffcorrect(FP, varargin)
% Applies different corrections to mutual information using information from
% the shuffle
%
% ------
% Inputs
% ------
%
% FP : struct
%
% futurepast struct
%
% FP.MI
% FP.shuffle
%
% --------------------
% Output modifications
% --------------------
%
% FP.MIcorrected
% FP.MIcorrectedFistribution
%
% Reference: 
% [[zotero://select/library/items/YT5PGHGJ][Dotson/Yartsev Supplement]]
%
% Uses non-spatial cell shuffle to obtain a curve that indicates
% avgMI = f( firing_rate, tau )
%
% And then we use this fitted curve to subtract/correct/normalize information
% changes due to shuffling/FR interaction.

ip = inputParser;
ip.addParameter('normalization', 'divide');
ip.addParameter('shufflereduce', @median);
ip.addParameter('frBinSize', []); % number of samples per bin
% --------- Statistics -------------
ip.addParameter('sigThreshold', 0.05); % number of samples per bin
ip.addParameter('bonferonni', true);
ip.parse(varargin{:});
Opt = ip.Results;

if iscell(FP.shuffle)
    FP.shuffle= cat(1,FP.shuffle{:});
end

% -------------------------------------
% Prepare the variables we will work on
% -------------------------------------

% Actual mutual information
mi.main = FP.MI.entropy;

% Shuffle mutual information
for i= 1:numel(FP.shuffle)
    mi.shuffle{i} = FP.shuffle(i).MI.entropy;
end
%MIshuffle = Opt.shufflereduce(MIshuffle, 1);

% And we will be needing firing rates for these time shifts per neuron
fr = struct();
fr.main = FP.spikeCount./FP.spikeTime;
receptiveFieldDims = [3:ndims(fr.main)];
fr.main = mean(fr.main, receptiveFieldDims); % mean out non-neuron and non-tau dims
spCountShuff = nd.fieldGetCell(FP.shuffle, 'spikeCount', 'shiftdim', -1);
spTimeShuff = nd.fieldGetCell(FP.shuffle, 'spikeTime', 'shiftdim', -1);
spCountShuff = cat(1, spCountShuff{:});
spTimeShuff  = cat(1, spTimeShuff{:});
receptiveFieldDims = 4:ndims(spCountShuff);
fr.shuff = squeeze(mean(spCountShuff./spTimeShuff, receptiveFieldDims));
% now we have FR for main data, and FR for shuffle

% Bin firing rates
% ----------------
if isempty(Opt.frBinSize)
    [fr.N, fr.edges] = histcounts(fr.shuff);
else
    [fr.N, fr.edges] = histcounts(fr.shuff, Opt.frBinSize);
end

keyboard
% Lookup structure
% ----------------
% logical and list
fr.lookup = struct();
for tau = 1:numel(spikes.dotson.shift)
for f   = 1:numel(fr_N)
    fr_edge = fr.edges(f, f+1);
    fr.lookup(f,tau).fr_edge = fr_edge;
    fr.lookup(f,tau).logical = util.logical.minmax(fr(tau,:,:,:,:,:), fr_edge);
end
end

% ---------------------------------------------
% Correction strategy 1: corr1 = MI( FR, tau ) : discover this average function
% ---------------------------------------------
for f = 1:numel(fr.N)
for tau = 1:numel(FP.shift)
    inds =fr.lookup(f, tau); %inds, might have to repmat this
    correct.shuff_frShift(f, tau) = (mean(mi.shuffle(inds),ndims(mi.shuffle))...
        * FP.spikeTime)/sum(FP.spikeTime,'all'); 
    % expectation( mi(i) * p(i) )
end
end

% ---------------------------------------------
% Correction strategy 2: corr2 = MI( neuron, tau, *receptive bins) : discover this average function
% ---------------------------------------------
correct.shuff_neuronShiftRF = mean(mi.shuffle, ndims(mi.shuffle));


% Which are our neurons with field information?
% ----------------------------------------------
% I believe this is neurons with where significance is determined by the MI
% versus the shuffle MI for that neuron
%
% Now, one discrepency. In the supplment, they caim 164 insignificant spatial
% info (SI) fields. In the text, 90% of the neurons had high SI, 183/204.
% Doesn't that leave only 21 cells insignificant?

% Obtain the average mutual information for all of the lookup vals
% ----------------------------------------------------------------
for tau = 1:numel(FP.shift)
for f   = 1:numel(fr.N)
    significance = mean(mi.main > mi.shuffle, ndims(mi.shuffle));
    neuron_dim = 2;
    sigunits = any(significance, setdiff(1:ndims(significance), neuron_dim));
end
end


% ---------------------------------------------
% Correction strategy 1b: corr1 = MI( FR, tau ) : discover this average function
% ---------------------------------------------
for f = 1:numel(fr.N)
for tau = 1:numel(FP.shift)
    inds =fr.lookup(f, tau) & sigunits; %inds, might have to repmat this
    correct.insig_frShift(f, tau) = (mean(mi.main(inds),ndims(mi.main))...
        * FP.spikeTime)/sum(FP.spikeTime,'all'); 
    % expectation( mi(i) * p(i) )
end
end


% -----------------
% Apply corrections
% -----------------
for correctionMethod = string(fieldnames(correct))'
    switch lower(char(Opt.normalization))
    case 'divide'
        FP.(correctionMethod).MI        = mi.main ./ correct.(correctionMethod);
        FP.(correctionMethod).shuffDist = mi.shuffle ./ correct.(correctionMethod);
    case 'subtract'
        FP.(correctionMethod).MI        = mi.main - correct.(correctionMethod);
        FP.(correctionMethod).shuffDist = mi.shuffle - correct.(correctionMethod);
    end
end
% ------------------


