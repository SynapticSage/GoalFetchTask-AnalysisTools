function F = shuffcorrect(F, varargin)
% Applies different corrections to mutual information using information from
% the shuffles
%
% ------
% Inputs
% ------
%
% F : struct
%
% futurepast struct
%
% F.MI
% F.shuffles
%
% --------------------
% Output modifications
% --------------------
%
% F.MIcorrected
% F.MIcorrectedFistribution
%
% Reference: 
% [[zotero://select/library/items/YT5PGHGJ][Dotson/Yartsev Supplement]]
%
% Uses non-spatial cell shuffles to obtain a curve that indicates
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
ip.addParameter('bonferonni', true):
ip.parse(varargin{:});
Opt = ip.Results;

% -------------------------------------
% Prepare the variables we will work on
% -------------------------------------

% Actual mutual information
MI = F.MI;

% Shuffle mutual information
MIshuffles = cellfun(@(x) x.MI, F.shuffles, 'UniformOutput', false);
MIshuffles = cat(1, MIshuffles{:});
%MIshuffles = Opt.shufflereduce(MIshuffles, 1);

% And we will be needing firing rates for these time shifts per neuron
fr = struct();
fr.main = F.spikeCount./F.spikeTime;
fr.main = mean(fr.main, [3:ndims(fr.main)]); % mean out non-neuron and non-tau dims
spCountShuff = nd.fieldGetCell(F.shuffles, 'spikeCount');
spTimeShuff = nd.fieldGetCell(F.shuffles, 'spikeTime');
spCountShuff = cat(ndims(fr)+1, spCountShuff);
spTimeShuff  = cat(ndims(fr)+1, spTimeShuff);
fr.shuff = mean(spCountShuff./spTimeShuff, [3:ndims(spCountShuff)]);
% now we have FR for main data, and FR for shuffles

% ---------------------------------------------
% F( FR, tau ) : discover this average function
% ---------------------------------------------

% Which are our neurons with field information?
% ----------------------------------------------
% I believe this is neurons with where significance is determined by the MI
% versus the shuffle MI for that neuron
%
% Now, one discrepency. In the supplment, they caim 164 insignificant spatial
% info (SI) fields. In the text, 90% of the neurons had high SI, 183/204.
% Doesn't that leave only 21 cells insignificant?



% Bin firing rates
% ----------------
if isempty(Opt.frBinSize)
    fr_N, fr_edges = histcount(frShuff);
else
    fr_N, fr_edges = histcount(frShuff, Opt.frBinSize);
end

% Lookup structure
% ----------------
% logical and list
fr.lookup = struct();
for tau = 1:numel(spikes.dotson.shift)
for f   = 1:numel(fr_N)
    fr_edge = fr_edges(f, f+1);
    fr.lookup(f,tau).fr_edge = fr_edge;
    fr.lookup(f,tau).logical = util.logical.minmax(fr(tau,:,:,:,:,:), fr_edge);
end
end

% Obtain the average mutual information for all of the lookup vals
% ----------------------------------------------------------------
for tau = 1:numel(spikes.dotson.shift)
for f   = 1:numel(fr_N)
    
end
end


% ---------------------------------------------
% Broad shuffle correction (not from the paper)
% ---------------------------------------------
switch lower(char(Opt.normalization))
case 'divide'
    broadCorrect.MI     = MI ./ Opt.shufflereduce(MIshuffles);
    broadCorrect.MIdist = MI ./ MIshuffles;
case 'subtract'
    broadCorrect.MI     = MI - Opt.shufflereduce(MIshuffles);
    broadCorrect.MIdist = MI - MIshuffles;
end
% ---------------------------------------------


% ---------------------------------
% Yartsev Dotson Shuffle Correction 
% ---------------------------------
switch lower(char(Opt.normalization))
case 'divide'
case 'subtract'
end
% ---------------------------------
