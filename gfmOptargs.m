function Opt = gfmOptargs(varargin)

ip = inputParser;

% Characteristics of spikes/behavior data
ip.addParameter('unit', 'multiunit'); % spikes (curated firings) / Multiunit (uncurated firings)
ip.addParameter('upsample', false);   % upsample behavior?

% Filtration
ip.addParameter('behFilter', 'abs($velVec) > 4');           % query to apply to all behavior
ip.addParameter('taskFilter', "string($type)==""run""");    % query epochs
ip.addParameter('cellFilter', 'ismember(string($area), ["CA1","PFC"]) & string($tag) ~= "artifcat" & $numspikes > 150'); % 

% Shuffle and null distributions
ip.addParameter('shuffleStruct', units.shuffle.optargs())   % struct of shuffling instructions (Reruns all analyses requested with sequences of shuffles specified here)
ip.addParameter('nullMethod', []);
% Shortcut varargin into shuffleStruct options
ip.addParameter('nShuffle', []);     % number of shuffles to generate
ip.addParameter('skipShuffled', []); % whether to skip already computed shuffles

% Analyses
ip.addParameter('analyses',   ["jercog", "sarel", "dotson"]);
ip.addParameter('shift', -1:0.1:1);

% Checkpointing
ip.addParameter('checkpoint', 10); % checkpoints if > 1, where number for loops
                                   % is how often,
                                   % if logical, skip checkpoints in loops, but
                                   % checkpoint everywhere else

ip.parse(varargin{:})
Opt = ip.Results;
if isempty(Opt.nShuffle) && ~isempty(Opt.shuffleStruct)
     Opt.nShuffle = Opt.shuffleStruct.nShuffle;
end
if isempty(Opt.skipShuffled) && ~isempty(Opt.shuffleStruct)
    Opt.skipShuffled = Opt.shuffleStruct.skipShuffled;
end
Opt.analyses = string(Opt.analyses);

