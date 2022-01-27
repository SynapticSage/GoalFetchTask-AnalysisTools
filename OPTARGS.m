function Opt = gfmOptargs(varargin)
% Main optional argumeent structure for my overaching analyses

ip = inputParser;

% Characteristics of spikes/behavior data
ip.addParameter('unit', 'multiunit'); % spikes (curated firings) / Multiunit (uncurated firings)
ip.addParameter('upsample', false);   % upsample behavior?

% Filtration
ip.addParameter('behFilter', 'abs($velVec) > 4');           % query to apply to all behavior
ip.addParameter('taskFilter', "string($type)==""run""");    % query epochs
ip.addParameter('cellFilter', 'ismember(string($area), ["CA1","PFC"]) & string($tag) ~= "artifcat" & $numspikes > 150'); % 

% Shuffle and null distributions
ip.addParameter('shuffleStruct', [])   % struct of shuffling instructions (Reruns all analyses requested with sequences of shuffles specified here)
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

% Shuffle struct specific modifications
if isempty(Opt.shuffleStruct)
    Opt.shuffleStruct = units.shuffle.inputParser();
    Opt.shuffleStruct.parse(varargin{:});
    Opt.shuffleStruct = Opt.shuffleStruct.Results;
end
% Shortcuts
if isempty(Opt.nShuffle) 
     Opt.nShuffle = Opt.shuffleStruct.nShuffle;
end
if isempty(Opt.skipShuffled) 
    Opt.skipShuffled = Opt.shuffleStruct.skipShuffled;
end

% Convert Analyses to string, if not already
Opt.analyses = string(Opt.analyses);



