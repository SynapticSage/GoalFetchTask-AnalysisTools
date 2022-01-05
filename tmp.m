
% ------------------
% Get the path setup
% ------------------
addpath(genpath('~/Code/analysis'))
%dbstop in units.shuffle.helper.preallocateBehaviorShifts at 11 if g==258
%dbstop in units.shuffle.conditional_time at 224

[animal, day, varargin] = deal('RY16', 36, {});

ip = inputParser;
ip.addParameter('unit', 'multiunit');                       % Spikes (curated firings) / Multiunit (uncurated firings)

ip.addParameter('shuffleStruct', []);                       % Struct of shuffling instructions (Reruns all analyses requested with sequences of shuffles specified here)
ip.addParameter('nullMethod', []);
ip.addParameter('upsample', false);                         % Upsample behavior?

% Filtration
ip.addParameter('behFilter', 'abs($velVec) > 4');           % Query to apply to all behavior
ip.addParameter('taskFilter', "string($type)==""run""");    % 
ip.addParameter('cellFilter', 'ismember(string($area), ["CA1","PFC"]) & string($tag) ~= "artifcat" & $numspikes > 150'); % 

% Analysees
ip.addParameter('analyses',   ["jercog", "sarel", "dotson"]);
ip.parse(varargin{:})
Opt = ip.Results;

%shuffleStructDefault

% -------------
% Get Unit data
% -------------
spikeDat = ndb.load(animal, Opt.unit, 'ind', day);
spikes   = units.getRateMatrix(animal, day,...
    'unit', Opt.unit,...
    'taskFilter', Opt.taskFilter,...
    'ploton',false,...
    'dense', false,...
    'cellFilter', Opt.cellFilter,...
    'spikeData', spikeDat);


beh = gbehavior.lookup(animal, [], day);
beh.x = beh.pos(:,1);
beh.y = beh.pos(:,2);
props = ["x", "y"];
beh = util.table.castefficient(beh,...
    'compressReals', true, 'negativeNan', true);
shift = -1:0.1:1;


% Prepare shuffle
kws = struct('groups', [], 'props', props, 'shift', shift,...
    'cacheToDisk', {{animal, day}},...
    'nShuffle', Opt.shuffleOpt.nShuffle,...
    'preallocationSize', 4);
groupby = ["epoch", "period"]; % properties in behavior table to shuffle within
profile on;
units.shuffle.conditional_time(beh, spikes, groupby, kws);
