% Test tensorization

[animal, day, varargin] = deal('RY16', 36, {});
false
ip = inputParser;
ip.addParameter('unit', 'multiunit');                       % Spikes (curated firings) / Multiunit (uncurated firings)

ip.addParameter('shuffleStruct', []);                       % Struct of shuffling instructions (Reruns all analyses requested with sequences of shuffles specified here)
ip.addParameter('nullMethod', []);
ip.addParameter('nShuffle', 50);                            % Number of shuffles to generate
ip.addParameter('upsample', false);                         % Upsample behavior?

% Filtration
ip.addParameter('behFilter', 'abs($velVec) > 4');           % Query to apply to all behavior
ip.addParameter('taskFilter', "string($type)==""run""");    % 
ip.addParameter('cellFilter', 'ismember(string($area), ["CA1","PFC"]) & string($tag) ~= "artifcat" & $numspikes > 150'); % 

% Analyses
ip.addParameter('analyses',   ["jercog", "sarel", "dotson"]);
ip.parse(varargin{:})
Opt = ip.Results;

% ------------------
% Get the path setup
% ------------------
addpath(genpath('~/Code/analysis'))

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

spikes = units.sparseToDense(spikes,...
            'startTime', spikes.timePeriods(1,1), ...
            'endTime', spikes.timePeriods(1,2), ...
            'interpMethod', 'linear',...
            'method', 'instantaneous');

beh = gbehavior.lookup(animal, [], day);
props = ["time","traj","correct","speed"];
beh = beh(:,props);
beh = util.table.castefficient(beh,...
    'compressReals', true,...
    'negativeNan', true,...
    'compressRealsExcept', {'time'});

[X, T, run_stats] = tensorize.rateBy(spikes, beh,...
    'vars_timeChanging',['traj'],...
    'vars_category',['correct'],...
    'vars_quantile',['speed']...
    );
