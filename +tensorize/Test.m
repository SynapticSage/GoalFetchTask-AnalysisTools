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

[X, T, data] = tensorize.rateBy(spikes, beh,...
    'vars_timeChanging', "traj",...
    'vars_category', "correct",...
    'vars_quantile', "speed"...
    );


reorderDims = 3; % trial
stackDims = 1:2; % speed x correct
nonstackDims = 4:5; % 4neuron x time
[Xreorder, run_stats] = tensorize.equalSamplesDims(X, reorderDims, stackDims, nonstackDims);


%% Now we wanna actually do some science
vars_timeChanging = "traj";
vars_category = ["correct", "cuemem"];
vars_quantile = [];

beh = gbehavior.lookup(animal, [], day);
props = ["time", vars_timeChanging, vars_category, vars_quantile];
beh = beh(:,props);
beh = util.table.castefficient(beh,...
    'compressReals', true,...
    'negativeNan', true,...
    'compressRealsExcept', {'time'});

[X, T, data{1}] = tensorize.rateBy(spikes, beh,...
    'vars_timeChanging', vars_timeChanging,...
    'vars_category', vars_category,...
    'vars_quantile', vars_quantile...
);

reorderDims = 3; % trial
stackDims = 1:2; % correct x cuemem
nonstackDims = 4:5; % neuron x time
[Xreorder, run_stats] = tensorize.equalSamplesDims(X, reorderDims, stackDims, nonstackDims);

Xhpc = X(:,:,:,:,spikes.cellTable.area=="CA1");
Xpfc = X(:,:,:,:,spikes.cellTable.area=="PFC",:);
