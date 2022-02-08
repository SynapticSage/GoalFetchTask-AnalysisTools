addpath('~/Code/analysis')

Opt = OPTARGS;
Opt.unit = 'spikes';
Opt.taskFilter = [];
Opt.cellFilter = [];

[animal, day] = deal('RY16', 36);

% ==========
% Input data
% ==========

% -------------
% Get Unit data
% -------------
spikeDat = ndb.load(animal, Opt.unit, 'ind', day);
spikes   = units.getRateMatrix(animal, day,...
                            'unit', Opt.unit,...
                            'taskFilter', Opt.taskFilter,...
                            'cellFilter', Opt.cellFilter,...
                            'spikeData', spikeDat);
clear spikeDat

beh = gbehavior.lookup(animal, [], day);
[spikes, beh_filtered, ~] = units.atBehavior(beh, spikes,...
    'useGPU', false,....
    'merge',  false,...
    'returnIndices', true, ...
    'query', Opt.behFilter);

% =========
% How do ceclls fire in reseponse to these?
% =========

% -----
% Combos
% -----
prop.xy = struct('prop', ["x", "y"], ...
                 'name', 'xy',...
                 'description', 'Place field');
prop.xygoalvec = struct('prop', ["x", "y", "currentAngle", "currentDistance"], ...
                        'name', 'xy',...
                        'description', 'Place field');
prop.xyheadvec = struct('prop', ["x", "y", "headdir"],...
                        'name', 'xyheadvec',  ...
                        'description', 'head place field');
prop.goalvec = struct('prop', ["currentAngle", "currentDistance"], ...
                      'name', 'goalvec',...
                      'description', 'goal vec field');

% ----------
% Get fields
% ----------

% -----------------------
% Get tidy firing  & save
% -----------------------


