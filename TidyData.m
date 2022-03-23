%                                     
% --.--o    |     ,--.      |         
%   |  .,---|,   .|   |,---.|--- ,---.
%   |  ||   ||   ||   |,---||    ,---|
%   `  ``---'`---|`--' `---^`---'`---^
%            `---'                    
%  Creates tidy data for a dataset, and drops to raw data section of project
%  folder setting.
addpath('~/Code/analysis')

Opt = OPTARGS;
Opt.unit = 'spikes';
Opt.taskFilter = [];
Opt.cellFilter = [];
Opt.behFilter = [];

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
spikes.beh = util.table.split(spikes.beh, "neuron");

% ===============================
% How do cells look (raster-wise) :
% Send tidy data to Julia
% ===============================
folder = coding.file.datafolder('exp_raw', 'visualize_raw_neural');

% Behavior
behavioral_properties = ["x","y","currentPathLength","currentAngle"];
beh_fn     = coding.file.datafolder('exp_raw', 'visualize_raw_neural', animal + "_" + day + "_" + 'beh.csv');
writetable(beh, beh_fn);

% Task tabble
task        = taskLib.taskTable(animal, day);
task_fn     = coding.file.datafolder('exp_raw', 'visualize_raw_neural', ...
                                         animal + "_" + day + "_" + 'task.csv');
writetable(task, task_fn);


% Spikes
labeledSpiking = units.labeledSpiking(spikes, ["block", "subblock", "blocktraj", "traj"], beh, task);
spiking_fn = coding.file.datafolder('exp_raw', 'visualize_raw_neural',...
                                        animal + "_" + day + "_" + 'labeled_spiking.csv');
writetable(labeledSpiking, spiking_fn);

% Cell properties
cell_fn    = coding.file.datafolder('exp_raw', 'visualize_raw_neural',...
                                        animal + "_" + day + "_" + 'cell.csv');
writetable(spikes.cellTable, cell_fn);

% Ripples (CA1 and cortical)
rippletime = ndb.load(animal, 'rippletime', 'inds', day);
rippletimepfc = ndb.load(animal, 'rippletimepfc', 'inds', day);
rip_fn     = coding.file.datafolder('exp_raw', 'visualize_raw_neural',...
                                        animal + "_" + day + "_" + 'ripple.csv');
rippletime    = ndb.toTidy(rippletime, 'labels', ["day","epoch"]);
rippletimepfc = ndb.toTidy(rippletimepfc, 'labels', ["day","epoch"]);
rippletime.area    = repmat("CA1", height(rippletime), 1);
rippletimepfc.area = repmat("PFC", height(rippletimepfc), 1);
writetable([rippletime; rippletimepfc], rip_fn);

% Rhythms (Just theta for now)
rhythm        = lfpLib.create.rhythmTable(animal, 'thetaref', day);
rhythm_fn     = coding.file.datafolder('exp_raw', 'visualize_raw_neural', ...
                                         animal + "_" + day + "_" + 'rhythm.nc');
util.table.netcdfwrite(rhythm, rhythm_fn);
util.notify.pushover("Finished writing rhythm table")

spikeDat = ndb.load(animal, Opt.unit, 'ind', day);
spikes   = units.getRateMatrix(animal, day,...
                            'unit', Opt.unit,...
                            'taskFilter', Opt.taskFilter,...
                            'cellFilter', Opt.cellFilter,...
                            'spikeData', spikeDat);
