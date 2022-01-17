function spikes = analysis(animal, day, spikes, Opt)
% Notes
% -------
% Tuning to angle and distance?
% Jercog Kandel and Abott Paper
% This method simultaneously computes over data splits
% and doesn't require for-looping

% Perform sarel analyses
beh = gbehavior.lookup(animal, [], day);
[spikes, beh_filtered, ~] = units.atBehavior(beh, spikes,...
                                    'merge', true,...
                                    'query', Opt.behFilter);

% For now I'm constraining this, just so that I can complete these analysis
% faster
sarelKws = struct(...
    'beh', beh,...
    'beh_filtered', beh_filtered, ...
    'split_by',    {{["neuron", "stopWell"]}},...
    'split_by_name', "stops",...
    'useGPU', true...
    );

% ----
% MAIN
% ----
% Acquire basic distributions of variables with respect to the goals
% and computes important indices wrt those tuning curves
X = coding.sarel.main(spikes, sarelKws);
if ~isfield(spikes, 'sarel')
    spikes.sarel = X;
else
    spikes.sarel = util.struct.update(spikes.sarel, X);
end
coding.file.save(animal, day, spikes,...
    'append', true,...
    'checkpointActive', {Opt.checkpoint});
util.notify.pushover('Sarel','Finished main data');
% Compute the following indices, too
% + Directionality index :: RV
% + optimal von mises tuning curve
% + maxMean indices for non-angular varialbes

% Any of the indices that do not work, for now, are computed outside of this function

% Computes goal indices in conjuction with shuffles
spikes.sarel.ground.goalplaceindex = ...
    coding.sarel.computerGoalPlaceIndex(spikes, beh, 'ground');
spikes.sarel.stops.goalplaceindex = ...
    coding.sarel.computerGoalPlaceIndex(spikes, beh, 'stops');

% ---------
% SHUFFLES
% ---------

% If the shuffle struct already exist and not a cell, remove it
if isfield(spikes.sarel, 'shuffle') && ~iscell(spikes.sarel.shuffle)
    spikes.sarel = rmfield(spikes.sarel, 'shuffle');
end
% If field non-existant, create an empty cell
if ~isfield(spikes.sarel, 'shuffle')
    spikes.sarel.shuffle = {}; 
end

% Process each shuffle!
Shuf = {animal, day};
for iS = progress(1:Opt.nShuffle, 'Title', 'Shuffles')
    % get
    [item, Shuf] = units.shuffle.get(Shuf, iS, 'debug', false, 'shiftless', true);
    % run
    spikes.sarel.shuffle{iS} = coding.sarel.main(item, sarelKws); 
    coding.file.save(animal, day, spikes,...
        'append', true,...
        'checkpointActive', {Opt.checkpoint, iS});
end
util.notify.pushover('Sarel','Finished analyzing shuffles');

