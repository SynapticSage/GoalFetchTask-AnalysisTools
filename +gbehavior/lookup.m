function [behavior, times_to_throwout] = lookup(animal, times, inds, varargin)
% Emulates the behavior.lookup function I wrote for wtrack
% data, execpt with the goalmaze
%
% TODO fix so not require dropfields

ip = inputParser;

% Shape of the data
ip.addParameter('valueOnly',false);    % Whether to not return ccolumns that function as indices (time, traj, day, epoch)
ip.addParameter('scaleVars',false);    % Whether to scale vars

% Data constraints
ip.addParameter('throwOutSleep',true); % Whether to scale vars
ip.addParameter('throwOutThresh',2);   % If 2 seconds from nearest point, toss it
ip.addParameter('dropFields',["goalVec", "egoVec", "angle","headEgoAngle", "euclidianDistance"]);   % If 2 seconds from nearest point, toss it

% USE PRELOADED inputs?
ip.addParameter('behavior', []);       % Preloaded behavior file
ip.addParameter('trajTable', []);       % Preloaded behavior file
ip.parse(varargin{:})
Opt = ip.Results;

if nargin < 2
    times = [];
end
if nargin < 3
    inds = [];
end


if Opt.throwOutSleep
    % Section being developed
    task    = ndb.load(animal,'task', 'inds', inds);
    tmp  = ndb.indicesMatrixForm(task);
    if ~isempty(inds)
        inds = tmp(ismember(tmp(:, 1:size(inds,2)), inds), :);
    else
        inds = tmp;
    end
    run_inds = evaluatefilter(task, 'string($type) == "run"');
    if ~isempty(inds)
        inds = intersect(inds, run_inds, 'rows');
    else
        inds = run_inds;
    end
end


% ==================
% Ego struct entries
% ==================

if ~isempty(Opt.behavior)
    behavior = Opt.behavior;
else
    behavior = ndb.load(animal, 'egocentric', 'inds', inds);
end

behavior_inds = ndb.indicesMatrixForm(behavior);
behavior      = ndb.toNd(behavior);

B = {};
for binds = inds'
    X = rmfield(nd.get(behavior, binds), ["gridTable", "pathX", "pathY", "trialTimes", Opt.dropFields]);
    B{end+1} = struct2table(X);
    B{end}.day = repmat(binds(1), height(B{end}), 1);
    B{end}.epoch = repmat(binds(2), height(B{end}), 1);
end

B = B(~cellfun(@isempty, B));
behavior         = cat(1, B{:});
behavior.time    = behavior.postime;
behavior.postime = [];

% Separate x and y
behavior.x = behavior.pos(:,1);
behavior.y = behavior.pos(:,2);


%% ===================
%% Traj Table Entriees
%% ===================

if ~isempty(Opt.trajTable)
    trajTable = Opt.trajTable;
else
    trajTable = ndb.load(animal, 'traj', 'inds', inds);
end

traj_inds = ndb.indicesMatrixForm(trajTable);
behavior.block = nan(height(behavior),1);
behavior.subblock = nan(height(behavior),1);
behavior.blocktraj = nan(height(behavior),1);
behavior.traj = nan(height(behavior),1);
for ind = traj_inds'

    epoch_table = ndb.get(trajTable, ind);
    for row = 1:height(epoch_table)

        row_table = epoch_table(row,:);
        start = row_table.start;
        stop = row_table.stop;

        matches = util.constrain.minmax(behavior.time, [start, stop]);
        if sum(matches) == 0
            continue
        end

        behavior.block(matches) = int16(row_table.block);       % block id : cue and memory
        behavior.subblock(matches) = int16(row_table.subblock); % subblock within a block : cue or memory
        behavior.blocktraj(matches) = int16(row_table.traj);    % traj within a block : there and back
        behavior.traj(matches) = int16(row_table.index);        % overall traj index : there or backjkk
    end
end




%% ============================
%% Derived calculation entriees
%% ============================

% behavior.periods : traj and inactive period index
traj = behavior.traj;
traj(isnan(traj)) = 0;
dTraj = diff(traj);
changes = [0; dTraj] ~= 0;
behavior.period = cumsum(changes);

% NOTE eventually, we might want to expand trial times to time length to have start stop

% NOTE missing start stop, superblock-trial and block-trials

% Acceleration

% Calculate percentage into each traj

%  idPhi

% Trajectory-correctness n-back and n-forward


if isempty(times)
    inds = true(size(behavior.time));
else
    inds = interp1(behavior.time, 1:height(behavior), times, 'nearest');
end
times_to_throwout = isnan(inds);

if Opt.throwOutThresh > 0
    disp("Throwing out times greater than " + Opt.throwOutThresh + " seconds away");
    chunksize = 5000;
    for t = progress(1:chunksize:numel(times),'Title','Throwing out bad times')
        stop = min(t+chunksize-1,numel(times));
        tinds = t:stop;
        time = times(tinds);
        toss_inds = all( abs( behavior.time - time(:)' ) > Opt.throwOutThresh,  1);
        times_to_throwout(tinds(toss_inds)) = true;
    end
end

% Return behavior
inds(times_to_throwout) = [];
behavior = behavior(inds,:);

% Set variables on same scale?
if Opt.scaleVars
    for field = string(fieldnames(behavior))'
        try
            behavior.(field) = behavior.(field) ./ sqrt(nanvar(behavior.(field)));
        catch
        end
    end
end

behavior = gbehavior.clean.removeduplicatetimes(behavior);
