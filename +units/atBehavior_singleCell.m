function [spikesBeh, spikeTimes] = ...
        atBehavior_singleCell(spikeTimes, beh, varargin)
% Handles the single cell time lookup into the behavior (either the behavior
% properties themselves or the indices of the behavior table that match). It
% can also handle time shifting those lookups efficiently.

% --------------
% Prepare inputs
% --------------
ip = inputParser;
ip.KeepUnmatched = true;
ip.addParameter('matchTimeTolerance', 1/27);      % tolerance for a spike time to match a behavior
ip.addParameter('violationHandling', []);
ip.addParameter('useGPU', false);
ip.addParameter('shift', 0);
ip.addParameter('prop', []);
ip.addParameter('returnIndices', false); % This triggers a behavior where we return a vector of indices instead of a table.
                                            % if we have shifts, then it's a matrix.
ip.addParameter('output', 'table'); % {table} | matrix ... ultimately determinse either straightup what we return, eeither a ceell of this type of the type itself
ip.addParameter('concatShifts', []); %wehther to concatonate the shifts we get at the end of the computation
ip.addParameter('annotateNeuron', []);
ip.addParameter('annotateShift', true);
ip.parse(varargin{:})
Opt = ip.Results;
Opt.violationHandling = lower(Opt.violationHandling);

if isempty(Opt.violationHandling)
    if Opt.returnIndices
        Opt.violationHandling = 'nan';
    else
        Opt.violationHandling = 'remove';
    end
end
if isempty(Opt.concatShifts)
    if Opt.returnIndices
        Opt.concatShifts = true;
    else
        Opt.concatShifts = false;
    end
end
if Opt.returnIndices && strcmp(Opt.violationHandling, 'remove')
    error('incomplatible modes violationHandling="remove" and Opt.returnIndices=true');
end

spikeTimesExist = ~isempty(spikeTimes);
multipleShifts  = length(Opt.shift) > 1;

% --------------------------
% Main : when we have spikes
% --------------------------
if spikeTimesExist

    % BRING TO GPU
    % ------------
    if Opt.useGPU
        spikeTimes = gpuArray(spikeTimes);
        beh        = util.table.table2GPUtable(beh);
    end

    % INTERPOLATE
    % -----------
    % Warning: sometimes using single here can lead to non-unique times
    try
        inds = interp1(beh.time, 1:numel(beh.time), ...
            spikeTimes(:) + Opt.shift(:)', ...
            'nearest', 'extrap');
    catch ME
        if any(util.getduplicates_logical(beh.time))
            error("Your beh.time entries have duplicate values");
        end
    end

    % MUlTIPLE SHIFTS 
    % ---------------
    if multipleShifts
        % JUST INDICES
        if Opt.returnIndices
            spikesBeh = arrayfun(...
                @(col) inds(:, col),...
                1:numel(Opt.shift),...
                'UniformOutput', false);
        % BEHAVIOR OF INDICES
        else
            spikesBeh = arrayfun(...
                @(col) beh(inds(:, col), :),...
                1:numel(Opt.shift),...
                'UniformOutput', false);
        end
        % Have to index and reshape, otherwise, always get raveled output for table column
        behavior_times = beh.time(inds(:));
        behavior_times = reshape(behavior_times, size(inds)); 
        violations = abs(behavior_times - (spikeTimes(:)+Opt.shift(:)')) > Opt.matchTimeTolerance; 
        for shift = 1:numel(Opt.shift)
            switch Opt.violationHandling
                case 'remove'
                spikesBeh{shift} = spikesBeh{shift}(~violations(:,shift),:);
                case 'nan'
                spikesBeh{shift}(violations(:,shift),:) = nan;
            end
        end
    % ONE SHIFT
    % ---------------
    else
        % JUST INDICES
        if Opt.returnIndices
            spikesBeh = inds;
        % BEHAVIOR OF INDICES
        else
            spikesBeh = beh(inds,:);
        end
        violations = abs(beh.time(inds) - spikeTimes(:)) > Opt.matchTimeTolerance;
        switch Opt.violationHandling
            case 'remove'
            spikeTimes(violations) = nan;
            case 'nan'
            spikeTimes = spikeTimes(~violations);
        end
    end

% --------------------------
% Main : No spikes provided
% --------------------------
else
    if multipleShifts
        spikesBeh = {};
    else
        spikesBeh = table();
    end
end

% --------------
% Prepare output
% --------------
if strcmp(Opt.output,'cellmatrix')

    if Opt.concatShifts
        spikesBeh = cat(2, spikesBeh{:});
    end
    if ~Opt.returnIndices
        for i = 1:numel(spikesBeh)
            spikesBeh{i} = table2array(spikesBeh);
        end
    end
    spikesBeh = util.type.castefficient(spikesBeh);
    
elseif strcmp(Opt.output, 'table')

    if Opt.returnIndices 
        indicesRange = [1, height(beh)];
        spikesBeh = util.type.castefficient(spikesBeh,...
            'negativeNan', true, 'minmaxOverrid', indicesRange);
        if iscell(spikesBeh) && ~istable(spikesBeh{1})
            for i = 1:numel(spikesBeh)
                spikesBeh{i} = table(spikesBeh{i}, 'VariableNames', "indices");
            end
        elseif ~iscell(spikesBeh) && ~istable(spikesBeh)
            spikesBeh = table(spikesBeh, 'VariableNames', 'indices');
        else
            error("What the?")
        end
    end

    if Opt.concatShifts
        for i = 1:numel(spikesBeh)
            spikesBeh{i}.shift = i * ones(height(spikesBeh{i}), 1);
        end
        spikesBeh = vertcat(spikesBeh{:});
        spikesBeh.shift = util.type.castefficient(spikesBeh.shift);
    end

    if Opt.useGPU
        spikeTimes = gather(spikeTimes);
        spikesBeh  = util.table.GPUtable2table(spikesBeh);
    end

else
    error("unrecognized return type");
end
