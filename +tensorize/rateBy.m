function [X, T, data] = rateBy(spikes, by, varargin)
% Creates tensorized rate matrix, where the tensor dimensions encode different
% categoricals or time variables
%
% used by many different analysees
%
% Major work horse of getting data ready to input into dpca and dsca

ip = inputParser;
%% Spikes
ip.addParameter('sparseOnTheFly', false); % ram friendly
%% COLUMN-CONTROL
ip.addParameter('vars_timeChanging', []); %Lowest level can be pure time or pure distance
ip.addParameter('vars_category', []);
ip.addParameter('vars_quantile', []);
ip.addParameter('var_time', 'time'); % only one. the core time variable (could be space instead of time)
%% Norm-lowest-timeChagning var?
ip.addParameter('norm_timeChangingLowest', true); % If your trials are not equal time/distance, this has to be true
%% MODES
ip.addParameter('dynamicReduceVar', []); % if off, just standard
%% OPTIONS
ip.addParameter('warpToFit', 'nope'); % {nope} | dynamic | static
% where dynamic = dtw and static is simple interpolation squishing to the same size
ip.addParameter('warpSize', 20); % number of time slots if we're warping
ip.addParameter('quantileSize', 3); % discretization of quantile variable
%
ip.addParameter('verbose', false);
ip.parse(varargin{:})
Opt = ip.Results;

Opt.vars_timeChanging = string(Opt.vars_timeChanging);
Opt.vars_quantile = string(Opt.vars_quantile);
Opt.vars_category = string(Opt.vars_category);
data = struct('skip', 0);

if isempty(Opt.vars_timeChanging)
    error('Must have timeChangingVars')
end

% Find categorical grroups
cVars     = Opt.vars_category;
cGroups = util.table.findgroups(by, cVars);

% Find quantile grroups
qVars   = Opt.vars_quantile;
if ~isempty(qVars)
    for qVar = qVars
        by.(qVar) = util.num.getQuantile(by.(qVar));
        by.(qVar) = discretize(by.(qVar), Opt.quantileSize);
    end
    qGroups = util.table.findgroups(by, qVars);
else
    by.unity = ones(height(by), 1);
    qGroups = util.table.findgroups(by, 'unity');
end

% Find temporal grroups
tVars   = setdiff(Opt.vars_timeChanging, Opt.var_time);
tGroups = util.table.findgroups(by, tVars);
if ~isempty(Opt.warpSize)
    smallest_time_dim = Opt.warpSize;
else
    smallest_time_dim = gbehavior.distributionOfMatchingTimes(beh, {cGroups, qGroups, tGroups});
end


nNeurons  = size(spikes.data, 2);
cDims = cellfun(@numel, cGroups.group.values);
qDims = cellfun(@numel, qGroups.group.values);
tDims = cellfun(@numel, tGroups.group.values);
tDimsFinal = [tDims, smallest_time_dim];
if ~isempty(Opt.dynamicReduceVar)
    tDimsFinal(Opt.dynamicReduceVar) = 1;
end

% Intitialize the data matrix
dimsFinal = [cDims, qDims, tDimsFinal, nNeurons];
X = nan(dimsFinal, 'single');
Xt = nan([dimsFinal(1:end-1) 1], 'single');
data.vars = [Opt.vars_category, Opt.vars_quantile, Opt.vars_timeChanging, "time", "neuron"];

for c = progress(cGroups.uGroups', 'Title', 'Categorical groups')
for q = progress(qGroups.uGroups', 'Title', 'Quantile groups')
for t = tGroups.uGroups'

    % Matrix address
    C = cGroups.group.addressByGroupNum{c};
    Q = qGroups.group.addressByGroupNum{q};
    T = tGroups.group.addressByGroupNum{t};

    % Find data corresponding to our cuts
    dataCut = cGroups.time.groups == c & ...
              qGroups.time.groups == q & ...
              tGroups.time.groups == t;
    
    matchingTimes = by.time(dataCut);
    if isempty(matchingTimes)
        continue
    end

    minT = min(matchingTimes);
    maxT = max(matchingTimes);
    if Opt.verbose
        disp(maxT - minT);
    end


    % if the type of the spike data is sparse,
    % then we should probably translate it to
    % a dense format for the following code.
    if strcmp(spikes.type,'sparse') && Opt.sparseOnTheFly
        spikes = units.sparseToDense(spikes,...
            'returnOutput', 'matrix',...
            'startTime', minT, ...
            'endTime', maxT, ...
            'nSamples', Opt.warpSize,...
            'method', 'instanteous');
        % TODO if this doesn't return regularly latticed time samples, then
        % i should pass in an interpolation arugment.
        x = spikes.data;
        time = spikes.time;
        if numel(t) < Opt.warpSize
            data.skip = data.skip + 1;
            continue
        end
    elseif strcmp(spikes.type,'dense')
        selection = util.constrain.minmax(spikes.time, [minT, maxT]);
        if ~sum(selection)
            data.skip = data.skip + 1;
            continue
        end
        x = spikes.data(selection,:);
        time = spikes.time(selection);
        tnew = linspace(minT, maxT, Opt.warpSize);
        x = util.interp.interp1(time, x, tnew);
        time = tnew;
    else
        error("Invalid")
    end



    % Take the little block that we recovered, and elaborate it onto
    % our tensor
    Xt(C, Q, T, :)   = time;
    X(C, Q, T, :, :) = x;

    if ~isempty(Opt.dynamicReduceVar)
        higherBlockChanged = ~all(T(0:end-1) ~= T(1:end-1));
        % Reduction

        Tprev = T;
    end

end
end
end

