function T = rateBy(spikes, by, varargin)
% Creates tensorized rate matrix, where the tensor dimensions encode different
% categoricals or time variables
%
% used by many different analysees
%
% Major work horse of getting data ready to input into dpca and dsca

ip = inputParser;
%% Spikes
%% COLUMN-CONTROL
ip.addParameter('col_timeChanging', []); %Lowest level can be pure time or pure distance
ip.addParameter('col_categoryVars', []);
ip.addParameter('col_quantileVars', []);
%% Norm-lowest-timeChagning var?
ip.addParameter('norm_timeChangingLowest', true); % If your trials are not equal time/distance, this has to be true
%% MODES
ip.addParameter('dynamicReduceVar', []); % if off, just standard
%% OPTIONS
ip.addParameter('warpToFit', 'nope'); % {nope} | dynamic | static
% where dynamic = dtw and static is simple interpolation squishing to the same size
ip.addParameter('warpSize', 20); % number of time slots if we're warping
ip.parse(varargin{:})
Opt = ip.Results;

if isempty(Opt.col_timeChangingVars)
    error('Must have timeChangingVars')
else
    tVars   = Opt.col_timeChangingVars;
    tGroups = util.group.findgroups(by, tVars);
end

% Find categorical grroups
cVars     = Opt.col_categoryVars;
cGroups = util.group.findgroups(by, cVars);

% Find quantile grroups
qVars   = Opt.col_quantileVars;
for qVar = qVars
    by.(qVar) = util.num.getQuantile(by.(qVar));
end
qGroups = util.group.findgroups(by, qVars);

% Find temporal grroups
tVars   = setdiff(Opt.col_timeChangingVars, 'time');
for tVar = tVars
    by.(tVar) = util.num.getQuantile(by.(tVar));
end
tGroups = util.group.findgroups(by, tVars);
if ~isempty(Opt.warpSize)
    smallest_time_dim = Opt.warpSize;
else
    smallest_time_dim = gbehavior.distributionOfMatchingTimes(beh, {cGroups, qGroups, tGroups});
end


nNeurons  = size(rate, 1);
cDims = cellfun(@numel, cGroups.group.values);
qDims = cellfun(@numel, qGroups.group.values);
tDims = cellfun(@numel, tGroups.group.values);
tDimsFinal = [tDims, smallest_time_dim];
if ~isempty(Opt.dynamicReduceVar)
    tDimsFinal(Opt.dynamicReduceVar) = 1;
end

% Intitialize the data matrix
dimsFinal = [cDims, qDims, tDimsFinal];
X = nan(dimsFinal, 'single');


for c = cGroups.uGroups'
for q = qGroups.uGroups'
for t = tGroups.uGroups'

    % Matrix address
    C = cGroups.group.addressByGroupNum{c};
    Q = cGroups.group.addressByGroupNum{q};
    T = cGroups.group.addressByGroupNum{t};

    % Find data corresponding to our cuts
    dataCut = cGroups.time.groups == c & ...
              qGroups.time.groups == q & ...
              tGroups.time.groups == t;
    
    matchingTimes = by.time(dataCut);


    % if the type of the spike data is sparse,
    % then we should probably translate it to
    % a dense format for the following code.
    if strcmpi(spikes.type,'sparse') && Opt.sparseOnTheFly
        minT = min(matchingTimes);
        maxT = min(matchingTimes);
        x = units.sparseToDense(spikes,...
            'returnOutput', 'matrix',...
            'startTime', minT, ...
            'endTime', maxT ...
            'method', 'instanteous'
            );
        % TODO if this doesn't return regularly latticed time samples, then
        % i should pass in an interpolation arugment.
        
    end


    % Take the little block that we recovered, and elaborate it onto
    % our tensor
    X(C, Q, T, :) = x;


    if Opt.dynamicReduceVar
        higherBlockChanged = ~all(T(0:end-1) ~= T(1:end-1));
        % Reduction
        Tprev = T;
    end

end
end
end





