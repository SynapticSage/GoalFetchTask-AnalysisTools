function out = jercog_allNeurons(spikestruct, behstruct, varargin)
%function out = jercog_allNeurons(firing, position, postimes, H, varargin)
% Reference point firing rate analysis from the columbia people: KANDEL/ABBOTT
% inspired by the Sarel paper.
%
% Inputs:
% -------
% Takes t (postimes), x, y (position) , and H at spike times for a cell
%
% Key questions:
% -------------
% - are reference points dynamic (change with homewell and barrier structure?)
% - does pfc fit with reference point?
% - references concentrate at goal locations?
% - first quarter epoch of the day (or first epoch of the day with static
% barrier) match the last quarter epoch of the day?



if ~iscell(firing)
    firing   = double(firing);
    assert(size(firing,1) > size(firing,2));
else
    inds = find(cellfun(@isempty, firing));
    for ind = inds
        firing{ind} = inf; % empty sets not allowed: inf ensures will not match any times
    end
end
position = double(position);
H        = double(H);
nNeurons = size(firing,2);

ip = inputParser;
ip.addParameter('acceptMethod', 'or'); % Conservative: and, Liberal, or
ip.addParameter('grid', 50);
ip.addParameter('gridH', 8);
ip.addParameter('regularizeBounds', false);
ip.addParameter('type', 'egocentric');
ip.addParameter('requiredDegreeSample', 50);
ip.addParameter('requiredSampleTime', 0.5); % Half a second
ip.addParameter('samprate', []); % Half a second
ip.addParameter('spikeTimeToBehaviorTolerance', []);
ip.addParameter('method', "fminsearch");
ip.parse(varargin{:})
Opt = ip.Results;
if isempty(Opt.spikeTimeToBehaviorTolerance)
    Opt.spikeTimeToBehaviorTolerance = median(diff(postimes));
end


if Opt.requiredSampleTime > 0
    assert(~isempty(Opt.samprate), 'Must provid samplerate if requiring a total duration of time for samples of x,y,H');
end

bounds   = [min(position); max(position)]';
edges    = arrayfun(@(x) linspace(bounds(x,1), bounds(x,2), Opt.grid+1), 1:size(bounds, 2), 'UniformOutput', false);
edges{3} = linspace(min(H), max(H), Opt.gridH);
centers  = arrayfun(@(x) mean([edges{x}(1:end-1)', edges{x}(2:end)'],2) ,1:3, 'UniformOutput', false);

% ----------------------------
% X-Y position binned measures
% ----------------------------
xy = struct('x',nan, 'y', nan, 'N', 0, 'avgFR', nan, 'sumFR', nan, 'FR_occNorm',  nan);
xy = repmat(xy, [nNeurons cellfun(@numel, centers(1:2))]);
X = position(:,1);
Y = position(:,2);
x = discretize(X, edges{1});
y = discretize(Y, edges{2});
[groups, x, y] = findgroups(x, y);
uGroups = unique(groups);
for group = progress(uGroups')
    if isnan(group)
        continue
    end
    inds = groups == group;
    
    if numel(sum(inds)) > 1
        keyboard
    end

    ix = x(group);
    iy = y(group);
    xcenter = centers{1}(ix);
    ycenter = centers{2}(iy);
    N = sum(inds);

    % Calculate aveerage firing for this bin and occupancy
    if iscell(firing) % Sparse form
        ptime = postimes(group == groups);
        parfor i = 1:numel(firing)
            times = firing{i};
            %spikesInTimes = ismembertol(times, ptime, Opt.spikeTimeToBehaviorTolerance);
            sumFR{i} = sum(...
                any(abs(times - ptime) < Opt.spikeTimeToBehaviorTolerance, 1)...
                );
            FR_occNorm{i} = sumFR{i}/N;
        end
        avgFR = num2cell(nan(1,numel(firing)));
    else % Dense form
        avgFR = nanmean(firing(inds,:),1);
        sumFR = nansum(firing(inds,:),1);
        FR_occNorm = sumFR/N;
        avgFR = num2cell(avgFR);
        sumFR = num2cell(sumFR);
        FR_occNorm = num2cell(FR_occNorm);
    end
    % Struct for storing f(x,y) data
    [xy(:, ix, iy).x]          = deal(xcenter);
    [xy(:, ix, iy).y]          = deal(ycenter);
    [xy(:, ix, iy).N]          = deal(sum(inds)); % timespent is proportional to this
    [xy(:, ix, iy).avgFR]      = deal(avgFR{:});
    [xy(:, ix, iy).sumFR]      = deal(sumFR{:});
    [xy(:, ix, iy).FR_occNorm] = deal(FR_occNorm{:});
end
    keyboard

% ----------------------------
% X-Y-H position binned measures
% ----------------------------
disp("Place-HD tuning")
xyh = struct('x',nan','y',nan','h',nan,'N',0,'avgFR',nan,'sumFR',nan,'FR_occNorm',nan, 'totaltime', 0, 'enoughtime', false);
xyh = repmat(xyh, [nNeurons, cellfun(@numel, centers)]);
x = discretize(X, edges{1});
y = discretize(Y, edges{2});
h = discretize(H, edges{3});
[groups, x, y, h] = findgroups(x, y, h);
uGroups = unique(groups);
no_inds = 0
for group = progress(uGroups')
    if isnan(group)
        continue
    end

    inds = groups == group;

    ix = x(group);
    iy = y(group);
    ih = h(group);
    xcenter = centers{1}(ix);
    ycenter = centers{2}(iy);
    hcenter = centers{3}(ih);

    if iscell(firing) % Sparse form
        ptime = postimes(group);
        N = sum(inds);
        for i = 1:numel(firing)
            times = firing{i};
            sumFR{i} = sum(...
                any(abs(times - ptime) < Opt.spikeTimeToBehaviorTolerance, 1)...
                );
            FR_occNorm{i} = sumFR/N;
        end
        avgFR = num2cell(nan(1,numel(firing)));
    else % Dense form
        % Calculate aveerage firing for this bin and occupancy
        avgFR = nanmean(firing(inds,:),1);
        sumFR = nansum(firing(inds,:), 1);
        FR_occNorm = sumFR/sum(inds);
        avgFR = num2cell(avgFR);
        sumFR = num2cell(sumFR);
        FR_occNorm = num2cell(FR_occNorm);
    end
    N = sum(inds);

    % Struct for storing f(x,y) data
    [xyh(:, ix, iy, ih).x]          = deal(xcenter);
    [xyh(:, ix, iy, ih).y]          = deal(ycenter);
    [xyh(:, ix, iy, ih).h]          = deal(hcenter);
    [xyh(:, ix, iy, ih).N]          = deal(N); % timespent is proportional to this
    [xyh(:, ix, iy, ih).avgFR]      = deal(avgFR{:});
    [xyh(:, ix, iy, ih).sumFR]      = deal(sumFR{:});
    [xyh(:, ix, iy, ih).FR_occNorm] = deal(FR_occNorm);
    [xyh(:, ix, iy, ih).totaltime]  = deal(N * 1/(Opt.samprate));
    [xyh(:, ix, iy, ih).enoughtime] = deal(N * 1/(Opt.samprate) > Opt.requiredSampleTime);
end

% Get place field centroids
% -------------------------
disp('Calculating place field centroids')
fr      = arrayfun(@(neuron) nd.fieldGet(xy(neuron, :,:), 'avgFR'), 1:nNeurons, 'UniformOutput', false);
fr      = cellfun(@squeeze, fr, 'UniformOutput', false);
centroid = [];
for i = 1:numel(fr)
    %[~, fr{i}]   = max(fr{i}, [], 'all', 'linear'); % find max average firing rate
    %[ii{i},jj{i}]    = ind2sub(size(xx{i}), I{i});
    %centroid{i} = [centers{1}(ii{i}), centers{2}(jj{i})];
    [Cx, Cy] = ndgrid(centers{1}, centers{2});
    xcoord = nanmean(fr{i} .* Cx, 'all')./nanmean(fr{i}, 'all');
    ycoord = nanmean(fr{i} .* Cy, 'all')./nanmean(fr{i}, 'all');
    centroid = [centroid; xcoord, ycoord];
end


% PRODUCE MODEL WITH UNIFIED REFERENCE POINT
% ------------------------------------------
disp('RH tuning: getting field data')
% Compute heading direction stats
h_prob  = nd.fieldGet(xyh, 'N');
h_prob  = bsxfun(@rdivide, h_prob, nansum(h_prob, 3));

N_xyh = nd.fieldGet(xyh, 'N');
N_xy  = nd.fieldGet(xy,  'N');
fr_xyh = nd.fieldGet(xyh, 'FR_occNorm');
fr_xy  = nd.fieldGet(xy,  'FR_occNorm');
keyboard;
h_set = nd.fieldGet(xyh, 'h');
x_set = nd.fieldGet(xyh, 'x');
y_set = nd.fieldGet(xyh, 'y');
h_width = 360 / Opt.gridH;
enoughtime = nd.fieldGet(xyh, 'enoughtime');
Hdim = 4;
degreeSampled = sum(enoughtime, Hdim) * h_width;
enoughdegree = degreeSampled > Opt.requiredDegreeSample;
acceptSample = bsxfun(@and, enoughdegree, enoughtime);
if Opt.acceptMethod == "and"
    acceptSample = prod(acceptSample, Hdim);
elseif Opt.acceptMethod == "or"
    acceptSample = sum(acceptSample, Hdim) > 0;
end
acceptSample = repmat(acceptSample, 1, 1, size(N_xyh, Hdim));

% Modulation ratio
disp('RH tuning: modulation ratio')
rxyh_over_rxy = bsxfun(@rdivide, fr_xyh, fr_xy);
rxyh_over_rxy(~acceptSample) = nan;

% Perform fminsearch, as in the paper, to find the reference points
disp('Fitting')
searchFunc        = @(neuron, X) coding.vectorCells.multiNeuronJercogModel(rxyh_over_rxy(neuron,:,:,:,:), x_set(neuron,:,:,:,:), y_set(neuron,:,:,:,:), h_set(neuron,:,:,:,:), h_prob(neuron,:,:,:,:), X(1), X(2), X(3:4));
%initialConditions = [repmat(0,nNeurons,1), repmat(0,nNeurons,1), centroid(:,1), centroid(:,2)]; % g=0, theta_p=0, ref=centroid of the place field
if Opt.regularizeBounds
    searchFunc = @(X) regularizeBounds(searchFunc(X), X(3:4), bounds);
end

% Ready fit process
disp('Fitting jercog et al. model via patternSearch mesh method');
if Opt.method == "fminsearch"
    options = optimset();
    options.MaxFunEvals = 1e6;
    options.OptimalityTolerance = 1e-9;
    options.CheckGradients = true;
    options.DerivativeCheck = true;
    options.PlotFcns = @optimplotfval;
    %options.TolX = 1e-8;
elseif Opt.method == "patternsearch"
    options = optimoptions('patternsearch');
    options.MeshExpansionFactor   = 1.25;
    options.MeshContractionFactor = 0.75;
    options.UseCompletePoll   = true;
    options.UseCompleteSearch = true;
    A = eye(4); b = [100, pi, 1000, 1000];
    Am = -A; bm = -b;
    %A = [A; Am]; 
    options.MeshTolerance = 1e-9;
end

%options.FunctionTolerance = 1e-8;
%options.MaxIter = 1e7;
%options.MaxFunctionEvaluations = 1e9;
%options.StepTolerance = 1e-8;
initialConditions = [repmat(0,nNeurons,1), repmat(0,nNeurons,1), centroid(:,1), centroid(:,2)]; % g=0, theta_p=0, ref=centroid of the place field
%initialConditions(:,3:4) = initialConditions(:,3:4) * 1000;

params = nan(nNeurons, 4);
for n = progress(1:nNeurons, 'Title', char(Opt.method))
    if any(isnan(centroid(n,:)))
        continue
    end
    if Opt.method == "patternsearch"
        params(n,:) = patternsearch(@(X)  searchFunc(n, X), initialConditions(n, :), [], [],[],[],bm,b,[],  options); 
    else
        [params(n,:), fval(n), exitflag(n), output(n)] = fminsearch(@(X)  searchFunc(n, X), initialConditions(n, :), options); 
    end
end

% Organize outputs
amplitude        = params(:,1);
theta_preference = mod(params(:,2)+pi,2*pi)-pi;
reference_point  = params(:,3:4);
params           = [amplitude, theta_preference, reference_point];


out = struct('xy', xy, 'xyh', xyh, 'rxyh_over_rxy', rxyh_over_rxy, 'centroid', centroid);
out.amplitude        = amplitude;
out.theta_preference = theta_preference;
out.reference_point  = reference_point;

function output = model(rxyh_over_rxy, x, y, h_vals, h_prob, g, theta_preference, ref)
% The actualy nitty gritty implementation of the model

    theta_preference  = mod(theta_preference + pi, 2*pi) - pi; % 1x1

    % LAMBDAS
    % -------
    ref_angle_func = @(h, ref) atan(y-ref(2))./(x-ref(1)); % Lambda : H(x,y,h) x 1                             % Angle of animal position relative to reference point
    theta_func = @(ref_angle, h_vals) ref_angle - h_vals;  % Lambda : H(x,y,h) x 1                             % Heading relative to reference
    F_func     = @(ref_angle, h_vals, theta_preference) cos(theta_func(ref_angle, h_vals) - theta_preference); % Circular modulation by a prefered heading angle to reference
                                                                                                               % Lambda : Hx1

    ref_angle = ref_angle_func(h_vals, ref);                 % H(x,y,h) x 1, ANGLE TO REFERENCE
    F         = F_func(ref_angle, h_vals, theta_preference); % H(x,y,h) x 1, ranges 0 to 1, FEATURE MATCH


    % Should be function of x any y
    F_avg  = nanmean(h_prob .* F, 3); % expectation of F over heading directions, AVERAGE FEATURE MATCH

    % model output and imposition of model fit in jercog, Rxyh_model > 0 => abs(output)
    output = abs(1 + g .* bsxfun(@minus, F, F_avg)); % abs() because cannot have negative fr, 1 + modulated * (feature match relative to average)
    output = nanmean((output - rxyh_over_rxy).^2, 'all') + 10*inverseHeavisde(g); % mean sqaured eerror, penalizing g < 0

function output = inverseHeavisde(x)
    if x < 0
        output = 1;
    else
        output = 0;
    end


function output = regularizeBounds(x, ref, bounds)
    % Regularization function to penalize for references outside of the the camera's
    % recorded bounds.
    
    XCOORD = 1;
    YCOORD = 2;
    LOW    = 1;
    HIGH   = 2;
    
    if ref(HIGH) > bounds(YCOORD,HIGH) || ref(HIGH) < bounds(YCOORD,LOW) ...
        || ref(LOW) > bounds(XCOORD,HIGH) || ref(LOW) < bounds(XCOORD,LOW)
        output = output + 1;
    end


function x = expandNeuronDimension(x)
    % expands neuron dimension of structs with vector field entries


