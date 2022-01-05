function [spikes] = jercog_allNeurons(spikes, beh, varargin)
% function out = jercog_allNeurons(firing, position, postimes, H, varargin)
% Reference point firing rate analysis from the columbia people: KANDEL/ABBOTT
% inspired by the Sarel paper.
%
% -------
% Inputs:
% -------
% Takes t (postimes), x, y (position) , and H at spike times for a cell
%
% --------------
% Key questions:
% --------------
% - are reference points dynamic (change with homewell and barrier structure?)
% - does pfc fit with reference point?
% - references concentrate at goal locations?
% - first quarter epoch of the day (or first epoch of the day with static
% barrier) match the last quarter epoch of the day?


if iscell(spikes.beh)
        spikes.beh = nd.dimLabel(spikes.beh(:), 1, 'neuron');
        spikes.beh = cat(1, spikes.beh{:});
end

ip = inputParser;
ip.addParameter('useGPU', true); % Conservative: and, Liberal, or
ip.addParameter('acceptMethod', 'or'); % Conservative: and, Liberal, or
ip.addParameter('grid', 28);
ip.addParameter('gridH', 8);
ip.addParameter('regularizeBounds', false);
ip.addParameter('type', 'egocentric');
ip.addParameter('requiredDegreeSample', []);
ip.addParameter('requiredSampleTime', 0.5); % Half a second
ip.addParameter('samprate', []); % Half a second
ip.addParameter('spikeTimeToBehaviorTolerance', []);
ip.addParameter('method', "fminsearch");
ip.addParameter('location', "jercog");
ip.parse(varargin{:})
Opt = ip.Results;
if isempty(Opt.spikeTimeToBehaviorTolerance)
    Opt.spikeTimeToBehaviorTolerance = median(diff(beh.time));
end

beh.X   = beh.pos(:,1);
beh.Y   = beh.pos(:,2);
beh.H   = angle(beh.velVec);
spikes.beh.H   = double(angle(spikes.beh.velVec));
beh.pos = double(beh.pos);
beh.H   = double(beh.H);

% -----------------
% Figure out bounds
% -----------------
bounds   = [min(beh.pos); max(beh.pos)]';
edges    = arrayfun(@(x) linspace(bounds(x,1), bounds(x,2), Opt.grid+1), 1:size(bounds, 2), 'UniformOutput', false);
edges{3} = linspace(min(beh.H), max(beh.H), Opt.gridH);
centers  = arrayfun(@(x) mean([edges{x}(1:end-1)', edges{x}(2:end)'],2) ,1:3, 'UniformOutput', false);
uNeurons = unique(spikes.beh.neuron);
nNeurons = max(uNeurons);
samplingPeriod = median(diff(beh.time));
nCenters = cellfun(@numel, centers);
nCentersCell = num2cell(nCenters);

% ---- 
% Bin XY
% ----
beh.xy_binx = discretize(beh.X, edges{1});
beh.xy_biny = discretize(beh.Y, edges{2});
[beh.xy_groups, ~, ~] = findgroups(beh.xy_binx, beh.xy_biny);
% ------------------
% Summarize behavior XY
% ------------------
beh_xy = struct();
euc = (beh.pos(:,1) - centers{1}').^2;
[~, beh.x_bin] = min(euc, [], 2);
euc = (beh.pos(:,2) - centers{2}').^2;
[~, beh.y_bin] = min(euc, [], 2);
beh.xy_bin = (beh.y_bin-1)*numel(centers{2}) + beh.x_bin;

beh_xy.x = centers{1};
beh_xy.y = centers{2};
[N, E] = histcounts(beh.xy_bin, unique(beh.xy_bin));
beh_xy.visits = zeros(numel(centers{2}), numel(centers{1}));
beh_xy.visits(E(1:end-1)) = N;
beh_xy.visit_time = beh_xy.visits * samplingPeriod;

% ----------------------------
% X-Y position neural measures
% ----------------------------
disp("Comuting XY measures");
tic
euc = (spikes.beh.pos(:,1) - centers{1}').^2;
[~,spikes.beh.x_bin] = min(euc, [], 2);
euc = (spikes.beh.pos(:,2) - centers{2}').^2;
[~,spikes.beh.y_bin] = min(euc, [], 2);
clear euc 
spikes.beh.xy_bin = (spikes.beh.y_bin-1)*numel(centers{2}) + spikes.beh.x_bin;

%spikes_xy = struct('x',nan, 'y', nan, 'N', 0, 'avgFR', nan, 'sumFR', nan, 'FR_occNorm',  nan);
spikes_xy = struct(...
    'spikeCount',zeros(max(uNeurons), nCentersCell{1:2}),...
    'spikeTime',zeros(max(uNeurons), nCentersCell{1:2})...
    );
neurons = uNeurons(:)';
for neuron = progress(neurons)
    nfilt = spikes.beh.neuron == neuron;
    subset = spikes.beh(nfilt,:);
    [N,E] = histcounts(subset.xy_bin, unique(beh.xy_bin));
    spikes_xy.spikeCount(neuron, E(1:end-1)) = N;
    
    % Bin specific calculations
    for bin = unique(beh.xy_bin)'
        filt = subset.xy_bin == bin;
        uTimes = unique(spikes.beh.time(filt));
        spikes_xy.spikeTime(neuron, bin) = sum(numel(uTimes)) * samplingPeriod; % Sum of sampling periods that have spikes
    end
end
spikes_xy.FR_occNorm = bsxfun(@rdivide, spikes_xy.spikeCount, shiftdim(beh_xy.visit_time,-1));
toc


disp("Compputing XYH measures")

beh.xyh_binx = discretize(beh.X, edges{1});
beh.xyh_biny = discretize(beh.Y, edges{2});
beh.xyh_binh = discretize(beh.H, edges{3});
[beh.xyh_groups, ~, ~] = findgroups(beh.xyh_binx, beh.xyh_biny, beh.xyh_binh);

euc = (beh.pos(:,1) - centers{1}').^2;
[~, beh.x_bin] = min(euc, [], 2);
euc = (beh.pos(:,2) - centers{2}').^2;
[~, beh.y_bin] = min(euc, [], 2);
euc = (beh.H - centers{3}').^2;
[~, beh.h_bin] = min(euc, [], 2);
beh.xyh_bin =  (beh.h_bin-1)*prod(nCenters(1:2)) + (beh.y_bin-1)*nCenters(1) + beh.x_bin;

beh_xyh = struct();
beh_xyh.x = centers{1};
beh_xyh.y = centers{2};
beh_xyh.h = centers{3};
[N, E] = histcounts(beh.xyh_bin, unique(beh.xyh_bin));
beh_xyh.visits = zeros(nCentersCell{:});
beh_xyh.visits(E(1:end-1)) = N;
beh_xyh.visit_time = beh_xyh.visits * samplingPeriod;

spikes_xyh = struct(...
    'spikeCount',zeros(max(uNeurons), nCentersCell{:}),...
    'spikeTime',zeros(max(uNeurons), nCentersCell{:})...
    );

euc = (spikes.beh.pos(:,1) - centers{1}').^2;
[~,spikes.beh.x_bin] = min(euc, [], 2);
euc = (spikes.beh.pos(:,2) - centers{2}').^2;
[~,spikes.beh.y_bin] = min(euc, [], 2);
euc = (spikes.beh.H - centers{3}').^2;
[~,spikes.beh.h_bin] = min(euc, [], 2);
clear euc
spikes.beh.xyh_bin =  (spikes.beh.h_bin-1)*prod(nCenters(1:2)) + (spikes.beh.y_bin-1)*nCenters(1) + spikes.beh.x_bin;
for neuron = progress(neurons)
    nfilt = spikes.beh.neuron == neuron;
    subset = spikes.beh(nfilt,:);
    [N,E] = histcounts(subset.xyh_bin, unique(beh.xyh_bin));
    spikes_xyh.spikeCount(neuron, E(1:end-1)) = N;
    
    % Bin specific calculations
    for bin = unique(beh.xyh_bin)'
        filt = subset.xyh_bin == bin;
        uTimes = unique(spikes.beh.time(filt));
        spikes_xyh.spikeTime(neuron, bin) = sum(numel(uTimes)) * samplingPeriod; % Sum of sampling periods that have spikes
    end
end
spikes_xyh.FR_occNorm = bsxfun(@rdivide, spikes_xyh.spikeCount, shiftdim(beh_xyh.visit_time,-1));
toc

% ----------------------------
% X-Y-H position binned measures
% ----------------------------
disp("Place-HD tuning")

% Get place field centroids
% -------------------------
disp('Calculating place field centroids')
clf
fr      = spikes_xy.FR_occNorm;
centroid = nan(size(fr,1),2);
for ii = 1:size(fr,1)
    %[~, fr{ii}]   = max(fr{ii}, [], 'all', 'linear'); % find max average firing rate
    %[ii{ii},jj{ii}]    = ind2sub(size(xx{ii}), I{ii});
    %centroid{ii} = [centers{1}(ii{ii}), centers{2}(jj{ii})];
    %[Cx, Cy] = ndgrid(centers{1}, centers{2});
    ith = @(ii) squeeze(fr(ii,:,:,:));
    %xcoord = nanmean(ith(ii) .* Cx, 'all')./nanmean(ith(ii), 'all');
    %ycoord = nanmean(ith(ii) .* Cy, 'all')./nanmean(ith(ii), 'all');
    [M, xx] = max(ith(ii), [], 1);
    [~, yy] = max(M, [], 2);
    xx = xx(yy);
    centroid(ii,:) = [centers{1}(xx), centers{2}(yy)];
    nexttile;
    imagesc(centers{2}, centers{1}, (squeeze(fr(ii,:,:))));
    hold on
    s =scatter(centroid(ii,2), centroid(ii,1), 20, 'white','filled');
    %if ii == 30
    %    keyboard
    %end
end


% PRODUCE MODEL WITH UNIFIED REFERENCE POINT
% ------------------------------------------
disp('RH tuning: getting field data')
% Compute heading direction stats
h_prob  = spikes_xyh.FR_occNorm;
h_axis  = ndims(h_prob);
h_prob  = bsxfun(@rdivide, h_prob, nansum(h_prob, h_axis)); % probability of firing at that angle given x,y,neuron

N_xyh = spikes_xyh.spikeTime;
N_xy  = spikes_xy.spikeTime;
fr_xyh = spikes_xyh.FR_occNorm;
fr_xy  = spikes_xy.FR_occNorm;

[x_set, y_set, h_set] = meshgrid(centers{:});

h_width = 360 / Opt.gridH;
enoughtime = N_xyh > Opt.requiredSampleTime;
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
disp("Fraction with acceptable time → " + sprintf('%2.2f', mean(enoughtime,'all')));
disp("Fraction with acceptable degree → " + sprintf('%2.2f', mean(enoughdegree,'all')));
disp("Fraction of accepted samples → " + sprintf('%2.2f', mean(acceptSample,'all')));
if mean(acceptSample, 'all') == 0
    error("You have no accepted samples");
end

% Compute modulation ratio
disp('RH tuning: modulation ratio')
rxyh_over_rxy = bsxfun(@rdivide, fr_xyh, fr_xy);
rxyh_over_rxy(~acceptSample) = nan;

% Ready fit process
if Opt.method == "fminsearch"
    options = optimset();
    options.MaxFunEvals = 1e6;
    options.OptimalityTolerance = 1e-9;
    options.CheckGradients = true;
    options.DerivativeCheck = true;
    %options.PlotFcns = @optimplotfval;
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

initialConditions = [zeros(max(uNeurons),1), zeros(max(uNeurons),1),...
    centroid(:,1), centroid(:,2)]; % g=0, theta_p=0, ref=centroid of the place field
if Opt.useGPU
    fitfun = @coding.jercog.gpumodel;
else
    fitfun = @coding.jercog.model;
end

searchFunc        = @(neuron, X) fitfun(...
    rxyh_over_rxy(neuron,:,:,:,:), ...
    shiftdim(x_set(:,:,:),-1), ...
    shiftdim(y_set(:,:,:),-1), ...
    shiftdim(h_set(:,:,:),-1), ...
    h_prob(neuron,:,:,:), ... % PROBABILITY of firing at H PER NEURON, X, Y 
    X(1), X(2), X(3:4));

if Opt.regularizeBounds
    searchFunc = @(X) regularizeBounds(searchFunc(X), X(3:4), bounds);
end


pc = parcluster('local');
progBar = ProgressBar(nNeurons, ...
    'IsParallel', true, ...
    'Title', char(Opt.method) ...
    );


params = nan(nNeurons, 4);

%% FIGURE OUT PARAMS and EVALUATE
for n = progress(1:nNeurons, 'Title', char(Opt.method))
    if any(isnan(centroid(n,:))) || ~any(isnan(params(n,:)))
        continue
    end
    if Opt.method == "patternsearch"
        %params(n,:) = patternsearch(@(X)  searchFunc(n, X), initialConditions(n, :), [], [],[],[],bm,b,[],  options); 
    else
        [params(n,:), fval(n), exitflag(n), output(n)] = fminsearch(@(X)  searchFunc(n, X), initialConditions(n, :), options); 
        %params(1:n,:)
    end

    [goodness_of_fit(n), tmp] = searchFunc(n, params(n,:));
    evaluation(n,:,:,:) = tmp;
end

% Organize outputs
amplitude        = params(:,1);
theta_preference = mod(params(:,2)+pi,2*pi)-pi;
reference_point  = params(:,3:4);

out                  = struct();

out.x                = centers{1};
out.y                = centers{2};
out.h                = transpose(exp(1i*centers{3}));

out.xy               = spikes_xy;
out.xyh              = spikes_xyh;

out.h_prob           = h_prob;
out.h_score          = nanmean(bsxfun(@times,...
                               shiftdim(out.h, -(h_axis-1)), out.h_prob), h_axis);

out.rxyh_over_rxy    = rxyh_over_rxy;
out.goodness_of_fit  = goodness_of_fit;
out.model_output     = evaluation;

out.centroid         = centroid;
out.amplitude        = amplitude;
out.theta_preference = theta_preference;
out.reference_point  = reference_point;

eval("spikes." + Opt.location + " = out");

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

