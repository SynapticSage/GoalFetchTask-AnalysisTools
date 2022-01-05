function out = jercog_allNeurons_lean(spikes, beh, varargin)
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
ip.addParameter('useGPU', true);                     % Should we try to use gpu (esp for fminsearch)
ip.addParameter('useGPU_calcfield', true);           % Should we try to use gpu for calculating xyh fields?
ip.addParameter('acceptMethod', 'and');              % Conservative: and, Liberal, or
ip.addParameter('grid', 28);                         % Size of the linear binning in xy
ip.addParameter('gridH', 8);                         % Sizee of the binning in h
ip.addParameter('regularizeBounds', false);          % Do we regularize bounds?
ip.addParameter('type', 'egocentric');               % What type of H to use?
ip.addParameter('requiredDegreeSample', []);         % Default = (360 / gridH) + 1, which gaurentees at least 2 bins
ip.addParameter('requiredSampleTime', 1/30 * 5);     % 7.5 frames worth of sample time
ip.addParameter('samprate', []);                     % Half a second
ip.addParameter('spikeTimeToBehaviorTolerance', []); % 
ip.addParameter('method', "fminsearch");             % 
ip.addParameter('location', "jercog");
ip.parse(varargin{:})
Opt = ip.Results;
if isempty(Opt.spikeTimeToBehaviorTolerance)
    Opt.spikeTimeToBehaviorTolerance = median(diff(beh.time));
end
if isempty(Opt.requiredDegreeSample)
    Opt.requiredDegreeSample = 360/Opt.gridH + 1;  % default to minimimally two bins of angle
end
if ismember('pos', fieldnames(beh))
    beh.x   = beh.pos(:,1);
    beh.y   = beh.pos(:,2);
end
beh.h   = angle(beh.velVec);
spikes.beh.h   = double(angle(spikes.beh.velVec));
beh.pos = double(beh.pos);
beh.h   = double(beh.h);

try

% -----------------
% Figure out bounds
% -----------------
bounds   = [min(beh.pos); max(beh.pos)]';
edges    = arrayfun(@(x) linspace(bounds(x,1), bounds(x,2), Opt.grid+1), 1:size(bounds, 2), 'UniformOutput', false);
edges{3} = linspace(min(beh.h), max(beh.h), Opt.gridH);
centers  = arrayfun(@(x) mean([edges{x}(1:end-1)', edges{x}(2:end)'],2) ,1:3, 'UniformOutput', false);
uNeurons = unique(spikes.beh.neuron);
nNeurons = max(uNeurons);
samplingPeriod = median(diff(beh.time));
nCenters = cellfun(@numel, centers);
nCentersCell = num2cell(nCenters);

% -------------------------
% Get binned FR_occNorm
% -------------------------

Grid = struct('x', Opt.grid, 'y', Opt.grid, 'h', Opt.gridH);
[spikes_xy, beh_xy, scaff_xy] = coding.field.calc(spikes.beh,...
    'props', ["x", "y"],...
    'beh', beh,...
    'useGPU', Opt.useGPU && Opt.useGPU_calcfield,...
    'grid', Grid);
[spikes_xyh, beh_xyh, scaff_xyh] = coding.field.calc(spikes.beh,...
    'props', ["x", "y", "h"],...
    'beh', beh,...
    'useGPU', Opt.useGPU && Opt.useGPU_calcfield,...
    'grid', Grid);

% -------------------------
% Get place field centroids
% -------------------------
disp('Calculating place field centroids')
fr      = spikes_xy.FR_occNorm;
centers = nd.fieldGetCell(scaff_xyh, 'center');
centroid = nan(size(fr,1),2);
for ii = 1:size(fr,1)
    ith = @(ii) squeeze(fr(ii,:,:,:));
    [M, xx] = max(ith(ii), [], 1);
    [~, yy] = max(M, [], 2);
    xx = xx(yy);
    centroid(ii,:) = [centers{1}(xx), centers{2}(yy)];
end

% ------------------------------------------
% PRODUCE MODEL WITH UNIFIED REFERENCE POINT
% ------------------------------------------
disp('RH tuning: getting field data')
% Compute heading direction stats
h_prob  = spikes_xyh.FR_occNorm;
h_axis  = ndims(h_prob);
h_prob  = bsxfun(@rdivide, h_prob, nansum(h_prob, h_axis));

N_xyh  = spikes_xyh.spikeTime;
N_xy   = spikes_xy.spikeTime;
fr_xyh = spikes_xyh.FR_occNorm;
fr_xy  = spikes_xy.FR_occNorm;

[x_set, y_set, h_set] = meshgrid(centers{:});

% ------------
% Constraints!
% ------------
% 1) Enough time spent in each bin?
% 2) Did we sample enough of the full spectrum of degrees in each (x,y)?
% 3) Do we jointly have both of these? -- if so, accept the sample
h_width = 360 / Opt.gridH;
enoughtime = N_xyh > Opt.requiredSampleTime;
Hdim = 4;
degreeSampled = sum(enoughtime, Hdim) * h_width;
enoughdegree  = degreeSampled > Opt.requiredDegreeSample;
acceptSample  = bsxfun(@and, enoughdegree, enoughtime);
if Opt.acceptMethod == "and"
    acceptSample = prod(acceptSample, Hdim);
elseif Opt.acceptMethod == "or"
    warning('Liberal OR-based accept')
    acceptSample = sum(acceptSample, Hdim) > 0;
end
acceptSample = repmat(acceptSample, 1, 1, size(N_xyh, Hdim));
out.fracEnoughTime       = mean(enoughtime,'all');
out.acceptedSamples      = mean(acceptSample,'all');
out.fracAcceptableDegree = mean(enoughdegree,'all');
disp("Fraction with acceptable time → " + sprintf('%2.2f',   out.fracEnoughTime));
disp("Fraction with acceptable degree → " + sprintf('%2.2f', out.fracAcceptableDegree));
disp("Fraction of accepted samples → " + sprintf('%2.2f',    out.acceptedSamples));
if mean(acceptSample, 'all') == 0
    error("You have no accepted samples");
end

% Compute modulation ratio
disp('RH tuning: modulation ratio')
rxyh_over_rxy = bsxfun(@rdivide, fr_xyh, fr_xy);
rxyh_over_rxy(~acceptSample) = nan;

% Ready fit process
if Opt.method == "fminsearch"
    disp('Fitting jercog et al. model via fminsearch mesh method');
    options = optimset();
    options.MaxFunEvals         = 2e6;
    options.OptimalityTolerance = 1e-9;
    options.CheckGradients      = true;
    options.DerivativeCheck     = true;
    %options.PlotFcns = @optimplotfval;
    %options.TolX = 1e-8;
elseif Opt.method == "patternsearch"
    disp('Fitting jercog et al. model via patternSearch mesh method');
    options = optimoptions('patternsearch');
    options.MeshExpansionFactor   = 1.25;
    options.MeshContractionFactor = 0.75;
    options.UseCompletePoll       = true;
    options.UseCompleteSearch     = true;
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

catch ME
    if contains(ME.message, 'gpu', 'IgnoreCase', true)
        disp("Gpu  error. Rerunning without GPU");
        Opt.useGPU = false;
        spikes = coding.jercog.allNeurons_sparse_lean(spikes, beh, Opt);
    else
        throw(ME);
    end
end

if Opt.useGPU
    out = util.struct.GPUstruct2struct(out);
end

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



