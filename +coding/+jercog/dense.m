function [xy, xyh, params, rxyh_over_rxy] = jercog(firing, position, H, varargin)
% Reference point firing rate analysis from the columbia people: KANDEL/ABBOTT
% inspired by the Sarel paper.
%
% Inputs:
% -------
% Takes x, y , and H at spike times for a cell
%
% Key questions:
% -------------
% - are reference points dynamic (change with homewell and barrier structure?)
% - does pfc fit with reference point?
% - references concentrate at goal locations?
% - first quarter epoch of the day (or first epoch of the day with static
% barrier) match the last quarter epoch of the day?

firing   = double(firing);
position = double(position);
H        = double(H);

ip = inputParser;
ip.addParameter('acceptMethod', 'or'); % Conservative: and, Liberal, or
ip.addParameter('grid', 50);
ip.addParameter('gridH', 8);
ip.addParameter('regularizeBounds', false);
ip.addParameter('type', 'egocentric');
ip.addParameter('requiredDegreeSample', 50);
ip.addParameter('requiredSampleTime', 0.5); % Half a second
ip.addParameter('samprate', []); % Half a second
ip.parse(varargin{:})
Opt = ip.Results;

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
xy = repmat(xy, cellfun(@numel, centers(1:2)));
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

    if sum(inds)

        % Calculate aveerage firing for this bin and occupancy
        avgFR = nanmean(firing(inds));
        sumFR = nansum(firing(inds));

        % Struct for storing f(x,y) data
        xy(ix, iy).x = xcenter;
        xy(ix, iy).y = ycenter;
        xy(ix, iy).N  = sum(inds); % timespent is proportional to this
        xy(ix, iy).avgFR = avgFR;
        xy(ix, iy).sumFR = sumFR;
        xy(ix, iy).FR_occNorm = sumFR/sum(inds);

    else
        % Struct for storing f(x,y) data
        xy(ix, iy).x = xcenter;
        xy(ix, iy).y = ycenter;
        xy(ix, iy).N = 0;
        xy(ix, iy).avgFR = nan;
        xy(ix, iy).sumFR = nan;
        xy(ix, iy).FR_occNorm = nan;
    end
end

% ----------------------------
% X-Y-H position binned measures
% ----------------------------
disp("Place-HD tuning")
xyh = struct('x',nan','y',nan','h',nan,'N',0,'avgFR',nan,'sumFR',nan,'FR_occNorm',nan);
xyh = repmat(xyh, cellfun(@numel, centers));
x = discretize(X, edges{1});
y = discretize(Y, edges{2});
h = discretize(H, edges{3});
[groups, x, y, h] = findgroups(x, y, h);
uGroups = unique(groups);
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
    if sum(inds)

        % Calculate aveerage firing for this bin and occupancy
        avgFR = nanmean(firing(inds));
        sumFR = nansum(firing(inds));

        % Struct for storing f(x,y) data
        xyh(ix, iy, ih).x = xcenter;
        xyh(ix, iy, ih).y = ycenter;
        xyh(ix, iy, ih).h = hcenter;
        xyh(ix, iy, ih).N  = sum(inds); % timespent is proportional to this
        xyh(ix, iy, ih).avgFR = avgFR;
        xyh(ix, iy, ih).sumFR = sumFR;
        xyh(ix, iy, ih).FR_occNorm = sumFR/sum(inds);
    else
    end

end

% Get place field centroids
% -------------------------
xx = nd.fieldGet(xy, 'avgFR');
[~, I] = max(xx, [], 'all', 'linear');
[i,j] = ind2sub(size(xx), I);
centroid  = [centers{1}(i), centers{2}(j)];


% PRODUCE MODEL WITH UNIFIED REFERENCE POINT
% ------------------------------------------
disp('RH tuning')

% Determine which x,y positions lack more than 50 degrees of coverage
xyh = nd.apply(xyh, 'lambda', @(x) x * 1/(Opt.samprate), 'field', "totaltime-N"); % assign a total time field from N indicating the amount of time sampled
xyh = nd.apply(xyh, 'lambda', @(x) x > Opt.requiredSampleTime, 'field', "enoughtime-totaltime"); % assign a binary choice (enoughtime) from total time


% Compute heading direction stats
h_prob  = nd.fieldGet(xyh, 'N');
h_prob  = bsxfun(@rdivide, h_prob, nansum(h_prob, 3));

N_xyh = nd.fieldGet(xyh, 'N');
N_xy  = nd.fieldGet(xy,  'N');
h_set = nd.fieldGet(xyh, 'h');
x_set = nd.fieldGet(xyh, 'x');
y_set = nd.fieldGet(xyh, 'y');
h_width = 360 / Opt.gridH;
enoughtime = nd.fieldGet(xyh, 'enoughtime');
Hdim = 3;
degreeSampled = sum(enoughtime, Hdim) * h_width;
enoughdegree = degreeSampled > Opt.requiredDegreeSample;
acceptSample = bsxfun(@and, enoughdegree, enoughtime);
if Opt.acceptMethod == "and"
    acceptSample = prod(acceptSample, 3);
elseif Opt.acceptMethod == "or"
    acceptSample = sum(acceptSample, 3) > 0;
end
acceptSample = repmat(acceptSample, 1, 1, size(N_xyh, 3));

% Modulation ratio
rxyh_over_rxy = bsxfun(@rdivide, N_xyh, N_xy);
rxyh_over_rxy(~acceptSample) = nan;

%samplesWithData = any(~isnan(rxyh_over_rxy), [3]);
%allSampleWithData= all(~isnan(rxyh_over_rxy), [3]);

% Perform fminsearch, as in the paper, to find the reference points
searchFunc        = @(X) model(rxyh_over_rxy, x_set, y_set, h_set, h_prob, X(1), X(2), X(3:4));
initialConditions = [0, 0, centroid(1), centroid(2)]; % g=0, theta_p=0, ref=centroid of the place field
if Opt.regularizeBounds
    searchFunc = @(X) regularizeBounds(searchFunc(X), X(3:4), bounds);
end

% Ready fit process
disp('Fitting jercog et al. model via patternSearch mesh method');
options = optimoptions('patternsearch');
options.MeshExpansionFactor   = sqrt(2);
options.MeshContractionFactor = 1/2;
options.UseCompletePoll   = true;
options.UseCompleteSearch = true;
options.MeshTolerance = 1e-7;
%params = patternsearch(searchFunc, initialConditions, [],[],[],[],[],[],[], options);
params = fminsearch(searchFunc, initialConditions);

% Organize outputs
amplitude        = params(1)
theta_preference = mod(params(2)+pi,2*pi)-pi
reference_point  = params(3:4)

params = [amplitude, theta_preference, reference_point];


function output = model(rxyh_over_rxy, x, y, h_vals, h_prob, g, theta_preference, ref)
% The actualy nitty gritty implementation of the model

    theta_preference  = mod(theta_preference + pi, 2*pi) - pi; % 1x1

    % LAMBDAS
    % -------
    ref_angle_func = @(h, ref) atan(y-ref(2))./(x-ref(1)); % Lambda : H(x,y,h) x 1                             % Angle of animal position relative to reference point
    theta_func = @(ref_angle, h_vals) ref_angle - h_vals;  % Lambda : H(x,y,h) x 1                             % Heading relative to reference
    F_func     = @(ref_angle, h_vals, theta_preference) cos(theta_func(ref_angle, h_vals) - theta_preference); % Circular modulation by a prefered heading angle to reference
                                                                                                               % Lambda : Hx1

    ref_angle = ref_angle_func(h_vals, ref);                 % H(x,y,h) x 1
    F         = F_func(ref_angle, h_vals, theta_preference); % H(x,y,h) x 1, ranges 0 to 1


    % Should be function of x any y
    F_avg  = nanmean(h_prob .* F, 3); % expectation of F over heading directions

    % model output and imposition of model fit in jercog, Rxyh_model > 0 => abs(output)
    output = abs(1 + g .* bsxfun(@minus, F, F_avg)); % abs() because cannot have negative fr
    output = nanmean((output - rxyh_over_rxy), 'all'); % mean sqaured eerror


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

