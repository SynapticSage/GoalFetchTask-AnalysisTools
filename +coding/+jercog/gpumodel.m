function [goodness_of_fit, evaluation] = gpumodel(rxyh_over_rxy, x, y, h_vals, h_prob, g, theta_preference, ref)
% The actualy nitty gritty implementation of the model
    
    debug = false;

    rxyh_over_rxy    = gpuArray(rxyh_over_rxy);
    x                = gpuArray(x);
    y                = gpuArray(y);
    h_vals           = gpuArray(h_vals);
    h_prob           = gpuArray(h_prob);
    g                = gpuArray(g);
    theta_preference = gpuArray(theta_preference);
    ref              = gpuArray(ref);
    if debug
        size(rxyh_over_rxy),
        size(x),
        size(y),
        size(h_vals),
        size(h_prob),
        size(g),
        size(theta_preference),
        size(ref)
    end
    % -------
    % LAMBDAS
    % -------
    bsxminus   = @(a,b) bsxfun(@minus, a, b);
    ref_angle_func = @(h, ref) atan(bsxminus(y,ref(:,2)) ...
        ./bsxminus(x,ref(:,1))); % Lambda : H(x,y,h) x 1                             % Angle of animal position relative to reference point
    theta_func = @(ref_angle, h_vals) bsxminus(ref_angle, h_vals);  % Lambda : H(x,y,h) x 1                             % Heading relative to reference
    F_func     = @(ref_angle, h_vals, theta_preference) cos(bsxminus(theta_func(ref_angle, h_vals), theta_preference)); % Circular modulation by a prefered heading angle to reference

    % Get macro quantities
    theta_preference  = mod(theta_preference + pi, 2*pi) - pi; % 1x1
    ref_angle = ref_angle_func(h_vals, ref);                 % nNeurons x H(x,y,h)
    F         = F_func(ref_angle, h_vals, theta_preference); % H(x,y,h) x 1, ranges 0 to 1


    % Should be function of X any Y
    headingDimension = ndims(rxyh_over_rxy);
    F_avg            = nanmean(h_prob .* F, headingDimension); % expectation of F over heading directions, measures the average circular modulation

    % model goodness_of_fit and imposition of model fit in jercog, Rxyh_model > 0 => abs(goodness_of_fit)
    evaluation = abs(1 + (g .* bsxfun(@minus, F, F_avg))); % abs() because cannot have negative fr
    goodness_of_fit = nanmean((evaluation - rxyh_over_rxy).^2, 'all'); % mean sqaured error

    % Bring back from the GPU
    evaluation      = gather(evaluation);
    goodness_of_fit = gather(goodness_of_fit);
    %if isnan(goodness_of_fit )
    %    goodness_of_fit = inf;
    %end

