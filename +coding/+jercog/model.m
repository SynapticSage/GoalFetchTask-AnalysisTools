function [goodness_of_fit, evaluation] = model(rxyh_over_rxy, x, y, h_vals, h_prob, g, theta_preference, ref)
% The actualy nitty gritty implementation of the model

%global num
%ref = ref / 1000;

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


% Should be function of x any y
F_avg  = nanmean(h_prob .* F, ndims(rxyh_over_rxy)); % expectation of F over heading directions

% model goodness_of_fit and imposition of model fit in jercog, Rxyh_model > 0 => abs(goodness_of_fit)
evaluation = abs(1 + (g .* bsxfun(@minus, F, F_avg))); % abs() because cannot have negative fr
goodness_of_fit = nanmean((evaluation - rxyh_over_rxy).^2, 'all'); % mean sqaured eerror
