function [C1, kappa, thetahat, C2, meanRateAngle] = ...
        vonmises(rate, binnedAngle, angledim, varargin)
% Fits von mises to rate and binneed angle

if nargin < 3
    angledim = 2;
end
if angledim == -1
    angledim = ndims(rate);
end

ip = inputParser;
ip.addParameter('allowNegativeAmplitude', false); % ie allow inhibtion
ip.parse(varargin{:})
Opt = ip.Results;

% obtain mean angle
% have to broadcast binnedAngle to rate because circ_mean made a dumb design
% choice
meanRateAngle = circ_mean(ones(size(rate)) .* binnedAngle, rate, angledim);

% Setup for size manipulations later. We have an inner-product style dimension
% and a set of outer-product style dimensions. And we want to know not only
% what those are, but how long the outer dimensions are if raveled.
szOuterDims  = [size(rate,1:ndims(rate)-1), 1];
szInnerOuter = [size(rate,ndims(rate)) prod(szOuterDims)];
nOuterVars   = szInnerOuter(2);

R = gather(rate(:,:)); 
B = gather(binnedAngle(:,:));


% ----
% FIT
% ----
% Setup von mises function
if Opt.allowNegativeAmplitude
    vmfun    = @(C1, kappa, thetahat, C2, anglesample)...
        C1 * exp(abs(kappa) .* cos(thetahat - anglesample)) + C2;
else
    vmfun    = @(C1, kappa, thetahat, C2, anglesample)...
        abs(C1) * exp(abs(kappa) .* cos(thetahat - anglesample)) + C2;
end

% Find across the outer variable
opt = optimset();
param_estimates = zeros(nOuterVars, 4);
fval = zeros(nOuterVars, 1);
for i = progress(1:nOuterVars)
    anglesample = B(:, i);
    ratesample  = R(:,   i);
    if all(ratesample == 0)
        param_estimates(:,i) = nan(4, 1);
    else
        optimfun    = @(X) gather(...
            sqrt(nansum((ratesample - vmfun(X(1), X(2), X(3), X(4), anglesample)).^2))...
            );
        params = double([0, 0, meanRateAngle(i), 0]);
        [param_estimates(i,:), fval(i)] = fminsearch(optimfun, params, opt);
    end
end

% Unravel outputs
% ---------------
unravel  = @(x) reshape(x, szOuterDims);
if Opt.allowNegativeAmplitude
    C1            = unravel(param_estimates(:,1));
else
    C1            = abs(unravel(param_estimates(:,1)));
end
kappa         = unravel(abs(param_estimates(:,2))); % kappa is positve in the function only
thetahat      = unravel(mod(param_estimates(:,3)+pi,2*pi) - pi); % wrap it onto a circle
C2            = unravel(param_estimates(:,4));
meanRateAngle = unravel(meanRateAngle(:));

