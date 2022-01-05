function [C1, kappa, thetahat, C2, meanRateAngle] = ...
        vonmises(rate, binnedAngle, angledim)
% Fits von mises to rate and binneed angle

    if nargin < 3
        angledim = 2;
    end
    if angledim == -1
        angledim = ndims(rate);
    end

    meanRateAngle = circ_mean(binnedAngle, rate, angledim);

    szOuterDims   = [size(rate,2:ndims(rate)), 1];
    szInnerOuter = [size(rate,1) prod(szOuterDims)];
    zeromat      = zeros(szInnerOuter(2), 1);
    nOuterVars   = szInnerOuter(2);

    % Fit a von mises function
    R = gather(rate(:,:));
    B = gather(binnedAngle(:,:));
    vmfun    = @(C1, kappa, thetahat, C2, anglesample)...
        C1 * exp(abs(kappa) .* cos(thetahat - anglesample)) + C2;
    opt = optimset();
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
            [param_estimates(:,i), fval(i)] = fminsearch(optimfun, params, opt);
        end
    end

    % Unravel outputs
    % ---------------
    unravel  = @(x) reshape(x, szOuterDims);
    C1       = unravel(param_estimates(1,:));
    kappa    = unravel(abs(param_estimates(2,:))); % kappa is positve in the function only
    thetahat = unravel(mod(param_estimates(3,:)+pi,2*pi) - pi); % wrap it onto a circle
    C2       = unravel(param_estimates(4,:));

    % Discover the estimated maximal angle
    szRate = size(rate);
    binnedAngle = bsxfun(@times, binnedAngle, ones(szRate));

    permuteDims = 1:ndims(rate); % Bring the inner dimension to the first position, and outers to the last two positions
    permuteDims = circshift(permuteDims, 1);

    rate     = permute(rate, permuteDims);
    binnedAngle   = permute(binnedAngle, permuteDims);
    angledim = find(permuteDims == angledim);

    meanRateAngle = circ_mean(binnedAngle, rate, angledim);

    szOuterDims   = [size(rate,2:ndims(rate)), 1];
    szInnerOuter = [size(rate,1) prod(szOuterDims)];
    zeromat      = zeros(szInnerOuter(2), 1);
    nOuterVars   = szInnerOuter(2);

    % Fit a von mises function
    R = gather(rate(:,:));
    B = gather(binnedAngle(:,:));
    vmfun    = @(C1, kappa, thetahat, C2, anglesample)...
        C1 * exp(abs(kappa) .* cos(thetahat - anglesample)) + C2;
    opt = optimset();
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
            [param_estimates(:,i), fval(i)] = fminsearch(optimfun, params, opt);
        end
    end

    % Unravel outputs
    % ---------------
    unravel  = @(x) reshape(x, szOuterDims);
    C1            = unravel(param_estimates(1,:));
    kappa         = unravel(abs(param_estimates(2,:))); % kappa is positve in the function only
    thetahat      = unravel(mod(param_estimates(3,:)+pi,2*pi) - pi); % wrap it onto a circle
    C2            = unravel(param_estimates(4,:));
    meanRateAngle = unravel(meanRateAngle(:));

