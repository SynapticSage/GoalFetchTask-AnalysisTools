function T = modelreal(Struct, varargin)
% Outputs a table comparing real values to model values

ip = inputParser;
ip.addParameter('useGPU', true);
ip.parse(varargin{:})
Opt = ip.Results;

inds      = nd.indicesMatrixForm(Struct);
marginals = coding.jercog.table.findmarginals(nd.get(Struct, inds(1,:)));
T = cell(size(inds,1), 1);

keyboard
count = 0;
for ind = inds'
    
    count = count + 1;

    one = nd.get(Struct, ind);

    % GRAB MODEL LENGTH VARIABLES
    % ---------------------------
    rxyh = one.xyh.FR_occNorm;
    rxyh_over_rxy = one.rxyh_over_rxy;
    [nNeurons, nX, nY, nH] = size(rxyh);
    
    % GRAB LESS THAN MODEL LENGTH
    % ---------------------------
    ravel = @(x) x(:);
    shape = @(x) ravel(repmat(x, 1, 1, 1, nH));
    rxy   = shape(one.xy.FR_occNorm);

    % GRAB NEURON LENGTH VARIABLES
    % ----------------------------
    shape             = @(x) ravel(repmat(x, [1, nX, nY, nH]));
    neurons           = shape(1:size(one.amplitude, 1));
    centroid_x        = shape(one.centroid(:, 1));
    centroid_y        = shape(one.centroid(:, 2));
    reference_point_x = shape(one.reference_point(:, 1));
    reference_point_y = shape(one.reference_point(:, 2));
    theta_preference  = shape(one.theta_preference(:, 1));
    amplitude         = shape(one.amplitude(:, 1));

    % GRAB DIMENSIONS
    % ---------------
    shape = @(x) ravel(repmat(x, [nNeurons, 1, nY, nH]));
    x = shape(one.x);
    shape = @(y) ravel(repmat(y, [nNeurons, nX, 1, nH]));
    y = shape(one.y);
    shape = @(h) ravel(repmat(h, [nNeurons, nX, nY, 1]));
    h = shape(one.h);

    % Let's ravel our other varaibles now
    rxyh = ravel(rxyh);

    % Model output
    if isfield(one, 'model_output')
        model_output = ravel(one.model_output);
    else % it's not store, we have to compute it!
        if Opt.useGPU
            [goodness_of_fit, model_output] = coding.jercog.gpumodel(rxyh_over_rxy, x, y, h_vals, h_prob, g, theta_preference, ref);
        else
            [goodness_of_fit, model_output] = coding.jercog.model(rxyh_over_rxy, x, y, h_vals, h_prob, g, theta_preference, ref);
        end
    end
    goodness_of_fit = repmat(one.goodness_of_fit, [1, nX, nY, nH]);
    goodness_of_fit = ravel(goodness_of_fit);

    rxyh_over_rxy = ravel(rxyh_over_rxy);

    % Place what we've learned into the table
    tab = table(neurons, x, y, h, centroid_x, centroid_y, reference_point_x, reference_point_y, theta_preference, amplitude, rxyh, rxy, rxyh_over_rxy, model_output);

    % ADD MARGINALS
    for marginal = marginals
        tab.(marginal) = repmat(one.(marginal), height(tab), 1);
    end

    T{count} = tab;
end

T = cat(1, T{:});
