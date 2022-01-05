function T = params(Struct, varargin)
% Returns a params table

ip = inputParser;
ip.addParameter('addCellTable', []);
ip.parse(varargin{:})
Opt = ip.Results;

inds      = nd.indicesMatrixForm(Struct);
marginals = coding.jercog.table.findmarginals(nd.get(Struct, inds(1,:)));
disp("Marginals = " + marginals);
T = cell(size(inds,1), 1);

ii = 0;
for ind = inds'
    ii = ii + 1;

    one = nd.get(Struct, ind);


    ravel = @(x) x(:);

    if isempty(one.amplitude)
        continue
    end

    neurons           = ravel(1:size(one.amplitude, 1));
    centroid_x        = ravel(one.centroid(:, 1));
    centroid_y        = ravel(one.centroid(:, 2));
    reference_point_x = ravel(one.reference_point(:, 1));
    reference_point_y = ravel(one.reference_point(:, 2));
    theta_preference  = ravel(one.theta_preference(:, 1));
    amplitude         = ravel(one.amplitude(:, 1));


    tab = table(neurons, centroid_x, centroid_y, reference_point_x, reference_point_y, theta_preference, amplitude);
    if isfield(one, 'acceptedSamples') && ~isempty(one.acceptedSamples)
        tab.acceptedSamples      = repmat(one.acceptedSamples      , height(tab) , 1);
        tab.fracEnoughTime       = repmat(one.fracEnoughTime       , height(tab) , 1);
        tab.fracAcceptableDegree = repmat(one.fracAcceptableDegree , height(tab) , 1);
    end

    % Let's grab the distance from the place center to reference center
    tab.ref_to_cent_x = tab.reference_point_x - tab.centroid_x;
    tab.ref_to_cent_y = tab.reference_point_y - tab.centroid_y;
    tab.ref_to_cent   = sqrt(tab.ref_to_cent_x.^2 + tab.ref_to_cent_y.^2);

    for marginal = marginals
        if isempty(one.(marginal))
            one.(marginal) = -1;
        end
        tab.(marginal) = repmat(one.(marginal), height(tab), 1);
    end

    if ~isempty(Opt.addCellTable)
        if height(tab) ~= height(Opt.addCellTable)
            warning('tab height not addCellTable height')
            continue
        end
        tab = [tab, Opt.addCellTable];
    end

    T{ii} = tab;
end

T = util.cell.icat(T,'removeEmpty',true,'fieldCombine','union');
