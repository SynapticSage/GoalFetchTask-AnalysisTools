function [tab, summary, Opt] = scaffold2binning(tab, scaff, Opt)
% Store overall centerjcounts : The field grid!

nCenters     = arrayfun(@(x) numel(x.center), scaff);
nCentersCell = num2cell(nCenters);
binfields = string(nan(numel(scaff,1)));

% Figure out the bin per field
% -----------------------------------------------------
for ii = 1:numel(scaff)
    binfields(ii) = scaff(ii).field + "_bin";
    values = tab.(scaff(ii).field)(:,scaff(ii).col);
    %try
    %    euc = (values(:) - scaff(ii).center(:)').^2;
    %    [~, tab.(binfields(ii))] = min(euc, [], 2);
    %    clear euc
    %catch
        tab.(binfields(ii)) = util.interp.interp1(scaff(ii).center(:), 1:numel(scaff(ii).center(:)), values(:), 'nearest');
    %end
end

% Compute the overall bins for individual property bins
% -----------------------------------------------------
tab.overall_bin = zeros(height(tab), 1, 'int32');
for jj = 1:numel(scaff)
    ii = numel(scaff) - (jj - 1);
    if ii == 1
        product_term = 1;
        addition_term =  0;
    else
        product_term = prod(nCenters(1:ii-1));
        addition_term = -1;
    end
    X = table2array(tab(:, binfields(ii)));
    X = int32(X);
    tab.overall_bin = tab.overall_bin + ...
        (X + addition_term) * product_term;
end

% Create summary
% --------------
if ~isempty(tab.overall_bin)
    uBins = unique(tab.overall_bin);
    uBins = uBins(uBins > 0);
    %uBins = uBins(~isnan(uBins));
else
    uBins = [];
end
if ~isempty(uBins)
    [N, E] = histcounts(tab.overall_bin, ...
        uBins);
    summary.visits             = zeros(nCentersCell{:});
    E = double(E);
    E = E(1:end-1);
    clean = E>0;
    N = N(clean);
    E = E(clean);
    summary.visits(E) = N;
    summary.visit_time = summary.visits * Opt.samplingPeriod;
else
    warning("Empty uBins or empty tab overall_bins");
    summary.visits     = zeros(nCentersCell{:});
    summary.visit_time = zeros(nCentersCell{:});
end

