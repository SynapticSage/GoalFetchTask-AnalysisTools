function [Y, run_stats] = equalSamplesDims(X, reorderDims, stackDims, nonstackDims, varargin)
% Other critical functoin that equalizes number of samples in the equalizeDimensions across preserved dimensions.
%
% issue
%
% stackDims breakdown into unstacked and stacked preserve dqequaqequaims

ip = inputParser;
ip.addParameter('ploton', false);
ip.parse(varargin{:})
Opt = ip.Results;


if nargin < 4
    nonstackDims = size(X, ndims(X));
end

% Solve for missing dimensions
if isempty(reorderDims)
    reorderDims = setdiff(1:ndims(X), [stackDims, nonstackDims]);
elseif isempty(stackDims)
    stackDims = setdiff(1:ndims(X), [reorderDims, nonstackDims]);
elseif isempty(nonstackDims)
    nonstackDims = setdiff(1:ndims(X), [reorderDims, stackDims]);
end

szX = size(X);
%c.reorderDims  = num2cell(reorderDims);
%c.stackDims    = num2cell(stackDims);
%c.nonstackDims = num2cell(nonstackDims);
%[~,D] = util.num.dimTensor(X);

% for eacch of the preserve dims, figure out which of the eequalization dims are full samples (where all values are not nan)
filled      = any(~isnan(X), nonstackDims);
validcount  = sum(filled, reorderDims);
commoncount = min(validcount, [], reorderDims); % smallest number of equal across all conditions
commoncountoverall = min(commoncount(:));
if Opt.ploton
    fig commoncount;
    util.plot.tile(@imagesc, 3, commoncount, 'limit', 100);
end
if commoncountoverall == 0
    partials      = mean(~isnan(X), nonstackDims);
    fig partials;
    util.plot.tile(@imagesc, reorderDims, partials, 'limit', 100);
    error("Empty");
end
run_stats.filled = filled;
run_stats.validcount = validcount;
run_stats.commoncount = commoncount;

%fillInds = D;
%for r = reorderDims
%    fillInds{r} = fillInds{r}(filled);
%end

goodSamples = zeros(szX(stackDims));
inds = util.indicesMatrixForm(size(goodSamples));
szY = szX;
szY(reorderDims) = commoncountoverall;
Y = nan(szY);
for ind = progress(inds','Title','Reordering')

    getSubset = num2cell(repmat(':', 1, ndims(X)));
    getSubset(stackDims) = num2cell(ind);

    % Figure our which of the equalization dims are valid
    v = validcount(getSubset{:});
    v = randperm(v);
    v = v(1:commoncount);

    % Get coordinates in the reorder dim we will grab TODO
    samples = find(filled(getSubset{:}));
    samples = samples(v);

    setSubset = getSubset;
    getSubset{reorderDims} = v;
    Y(setSubset{:}) = X(getSubset{:});

end
