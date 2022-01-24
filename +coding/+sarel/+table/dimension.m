function D = dimension(raw, Binning, dimensions, tuning_curve)
%

if nargin < 4
    tuning_curve = "";
end
if isstruct(raw)
    sizes = structfun(@(x) numel(x), raw, 'UniformOutput', false);
    [fieldmax, max_] = deal("", 0);
    for field = string(fieldnames(sizes))'
        if sizes.(field) > max_; fieldmax = field; end
    end
    raw = raw.(fieldmax);
end

D = struct();
d = 0;
for dim = dimensions
    d = d + 1;
    tmp = 1:size(raw, d);
    tmp = shiftdim(tmp(:), -(d-1));
    D.("d_" + dim) = tmp;
    % dimension has an alias?
    if dim == "bins" && isfield(Binning.centers, tuning_curve)
        D.("bins") = reshape(Binning.centers.(tuning_curve), size(D.("d_" + dim)));
    end
end
