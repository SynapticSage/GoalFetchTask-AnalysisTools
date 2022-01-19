function D = dimension(raw, Binning, dimensions)
%

D = struct();
d = 0;
for dim = Dimensions
    d = d + 1;
    D.("d" + dim) = 1:size(raw, d);
    % dimension has an alias?
    if isfield(Binning.center.(dim))
        D.(dim) = Binning.center.(dim);
    end
end
