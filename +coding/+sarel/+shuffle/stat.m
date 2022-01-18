function out = stat(X, nestedFields, varargin)
% Which stats to pull from our main and shuffle data

ip = inputParser;
ip.addParameter('includeComponents', true);
ip.parse(varargin{:})
Opt = ip.Results;

main    = nd.nestedFieldCat(X, nestedFields);
shuffle = nd.nestedFieldCat(X.shuffle, ["shuffle", nestedFields]);
out.differences = main - shuffle;
out.meanDifferences = mean(differences, ndims(differences));
out.varDifferences  = var(differences,  ndims(differences));
if Opt.includeComponents
    out.main = main;
    out.shuffle = shuffle;
end
