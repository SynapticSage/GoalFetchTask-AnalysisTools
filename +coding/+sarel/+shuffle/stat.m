function out = stat(X, nestedFields, varargin)
% Which stats to pull from our main and shuffle data

ip = inputParser;
ip.addParameter('includeComponents', true);
ip.addParameter('onlyDifferences', false);
ip.addParameter('onlyShuffle', false);
ip.parse(varargin{:})
Opt = ip.Results;

main    = nd.nestedFieldCat(X, nestedFields);
out = struct();
if ~isstruct(main) && isnumeric(main)
    shuffle = nd.nestedFieldCat(X.shuffle, nestedFields);
    if Opt.onlyShuffle
        out = shuffle;
    elseif Opt.onlyDifferences
        out = main - shuffle;
    else
        out.differences = main - shuffle;
        out.meanDifferences = mean(out.differences, ndims(out.differences));
        out.varDifferences  = var(out.differences, 0,  ndims(out.differences));
        if Opt.includeComponents
            out.main    = main;
            out.shuffle = shuffle;
        end
    end
elseif isstruct(main)
    for field = string(fieldnames(main))'
        out.(field) = coding.sarel.shuffle.stat(X, [nestedFields, field], varargin{:});
    end
else
    warning('Unprocessable %s', join(nestedFields, '.'))
end
