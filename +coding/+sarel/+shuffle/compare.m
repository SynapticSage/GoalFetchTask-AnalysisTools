function comparison = compare(sarel, varargin)
% compares shuffle to main data in the sarel struct

ip = inputParser;
ip.addParameter('method', @minus);
ip.addParameter('measures', ["vm", "occNorm", "rayleigh", "maxmean_indices"]);
ip.parse(varargin{:})
Opt = ip.Results;

first = struct('type', '()', 'subs', {{[1]}});
datafields = @(x) util.struct.matchingfields(x,...
    'logic', 'or',...
    'keyConditions', {@(y) subsref(isstrprop(y, 'lower'), first)},...
    'valueConditions', {@(y) ~isstruct(y)});

measures = @(x) util.struct.matchingfields(x,...
    'logic', 'or',...
    'keyConditions', {@(y) all(isstrprop(y, 'lower'))},...
    'valueConditions', {@(y) isstruct(y)});

comparison = struct();

% Iterate the splits
splits = datafields(sarel);
for split = splits

    S = sarel.(split);
    metrics = intersect(measures(S), Opt.measures);

    for datafield = datafields(S)
        compare.(datafield) = shufstat(sarel, [split, datafield]);
        for metric = metrics
            compare.(metric).(datafield) = shufstat(sarel, [split, metric, datafield]);
        end
    end

end
