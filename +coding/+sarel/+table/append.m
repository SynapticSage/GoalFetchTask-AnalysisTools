function tab = append(previous_tab, Metrics, Description, varargin)
% Appends to a split table
%
% see coding.sarel.table()
%
% Inputs
% -------
%
% tab : table
%
% Metrics : struct
%   fields are metric_stat data columns of coding.sarel.table
%
% Description : struct
%   fields are constants on all of the rows
%   with the exception of a special field called
%   .Dimension
%   
%   .Dimension - each field is a dimension meter of
%                   the same sizez as Metrics fields
%                   telling us the meaning of the tensor
%                   grid of these variables

ip = inputParser;
ip.addParameter('includeMetricName', false);
ip.parse(varargin{:})
Opt = ip.Results;

for d = string(fieldnames(Description.Dimensions))'
    tab.(d) = Description.Dimensions.(d);
end

for d = string(setdiff(fieldnames(Description), "Dimensions"))'
    if ischar(Description.(d)) || iscellstr(Description.(d))
        tab.(d) = string(Description.(d));
    else
        tab.(d) = Description.(d);
    end
end


if Opt.includeMetricName
    Metrics = nd.fullyUnnest(Metrics);
end
for metric = string(fieldnames(Metrics))'
    if ~isstruct(Metrics.(metric))
        if isnumeric(Metrics.(metric))
            tab.(metric) = Metrics.(metric);
        elseif ischar(Metrics.(metric)) || iscellstr(Metrics.(metric))
            tab.(metric) = string(Metrics.(metric));
        end
    else
        for field = string(fieldnames(Metrics.(metric)))'
            tab.(metric + "_" + field) = Metrics.(metric).(field);
        end
    end
end

% broadcast the fields (in a pythonic sense)
temp = tab;
tab = nd.broadcast(tab);
% and apply the ravel function
tab = structfun(@(x) x(:), tab, 'UniformOutput', false);
tab = struct2table(tab);
tab = [previous_tab; tab];
