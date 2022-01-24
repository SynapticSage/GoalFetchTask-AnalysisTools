function tab = append(previous_tab, Metrics, Description)
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



for d = fieldnames(Description.Dimensions)
    tab.(d) = Description.Dimensions.(d);
end

for d = setdiff(fieldnames(Description), "Dimensions")'
    tab.(d) = Description.(d);
end

for metric = string(fieldnames(Metrics)')
    if ~isstruct(Metrics.(metric))
        tab.(metric) = Metrics.(metric);
    else
        for field = string(fieldnames(tab.(metric)))'
            tab.(metric + "_" + field) = Metrics.(metric);
        end
    end
end

% broadcast the fields (in a pythonic sense)
tab = nd.broadcast(tab);
% and apply the ravel function
tab = nd.apply(tab, @(x) x(:));
tab = struct2table(tab);
tab = [previous_tab; tab];