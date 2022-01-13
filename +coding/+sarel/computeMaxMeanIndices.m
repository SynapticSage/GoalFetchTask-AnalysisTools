function out = computeMaxMeanIndices(tuningStruct, Binning)
% Computes distance and angle indices from
%
% These can be seen in the supplementary sections of
% Sarel et al. 2017
%
% Indices : Max(firing) - Mean(firing)
%
% Using to define the euclidean and path distance
% indices

fields = coding.sarel.table.field.standardBin();

% Each of the fields
for field = fields(:)'

    out.(field + "_index") = ...
        max(tuningStruct.(field), [], ndims(tuningStruct.(field))) ...
        ./ mean(tuningStruct.(field), ndims(tuningStruct.(field)));
    
end
