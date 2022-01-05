function fieldStruct = combinefields(fieldStructs)
% function fieldStruct = combinefields(fieldStructs)
%
% Combines all of the individual field structures into one field sstructure

fieldStruct = struct();
for prop = string(fieldnames(fieldStructs))'
    if isnumeric(fieldStructs(1).(prop))
        values = arrayfun(@(x) shiftdim(x.(prop), -1), fieldStructs, 'UniformOutput', false);
        fieldStruct.(prop) = cat(1, values{:});
    else
        values = arrayfun(@(x) x.(prop), fieldStructs, 'UniformOutput', false);
        fieldStruct.(prop) = values;
    end
end

