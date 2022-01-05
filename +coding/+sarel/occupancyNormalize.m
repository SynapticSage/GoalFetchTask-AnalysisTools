function tuningStruct = occupancyNormalize(tuningStruct, occupancyStruct)
disp("Occcupancy normalizing results");
for field = string(fieldnames(tuningStruct))'

    if ~isnumeric(tuningStruct.(field))
        continue
    end


    % We will want to place the binning dimension last for occupacny data,
    % and create singular dimensions up to that point so that we can
    % divide one by the other
    Nt = ndims(tuningStruct.(field));
    Nc = ndims(occupancyStruct.(field));

    % Do occcuppancy normalize
    shift = Nt - Nc;
    tuningStruct.occNorm.(field) = bsxfun(@rdivide, ...
        tuningStruct.(field), ...
        shiftdim(occupancyStruct.(field), -shift));

end

