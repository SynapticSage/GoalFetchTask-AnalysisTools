function Contrast = computeContrastIndex(tuningStruct, dim, varargin)
% Computes contrast index from sarel

ip = inputParser;
ip.addParameter('nullIndex', 1); % which index in the dimension is the null value to compare against?
ip.addParameter('vmField', 'C1'); % which field for von mises to contrast?
ip.addParameter('contrastType', 'vonMises'); % which type of contrast?
ip.addParameter('Dimensions', []); % name of dimensions
ip.parse(varargin{:})
Opt = ip.Results;

if ischar(dim) || isstring(dim)
    if isempty(Opt.Dimensions)
        Opt.Dimensions = tuningStruct.Dimensions;
    end
    dim = find(strcmp(Opt.Dimensions, dim));
end
if dim <= 0
    error('Value error');
end

switch Opt.contrastType
    % Regular mode seen in the paper
    case 'vonMises'
    vals = tuningStruct.vm.(Opt.vmField);
    contrast = coding.sarel.get.contrast(vals, dim, Opt);
    Contrast.vm = contrast;
    % Compare all fields
    case 'field'
        for field = string(fieldnames(tuningStruct))'
            if isnumeric(tuningStruct.(field))
                vals = tuningStruct.(field);
                contrast = coding.sarel.get.contrast(vals, dim, Opt);
                Contrast.(field) = contrast;
            end
        end
    case 'diag_field'
        for field = string(fieldnames(tuningStruct))'
            if isnumeric(tuningStruct.(field))
                vals = diag(tuningStruct.(field));
                contrast = coding.sarel.get.contrast(vals, dim, Opt);
                Contrast.(field) = contrast;
            end
        end
    % Compare specific field
    otherwise
        vals = tuningStruct.(Opt.contrastType);
        contrast = coding.sarel.get.contrast(vals, dim, Opt);
        Contrast.(Opt.contrastType) = contrast;
end
