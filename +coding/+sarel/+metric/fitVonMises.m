function outStruct = fitVonMises(tuningStruct, angles, angledim, contain_fields, varargin)
% function outStruct = computeDirectionalityIndexAndThetaHat(tuningStruct, angles, angledim, contain_fields)
% 
% Computes directionality index and the theta hat metric

ip = inputParser;
ip.addParameter('allowNegativeAmplitude', false);
ip.parse(varargin{:})
Opt = ip.Results;

if isfield(tuningStruct, 'occNorm')
    tuningStruct = tuningStruct.occNorm;
end

disp("Rayleigh tuning scores");
outStruct = struct();
for field = string(fieldnames(tuningStruct))'

    if ~all(arrayfun(@(cf) contains(field, cf), contain_fields))
        continue
    end

    [outStruct.(field).C1, ...
    outStruct.(field).kappa, ...
    outStruct.(field).thetahat, ...
    outStruct.(field).C2, ...
        outStruct.(field).meanRateAngle] = ...
        coding.sarel.computation.vonmises(tuningStruct.(field), angles, angledim, Opt);
end

