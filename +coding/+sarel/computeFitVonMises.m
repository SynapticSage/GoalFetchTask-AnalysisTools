function outStruct = computeVonMises(tuningStruct, angles, angledim, contain_fields, varargin)
% function outStruct = computeDirectionalityIndexAndThetaHat(tuningStruct, angles, angledim, contain_fields)
% 
% Computes directionality index and the theta hat metric

ip = inputParser;
ip.addParameter('allowNegativeAmplitude', false);
ip.parse(varargin{:})
Opt = ip.Results;

disp("Rayleigh tuning scores");
outStruct = struct();
for field = string(fieldnames(tuningStruct))'

    if ~all(arrayfun(@(cf) contains(field, cf), contain_fields))
        continue
    end

    [outStruct.C1, ...
    outStruct.kappa, ...
    outStruct.thetahat, ...
    outStruct.C2, ...
        outStruct.meanRateAngle] = ...
        coding.sarel.get.vonmises(tuningStruct.(field), angles, angledim, Opt);
end

