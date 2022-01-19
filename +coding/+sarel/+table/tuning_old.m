function data = tuning(Struct, varargin)
% Outputs a table comparing real values to model values

ip = inputParser;
ip.addParameter('useGPU', true);
ip.addParameter('valueTag', "value");
ip.addParameter('nestedMarginalMeasurements', ["occNorm", "rayleigh"]);
ip.addParameter('addMargTag', true);
ip.addParameter('addCellTable', []);
ip.parse(varargin{:})
Opt = ip.Results;

inds        = nd.indicesMatrixForm(Struct);
T           = cell(size(inds,1), 1);
nonNeuronal = ["occupancy", "ground"];
Binning = Struct.Binning;

for margStName = getFieldnames(Struct)

    disp(margStName);

    % Grab the field and geet its marginals
    marginalStruct = Struct.(margStName);
    marginals = coding.sarel.table.field.findmarginals(margStName);
    if isempty(marginals)
        continue
    end

    % Grab data table
    data.(margStName) = marginalStructTable(marginalStruct, marginals, Binning, Opt);

end


%---------- Helper functions ------------------------------------------
function out = marginalStructTable(marginalStruct, marginals, Binning, Opt)

    out = struct();
    ravel = @(x) x(:);
    for hf = getFieldnames(marginalStruct)
        if isstruct(marginalStruct.(hf))
            continue
        end
        data = getBinSize(marginalStruct.(hf), marginals, hf, Binning, Opt);
        data.(Opt.valueTag) = marginalStruct.(hf);
        data.(Opt.valueTag + "_occNorm") = marginalStruct.occNorm.(hf);
        if isfield(marginalStruct.rayleigh, hf)
            for f = getFieldnames(marginalStruct.rayleigh.(hf))
                data.(Opt.valueTag + "_" + f) = marginalStruct.rayleigh.(hf).(f);
            end
        end
        data = nd.simplebroadcast(data);
        data = structfun(@(x) ravel(x), data, "UniformOutput", false);
        data = struct2table(data);
        if ~isempty(Opt.addCellTable)
            data = [data, Opt.addCellTable(data.neuron,:)];
        end
        out.(hf) = data;
    end


function out = getBinSize(item, marginals, hf, Binning, Opt)
    % getBinSize grabs the size of the indivudal bins

    %Setup a mapping from field to it's relevent binning dimension
    histfield = containers.Map();
    histfield('stopWell')        = 'possibleStops';
    histfield('startWell')       = 'possibleStarts';
    histfield('currentDistance') = 'distCenters';
    histfield('currentAngle')    = 'angleCenters';
    histfield('cuemem')          = 'possibleCuemem';
    marginalfield = containers.Map();
    marginalfield('stops')  = 'possibleStops';
    marginalfield('starts') = 'possibleStarts';
    marginalfield('cuemem') = 'possibleCuemem';
    marginalfield('stop')   = 'possibleStops';
    marginalfield('start')  = 'possibleStarts';

    count = 0;
    szItem = size(item);
    endItem = ndims(item);
    for s = szItem
        count = count + 1;
        if count == 1
            out.neuron = (1:s)';
        elseif count == endItem
            out.(hf) = shiftdim(Binning.(histfield(hf)), -count);
        else
            marginal = marginals(count-1);
            if Opt.addMargTag
                tag = "marg_";
            else
                tag = "";
            end
            out.(tag + marginal) = shiftdim(Binning.(marginalfield(marginal)), -count);
        end
    end

function out = getFieldnames(Struct, Opt)
    out = setdiff(string(fieldnames(Struct)), ["Binning", "findgroups"])';
