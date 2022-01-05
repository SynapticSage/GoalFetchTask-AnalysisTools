function out = fieldMI(field, behfield, varargin)
% function entropy = fieldMI(field, varargin)
% Calculates the entropy of a field.


ip = inputParser;
ip.addParameter('shiftdim', []); % dimension of time shifts
ip.addParameter('celldim', []); % dimension of cells
ip.addParameter('outerdim', []); % outer dimension where operations are replicated over
ip.addParameter('nBins', []); % TODO
ip.addParameter('valuesPerBin', 10); % TODO
ip.parse(varargin{:})
Opt = ip.Results;

if isempty(Opt.outerdim)
    Opt.outerdim = [Opt.shiftdim, Opt.celldim];
end

% Numeric input
if isnumeric(field)

    % -----------------------
    % Fill in bins in ungiven
    % -----------------------
    if isempty(Opt.nBins) && ~isempty(Opt.valuesPerBin)
        values = numel(field);
        Opt.nBins = floor(values/Opt.valuesPerBin);
        if ~isempty(Opt.outerdim)
            Opt.nBins = Opt.nBins / size(field, Opt.outerdim);
        end
    elseif isempty(Opt.nBins)
        error("Cannot determine bins.")
    end

    % ---------------
    % Compute entropy
    % ---------------
    if isempty(Opt.outerdim)
        if Opt.nBins == 0
            out.entropy  = nan;
        else
            field(isinf(field)) = nan;

            occupancy    = behfield;
            occdistr = occupancy./sum(occupancy, 'all');
            out.occdistr    = occdistr;

            r        = squeeze(field);
            ravg     = sum( occdistr .* r );
            out.ravg = ravg;
            out.rNorm   = r./ravg;
            out.entropy = nansum( out.occdistr .* out.rNorm .* log(out.rNorm), 'all');
        end

    else

        field = num2cell(field, setdiff(1:ndims(field), Opt.outerdim));
        opt = rmfield(Opt, ["outerdim", "shiftdim", "celldim"]);
        tmp = cellfun(@(x) coding.futurepast.fieldMI(x, behfield, opt), field, "UniformOutput", true);

        entropy = nd.fieldGet(tmp, "entropy");
        entropy(isinf(entropy)) = nan;
        out.entropy = entropy;

        % Let's do some stats to be sure this is meaningful
        if ~isempty(Opt.shiftdim)
            [~, out.tauInd] = min(entropy, [], Opt.shiftdim);
            out.std = std(entropy,  Opt.shiftdim)/sqrt(numel(entropy));
            out.mu  = mean(entropy, Opt.shiftdim);
        end

    end
    
% Struct mode
elseif isstruct(field)

    out = coding.futurepast.fieldMI(field.spikeCount, behfield, Opt);
    out = util.struct.update(field, out);
    if isfield(field, 'shift')
        out.shift = field.shift(field.tauInd);
    end

else
    error("field input is wrong type")
end
