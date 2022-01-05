function out = tuning(spikes, fields, Binning, varargin)
% Implements grabbing tunning curves for some
% marginalization scheme determined by fields

ip = inputParser;
ip.KeepUnmatched = true;
ip.addParameter('anglebins', 18);
ip.addParameter('distbins',  20);
ip.addParameter('anglefield', 'currentAngle');
ip.addParameter('distancefield', 'currentEucDist');
ip.parse(varargin{:});
Opt = ip.Results;

fields = string(fields); 
if isempty(fields)
    kws = {};
else
    kws = {'Title', char("fields = " + join(fields, "-"))};
end

if istable(spikes)
    where = "spikes.";
else
    where = "spikes.beh.";
end

if ~isempty(fields)
    gs = util.table.findgroups(spikes.beh, fields);
else
    gs.uGroups = 1;
end

for g = progress(gs.uGroups', kws{:})

    if ~isempty(fields)
        filtstr = replace(gs.evalstring(g),"$", where);
        filt    = eval(filtstr);
        subset  = spikes.beh(filt,:);
        address = num2cell(gs.group.addressByGroupNum{g});
    elseif istable(spikes)
        subset = spikes;
        address = {1};
    else
        subset = spikes.beh;
        address = {1};
    end

    [N, ~] = histcounts(subset.(Opt.anglefield),    Binning.angleEdges);
    currentAngle(address{:},:) = N;

    [N, ~] = histcounts(subset.(Opt.distancefield), Binning.distEdges); %TODO examine if these bin edges are correect
    currentDistance(address{:},:) = N;

    edges = union(Binning.possibleStarts, max(Binning.possibleStarts) + 1);
    [N, ~] = histcounts(subset.stopWell,  edges);
    stopWell(address{:},:) = N;

    [N, ~] = histcounts(subset.startWell, union(Binning.possibleStarts, max(Binning.possibleStarts) + 1));
    startWell(address{:},:) = N;

    [N, ~] = histcounts(subset.cuemem, union(Binning.possibleCuemem, max(Binning.possibleCuemem)+1));
    cuemem(address{:},:) = N;

    % Add histo to all points here!

end
disp('');
disp('');
disp('');


out = struct();
out.stopWell        = stopWell;
out.startWell       = startWell;
out.currentDistance = currentDistance;
out.currentAngle    = currentAngle;
out.cuemem          = cuemem;
if numel(gs.uGroups) > 1
    out = util.struct.appendgroupfield(out, gs);
end
