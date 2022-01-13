function [out, Binning] = main(spikes, beh, varargin)
% Implements the sarel paper tuning curves
% and rayleigh scoring

ip = inputParser;
ip.KeepUnmatched = true;
ip.addParameter('useGPU', true);
ip.addParameter('splitSelect', ["stops", "starts", "startstops",...
    "cuememstops", "cuememstarts", "cuememstartstops"]);
ip.addParameter('split_by_name', []);
ip.addParameter('split_by', {});
ip.addParameter('add_maginal_to_all_splits', []);
ip.parse(varargin{:})
Opt = ip.Results;

% If spikes.beh are indices(beh) for spikeTiems instead of actual behavior
% variables for spike time
if spikes.behtype == "indices"
    spikes.beh = units.indexToBehavior(spikes.beh, beh);
end

[Binning, optNew] = coding.sarel.binning(beh, varargin{:});
Opt = util.struct.update(Opt, optNew);

out.occupancy = coding.sarel.tuning(beh, [], Binning, Opt); % Returns behavior occupacny
out.ground    = coding.sarel.tuning(spikes, [], Binning, Opt); % Returns behavior occupacny
Opt.occupancy = out.occupancy;
Opt.ground    = out.ground;

%%%%%%%%%%%%%%%%%%%% TUNING CURVES per neuron %%%%%%%%%%%%%%%%%%%%%%%%
%% UNCONDITIONAL :  Capture distribution for a given neuron to any well
% and it's distribution of spiking to different start and end conditions

%%%%%%%%%%%%%%%%%%%% TUNING CURVES per neuron marginalized %%%%%%%%%%%
%% CONDITIONAL :  Capture the distributions when conditioned on well
% or maze locations
if isempty(Opt.split_by_name)
    Opt.split_by_name = ["tuning", "stops", "starts",...
        "startstops", "cuememstops", "cuememstarts", "cuememstartstops"];
end
if isempty(Opt.split_by)
    Opt.split_by = {"neuron",...
        ["neuron", "stopWell"], ...
        ["neuron", "startWell"], ...
        ["neuron", "stopWell", "startWell"], ...
        ["neuron", "cuemem", "stopWell"],...
        ["neuron", "cuemem", "startWell"], ...
        ["neuron", "cuemem", "stopWell", "startWell"]};
end
if Opt.useGPU
    workflow_func = @coding.sarel.helper.GPUworkflow;
else
    workflow_func = @coding.sarel.helper.workflow;
end

for i = 1:numel(Opt.split_by)

    split      = Opt.split_by{i};
    split_name = Opt.split_by_name(i);
    if ~ismember(split_name, Opt.splitSelect)
        disp("Skipping " + join(split_name, " "))
        keyboard
        continue
    end
    out.(split_name) = workflow_func(spikes, split, Binning, Opt);

end

out.Binning = Binning;
