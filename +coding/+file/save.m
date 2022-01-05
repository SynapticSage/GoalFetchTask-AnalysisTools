function save(animal, index, spikes, varargin)
% function save(spikes, varargin)
% Checkpoints the save data for the spikes struct

ip = inputParser;
ip.addParameter('rawFields', ["data", "spikeTimes", "time", "binEdges", "beh", "shuffle"]);
ip.addParameter('removeAllBut', []);
ip.addParameter('append', false);
ip.addParameter('withRaw', false);
ip.parse(varargin{:})
Opt = ip.Results;

if isempty(index)

    if isfield(spikes, 'beh')
        if iscell(spikes.beh)
            spikes.beh = cat(1, spikes.beh);
        end
        index = unique(spikes.beh.day);
    else
        error("Attempted to infer days of data being saved. But no behavioral/task data in the spikes structure");
    end

end

if ~Opt.withRaw
    fields = intersect(Opt.rawFields, string(fieldnames(spikes)));
    spikes = rmfield(spikes, fields);
end
if ~isempty(Opt.removeAllBut)
    spikes = rmfield(spikes, ...
        setdiff(string(fieldnames(spikes)), string(Opt.removeAllBut)));
end

if Opt.append
    appendKws = {'-append'};
else
    appendKws = {};
end

filename_full = coding.file.filename(animal, index, varargin{:});
fprintf("\nSaving spikes to %s\n", filename_full);
save(filename_full, '-struct', 'spikes', '-v7.3', appendKws{:}); 
