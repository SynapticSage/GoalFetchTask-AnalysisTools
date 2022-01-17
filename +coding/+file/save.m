function save(animal, index, spikes, varargin)
% function save(spikes, varargin)
% Checkpoints the save data for the spikes struct

ip = inputParser;
ip.addParameter('rawFields', ["data", "spikeTimes", "time", "binEdges", "beh", "shuffle"]);
ip.addParameter('addRawField', []);
ip.addParameter('removeAllBut', []);
ip.addParameter('append', false);
ip.addParameter('withRaw', false);
ip.addParameter('filename_full', []);
ip.addParameter('checkpointActive', {}, @iscell);
ip.parse(varargin{:})
Opt = ip.Results;

% Do not save if this isn't an active checkpoint
if ~coding.file.checkpointActive(Opt.checkpointActive{:})
    return;
else
    disp('Checkpointing...');
end

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
    if Opt.addRawField % put a set of raw fields back into the set of fields
        fields = setdiff(fields, Opt.addRawField);
    end
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

if isempty(Opt.filename_full)
    filename_full = coding.file.filename(animal, index, varargin{:});
else
    filename_full = Opt.filename_full;
end
fprintf("\nSaving spikes to %s\n", filename_full);
save(filename_full, '-struct', 'spikes', '-v7.3', appendKws{:}); 
