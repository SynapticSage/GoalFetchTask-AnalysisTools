function spikes = jercog(spikes, beh, varargin)
%function spikes = jercog(spikes, beh, varargin)
% 
% Runs jercog analysis

ip = inputParser;
ip.KeepUnmatched = true;
ip.addParameter('name', 'jercog');
ip.addParameter('parallel', false);
ip.addParameter('checkpointVar', []);
ip.parse(varargin{:})
Opt = ip.Results;


% Jercog Kandel and Abott Paper
tic
spikes.(Opt.name) = coding.jercog.allNeurons_sparse_lean(spikes, beh, 'useGPU', true);
toc

if Opt.checkpointVar
    disp("Checkpointing");
    assignin('base', Opt.checkpointVar, spikes.(Opt.name));
end

clear job

% Condition-wise
conditionwise_kws = {'useGPU_calcfield', true, 'useGPU', true, 'grid', 10, 'gridH', 6, 'acceptMethod', 'and'};
gs = util.table.findgroups(beh, ["stopWell", "correct"]);
address = cellfun(@num2cell, num2cell(cat(2,gs.group.address{:}), 2), 'UniformOutput', false);
if Opt.parallel
    cluster = parcluster;
end
spikes.(Opt.name).byGoal = struct();
for g = progress(setdiff(gs.uGroups(3:end),2)', 'Title', 'Submitting jercog conditions')
    disp(g)
    if exist('job','var') && numel(job) >= g && ~isempty(job{g}) 
        continue
    end
    filtstr = "$" + gs.conditionLabels' + " == " + cellfun(@(x) string(x(g)), gs.group.values)';
    clear sub
    [sub.spikes, sub.beh] = tabfilt.spikesbeh(spikes, beh, filtstr);
    if Opt.parallel
        util.job.waitQueueSize(4, 20, cluster);
        job{g} = batch(@coding.jercog.allNeurons_sparse_lean, 1, ...
            {sub.spikes, sub.beh, conditionwise_kws{:}});
    else
        try
            sub = coding.jercog.allNeurons_sparse_lean(sub.spikes, sub.beh, conditionwise_kws{:});
            job{g} = sub;
            sub = util.struct.update(sub, gs.group.field, @(x) x(g)); % Add field information
            spikes.(Opt.name).byGoal = nd.set(spikes.(Opt.name).byGoal, address{g}, sub);
        catch ME
            if contains(ME.message,"no accepted samples")
                job{g} = struct();
            else
                throw(ME)
            end
        end
    end

    if Opt.checkpointVar
        assignin('base', Opt.checkpointVar, spikes.jercog);
        assignin('base', 'g', g);
    end
end

if Opt.parallel
    for g = progress(gs.uGroups(end:-1:1)', 'Title', 'Fetching outputs')
        tmp = fetchOutputs(job{g});
        sub.spikes = tmp{1};
        sub.spikes.(Opt.name) = util.struct.update(sub.spikes.(Opt.name), gs.group.field, @(x) x(g)); % Add field information
        spikes.(Opt.name).byGoal(address{g}{:}) = sub.spikes.(Opt.name);
    end

    for g = progress(gs.uGroups(end:-1:1)', 'Title', 'Deleting jobs')
        delete(job{g});
    end
end

