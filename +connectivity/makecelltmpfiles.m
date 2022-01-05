function [spikesfile, outfile] = makecelltmpfiles(spikes, method)
% Creates the spikes files in the format expected by the glmcc.py
% file

if nargin == 1
    method = "connect";
end

firingFolder = '/tmp/cellFiring';
if ~exist(firingFolder,'dir')
    mkdir(firingFolder);
end
outfile = fullfile(firingFolder, 'estimated.csv');

currdir = pwd;
destruct = onCleanup(@() cd(currdir));
cd('/tmp/cellFiring');

if method == "glmcc"

    for iCell = progress(1:numel(spikes.spikeTimes),'Title','Making GLMCC temp files')
        tab = table(spikes.spikeTimes{iCell}');
        spikesfile{iCell} = sprintf('cell%d.txt', iCell));
        writetable(tab, fullfile(firingFolder, spikesfile{iCell},...
            'WriteVariableNames', false);
    end

elseif method == "connect"

    spikesfile = fullfile(firingFolder, 'spikesfile.txt');
    delete(spikesfile);
    for iCell = progress(1:numel(spikes.spikeTimes),'Title','Making CoNNECT temp files')
        tab = table(spikes.spikeTimes{iCell}');
        writetable(tab, spikesfile,...
            'WriteVariableNames', false, 'WriteMode', 'append');
        if iCell ~= numel(spikes.spikeTimes)
            fid = fopen(spikesfile, 'a');
            fwrite(fid, ';');
            fclose(fid);
        end
    end

end
