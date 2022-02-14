function tab = labeledSpiking(spikes, props, beh)
% Created a labeled tidy data structure for raw spiking events

tab = cell(numel(spikes.spikeTimes), 1);
for iCell = 1:numel(spikes.spikeTimes)
    time = spikes.spikeTimes{iCell};
    unit = repmat(iCell, numel(time), 1);
    spiketable = repmat(spikes.cellTable(iCell, :), numel(time), 1);
    spiketable.unit = unit(:);
    spiketable.time = time(:);
    % Add properties from behavior table?
    if ~isempty(props)
        for prop  = string(props(:))'
            inds = spikes.beh{iCell}.indices;
            nonzero = inds>0;
            x = nan(height(spiketable), 1);
            x(nonzero) = beh.(prop)(inds(nonzero)); 
            spiketable.(prop) = x;
        end
        spiketable.index = spikes.beh{iCell}.indices;
    end
    tab{iCell} = spiketable;
end

tab = util.cell.icat(tab);
tab = util.table.castefficient(tab);
