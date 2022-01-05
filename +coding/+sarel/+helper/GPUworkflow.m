function out = GPUworkflow(spikes, fields, Binning, Opt)
    % This workflow repeats, but marginalized in different ways
    
    preserve_fields = unique([fields, Opt.anglefield, Opt.distancefield, "stopWell", "startWell", "cuemem"]);
    spikes.beh = util.table.select(spikes.beh, preserve_fields);

    % Send all of the variables that would be usede to the
    spikes.beh = util.table.table2GPUtable(spikes.beh);
    Binning    = util.struct.struct2GPUstruct(Binning);

    out = coding.sarel.helper.workflow(spikes, fields, Binning, Opt);

    % Bring results back from the GPU
    out = nd.apply(out, "**", @gather, 'recurseCell', true);
    spikes.beh  = util.table.GPUtable2table(spikes.beh);
    reset(gpuDevice());
