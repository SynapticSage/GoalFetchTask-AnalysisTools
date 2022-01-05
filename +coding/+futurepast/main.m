function D = main(spikes, beh, props, varargin)
% coding.futurepast.main
%
% Calculates optimally shifted fields per neurons
% based on behavioral properties props
%
% -----
% Input
% -----
%
% spikes : struct
%     expects a normal spikes struct with .beh field element. This field elememnt should be a cell of tables.
%
% ------
% Output
% ------
%
% D : struct
% contains fields regarding the futuriness/pastiness of each cells fields
% for the chosen properties of interest

ip = inputParser;
ip.addParameter('grid', 15);
ip.addParameter('useGPU', false);
ip.parse(varargin{:})
Opt = ip.Results;

if ~Opt.useGPU
    warning('Speed without GPU is 80x slower');
end

% --------------------------------------
% -- Get the shifted fields ------------
% --------------------------------------

if iscell(spikes.beh)
    shufset = 1:numel(spikes.shift);
else
    shufset = unique(spikes.beh.shift);
end

        
tic
for i = progress(1:numel(shufset(:)), 'Title', 'Computing each shift')
    
    % What type of input were we passed?
    % ----------------------------------
    if iscell(spikes.beh)
        val = spikes.beh{i};
    else
        good_shift_indices = spikes.beh.shift==shufset(i);
    end

    % Grab indices
    % ------------
    if any(strcmp(fieldnames(spikes), 'behtype')) &&...
            spikes.behtype == "indices"
        % Type  = index

        % Pull data from our shuffle index structure
        good_shift_indices = good_shift_indices & spikes.beh.indices > 0;
        indices_into_behavior_data = spikes.beh.indices(good_shift_indices);
        neuron_labels = spikes.beh.neuron(good_shift_indices);

        % Use those to index out behavior at spike times and label those times
        % for corresponding neurons
        val = beh(indices_into_behavior_data, :);
        val.neuron = neuron_labels;
    else
        % Type  = actual data, not an index
        val = val(good_shift_indices);
    end

    [fields(i), behField] = ...
        coding.field.calc(val, ...
        'props', props, ...
        'beh', beh, ...
        'useGPU', Opt.useGPU, ... % about 80x faster with GPU
        'grid', Opt.grid);

    assert(util.num.fractionQueryMatch(fields(i).FR_occNorm, @(x) ~isnan(x) & ~isinf(x)) > 0,...
        'assertion failed: something is up');
end
toc

% Combine the fields into a single object
D = coding.field.combine(fields);
D.behaviorfield = behField;
D.shift = spikes.shift;
D.MI = coding.futurepast.fieldMI(D.FR_occNorm, ...
                                 D.behaviorfield.visits,...
                                 'shiftdim', 1,... % dimension of time shifts       
                                 'celldim', 2,...
                                 'nBins', 10);
