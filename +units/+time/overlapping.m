function [t_midpoints, t_startends] = overlapping(epochPeriods, samplingPeriod, window)

% Acquire all epoch starts stops
t_midpoints = arrayfun(...
    @(i) epochPeriods(i,1):samplingPeriod:epochPeriods(i,2),...
    1:size(spikes.timePeriods, 1),...
    'UniformOutput', false);
t_midpoints = cat(2, t_midpoints{:});

if isscalar(window)
    window = [-window(1)/2, window(2)/2];
end
assert(window(2) > window(1), 'End of window must be after start');

% get window edges around those time points (different knob for bin size here than samplingRate)
windowStarts = t_midpoints + window(1)/2;
windowStops =  t_midpoints + window(2)/2;
t_startends = [windowStarts; windowStops];
