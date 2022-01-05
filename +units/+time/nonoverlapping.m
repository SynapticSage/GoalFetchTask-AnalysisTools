function [t_midpoints, t_startends] = nonoverlapping(epochPeriods, samplingPeriod)

if isa(epochPeriods,'single')
    samplingPeriod = single(samplingPeriod);
end

t_startends = [];
for epoch = 1:size(epochPeriods, 1)
    start = epochPeriods(epoch, 1);
    stop =  epochPeriods(epoch, 2);
    t_startends = [t_startends, start:samplingPeriod:stop];
end

t_midpoints = [t_startends(1:end-1)',...
    t_startends(2:end)'];
t_midpoints = mean(t_midpoints,2)';

if all(t_midpoints==0)
    error('Fuck'); 
end

if ~util.isunique(t_midpoints)
    warning('Spike times not unique')
    keyboard;
end

