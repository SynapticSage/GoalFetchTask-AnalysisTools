function [spikes, beh] = spikesbeh(spikes, beh, filtstr)
% Filter the spikes at behavior struct with a filter
% string

spikes = tabfilt.spikes(spikes, filtstr);
beh    = tabfilt.behavior(beh, filtstr);
