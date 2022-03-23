% ------------------
% Get the path setup
% ------------------
addpath(genpath('~/Code/analysis'))
[animal, day, varargin] = deal('RY16', 36, {});
Opt = OPTARGS('unit', 'spikes');
Opt.shuffleStruct;
%Opt.shuffleStruct = []; % this makes program skip caching new shuffles, while Opt.nShuffle still acknowledges an ability to read out 50 of those cached items
Opt.analyses = [];
disp(newline)
disp("-------")
disp("Options")
disp("-------")
disp(Opt)
disp(newline)
disp(newline)
disp("-------")
disp("ShuffleOpt")
disp("-------")
disp(Opt.shuffleStruct)

% TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO 
%%  HARD CODE  : remove this
% TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO 

% -------------
% Get Unit data
% -------------
spikeDat = ndb.load(animal, Opt.unit, 'ind', day);
spikes   = units.getRateMatrix(animal, day,...
    'unit', Opt.unit,...
    'taskFilter', Opt.taskFilter,...
    'ploton',false,...
    'dense', false,...
    'cellFilter', Opt.cellFilter,...
    'spikeData', spikeDat);
clear spikeDat

% ------------------------------------
% Load previous analyses if they exist
% ------------------------------------
if coding.file.exist(animal, day)
    spikes = coding.file.load(animal, day, 'appendTo', spikes);
end

% ------------------------------------------------------------
% Pregenerate our shuffled indices we grab per neuron from beh
% (these shuffles are used for all analyses below)
% ------------------------------------------------------------
%if ~isempty(Opt.shuffleStruct)
%    disp("Shuffling")
%    tmp = struct('groups', [],...
%        'cacheToDisk', {{animal, day}},...
%        'shift', Opt.shift,...
%        'returnIndices', 1);
%    Opt.shuffleStruct = util.struct.update(Opt.shuffleStruct, tmp); clear tmp;
%    groupby = ["epoch", "period"]; % properties in behavior table to shuffle around
%beh = gbehavior.lookup(animal, [], day);
%    units.shuffle.shuffleWithinConditions(beh, spikes, groupby, Opt.shuffleStruct);
%else
%    disp("Skip shuffling")
%end

tmp = struct('groups', [],...
'cacheToDisk', {{animal, day}},...
'shift', Opt.shift,...
'query', Opt.behFilter,...
'returnIndices', 1);
Opt.shuffleStruct = util.struct.update(Opt.shuffleStruct, tmp); clear tmp;
