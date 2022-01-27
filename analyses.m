function [spikes, Opt] = analyses(animal, day, varargin)
% function stats = coordcentric(animal, day, varargin)
%
% This function organizes all of the goal- and future- coding measurements
% doing for R Young's paper.
%
% Until this file works cradle to grave, I'm running it as if it were a script.

Opt = OPTARGS(varargin{:});

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
if ~isempty(Opt.shuffleStruct)
    tmp = struct('groups', [],...
        'cacheToDisk', {{animal, day}},...
        'shift', Opt.shift,...
        'returnIndices', 1);
    Opt.shuffleStruct = util.struct.update(Opt.shuffleStruct, tmp); clear tmp;
    groupby = ["epoch", "period"]; % properties in behavior table to shuffle around
    beh = gbehavior.lookup(animal, [], day);
    % Any unmatched args are passed to atBehavior_singleCell
    units.shuffle.shuffleWithinConditions(beh, spikes, groupby, Opt.shuffleStruct);
end

% ========================================
% ,---.          |                        
% |---|,---.,---.|    ,   .,---.,---.,---.
% |   ||   |,---||    |   |`---.|---'`---.
% `   '`   '`---^`---'`---|`---'`---'`---'
%                     `---'               
% ========================================
if ismember("dotson", Opt.analyses)
    spikes = coding.futurepast.analysis(animal, day, spikes, Opt);
end %dotson


% --------------------------------------------------------------
% ------ JERCOG/ABOTT/KANDEL PAPER -----------------------------
% --------------------------------------------------------------
if ismember("jercog", Opt.analyses)
    spikes = coding.jercog.analysis(animal, day, spikes, Opt);
end%jercog


% -----------------------------------------------------------
% ------ SAREL/OKEEFE PAPER ---------------------------------
% -----------------------------------------------------------
if ismember("sarel", Opt.analyses)
    spikes = coding.sarel.analysis(animal, day, spikes, Opt);
end %sarel

% --------------------------------------------------------
% ------ IZTHAK/FRIED paper ------------------------------
% --------------------------------------------------------
%  TODO LIST
%  ---------
%  0) Method that identifies periods of heavy noise via
%  times of extreme broadband coherence 🚧
%   - create matfile for these in preprocess 🚧
%  1) Method that drops top 3 examples of bands per day
%  into lfp metadata table files, with noise appended
%   - create matfile for these in preprocess
%  2) Precession scales and bands:
%   a) theta + ripple
%   b) scales
%       - start-end visit
%       - out-and-back to home
%       - block of trials
if ismember("izthak", Opt.analyses)

    lfp_metadata = ndb.load(animal, 'lfpmeta'); % TODO
    lfp_metadata = lfp_metadata(...
        lfp_metadata.day == day,:);

    % PRECESSION
    for band = ["theta", "ripple"]
        best = lfp_metadata.("best_" + band);
        lfp = ndb.load(animal, band, 'inds', best);
        spikes = coding.precesssion.izthak(spikes, beh, lfp);
    end

end%fried

if ~isempty(Opt.analyses)
    coding.file.save(animal, day, spikes,...
        'append', true, 'checkpointActive', Opt.checkpoint);
end
