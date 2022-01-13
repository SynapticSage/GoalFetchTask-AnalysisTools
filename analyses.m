function spikes = analyses(animal, day, varargin)
% function stats = coordcentric(animal, day, varargin)
%
% This function organizes all of the goal- and future- coding measurements
% doing for R Young's paper.
%
% Until this file works cradle to grave, I'm running it as if it were a script.

[animal, day, varargin] = deal('RY16', 36, {});

ip = inputParser;
ip.addParameter('unit', 'multiunit');                       % Spikes (curated firings) / Multiunit (uncurated firings)

ip.addParameter('shuffleStruct', []);                       % Struct of shuffling instructions (Reruns all analyses requested with sequences of shuffles specified here)
ip.addParameter('nullMethod', []);
ip.addParameter('nShuffle', 50);                            % Number of shuffles to generate
ip.addParameter('upsample', false);                         % Upsample behavior?

% Filtration
ip.addParameter('behFilter', 'abs($velVec) > 4');           % Query to apply to all behavior
ip.addParameter('taskFilter', "string($type)==""run""");    % 
ip.addParameter('cellFilter', 'ismember(string($area), ["CA1","PFC"]) & string($tag) ~= "artifcat" & $numspikes > 150'); % 

% Analyses
ip.addParameter('analyses',   ["jercog", "sarel", "dotson"]);
ip.parse(varargin{:})
Opt = ip.Results;

%shuffleStructDefault

% ------------------
% Get the path setup
% ------------------
addpath(genpath('~/Code/analysis'))

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

%% -------------------
%% Null Distributions? 
% This will create a null distribution and rerun all major analyses with it
%% -------------------
if ~isempty(Opt.nullMethod)
    for method = string(Opt.nullMethod)'
        null.(method) = coding.nullFiring(spikes.data, 'method', method);
    end
    disp(null)
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

coding.file.save(animal, day, spikes, 'append', true);
