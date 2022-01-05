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

% ------------------------------------
% Load previous analyses if they exist
% ------------------------------------
if coding.file.exist(animal, day)
    spikes = coding.file.load(animal, day,...
        'appendTo', spikes);
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


% --------------------------------------------------------------
% ------ JERCOG/ABOTT/KANDEL PAPER -----------------------------
% --------------------------------------------------------------
if ismember("jercog", Opt.analyses)
    % Get sparse behavior at spike times for the next few analyses
    beh = gbehavior.lookup(animal, [], day);
    [spikes, beh, ~] = units.atBehavior(beh, spikes,...
        'useGPU', false,....
        'merge', true,...
        'query', Opt.behFilter);
    Opt.checkpointVar = 'jercogtmp';
    spikes = jercog(spikes, beh, Opt);
end%jercog


% -----------------------------------------------------------
% ------ SAREL/OKEEFE PAPER ---------------------------------
% -----------------------------------------------------------
% Tuning to angle and distance?
% Jercog Kandel and Abott Paper
% This method simultaneously computes over data splits
% and doesn't require for-looping
if ismember("sarel", Opt.analyses)

    % Perform sarel analyses
    beh = gbehavior.lookup(animal, [], day);
    profile on;
    [spikes, beh, ~] = units.atBehavior(beh, spikes,...
                                        'merge', true,...
                                        'query', Opt.behFilter);
    profile viewer;
    spikes.sarel = coding.sarel.sarel(spikes, beh);

    % Now let's get some shuffles to test some things about what we measured
    spikes = units.shuffle.conditional_time(beh, spikes, ...
        'within_condition', "trial",... % shuffle within individual trials
        'atBehavior_kws', {'query', Opt.behFilter})
    spikes.sarel.ground.goalplaceindex = ...
        coding.sarel.computerGoalPlaceIndex(spikes, beh, 'ground');
    spikes.sarel.stops.goalplaceindex = ...
        coding.sarel.computerGoalPlaceIndex(spikes, beh, 'stops');

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

% -------------------------------------------------------------------
% ------ DOTSON/YARTSEV nonlocal  PAPER -----------------------------
% -------------------------------------------------------------------
%  TODO LIST
%  ---------
% 1) Shifter function: shift spike timse and recompute beh table ✅
% 2) GPU-optomized place field compute ✅
% -) Test thus far 
% 3) Field-post-processors
%   - Entropy per shift
%   - Convex hulls of best and time 0
% 4) Connectivity of ca1 - pfc
%   Measures
%       - Cross-correlation
%       - GLM (how well pfc ensemble predicts a given ca1 cell)
% 5) Relate future/past/present of ca1 cells to (4)
%   - Cells with higher connectivity more future?
%   - PFC more future?
%   - Reference points along goal-directed paths also more future fields?

if ismember("dotson", Opt.analyses)

    futurepastKws = struct('useGPU', true, 'grid', 15);

    % -----------------------------------------------
    %     Acquire slice of behavior we'll cook with
    % -----------------------------------------------
    beh = gbehavior.lookup(animal, [], day);
    beh.x = beh.pos(:,1);
    beh.y = beh.pos(:,2);
    props = ["x", "y"];
    beh = util.table.castefficient(beh,...
        'compressReals', true,...
        'negativeNan', true,...
        'compressRealsExcept', {'time'});
    shift = -1:0.1:1;

    % -------------------------
    % Hard part: the shuffle
    % -------------------------
    
    % Pregenerate our shuffled indices we grab per neuron from beh
    % ------------------------------------------------------------
    shuffKws = struct('groups', [], 'props', props, 'shift', shift,...
        'cacheToDisk', {{animal, day}},...
        'skipShuffled',  1, ...
        'returnIndices', 1,...
        'startShuffle', 1,...
        'nShuffle', 50,...
        'preallocationSize', 50);
    groupby = ["epoch", "period"]; % properties in behavior table to shuffle around
    units.shuffle.conditional_time(beh, spikes, groupby, shuffKws);

    %    Get shifted fields for shuffles      
    % ----------------------------------------
    for iS = progress(30:Opt.nShuffle, 'Title', 'Shuffles')
        % get
        Shuf = {animal, day};
        [item, Shuf] = units.shuffle.get(Shuf, iS, 'debug', false);
        item = struct('beh', item);
        item.shift = shift;
        item.behtype = 'indices';
        % run
        spikes.dotson.shuffle{iS} = coding.futurepast.main(item, beh, props, ...
            futurepastKws); 
        % takes about 24m with GPU for 15^2 grid
        % means 100x for shuffles, so 41 hours with GPU
        % (136 days without GPU support)
    end

    % --------------------------------------------------------
    % Easy part: measure the real data and debias with shuffle
    % --------------------------------------------------------

    % Meausre the real data
    [spikes, beh, ~] = units.atBehavior(beh, spikes,...
        'useGPU', false,....
        'merge',  false,...
        'returnIndices', true, ...
        'shift', shift,...
        'query', Opt.behFilter);
    spikes.doston = coding.futurepast.main(spikes, beh, props, futurepastKws);

    % We debias the result by substracting the shuffles
    spikes.doston = coding.futurepast.shuffcorrect(spikes.dotson);

end %dotson

coding.file.save(animal, day, 'append', true);
