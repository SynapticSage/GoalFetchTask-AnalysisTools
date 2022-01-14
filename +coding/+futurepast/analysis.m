function spikes = analysis(animal, day, spikes, Opt)

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


futurepastKws = struct('useGPU', true, 'grid', 15);
futurepastKws

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
    'skipShuffled',  0, ...
    'returnIndices', 1,...
    'startShuffle', 1,...
    'nShuffle', Opt.nShuffle,...
    'preallocationSize', Opt.nShuffle);
groupby = ["epoch", "period"]; % properties in behavior table to shuffle around
units.shuffle.conditional_time(beh, spikes, groupby, shuffKws);
util.notify.pushover('FuturePast','Finished conditional shuffle creation');

%    Get shifted fields for shuffles      
% ----------------------------------------
if ~isfield(spikes, 'fp')
    spikes.fp = struct();
end
if isfield(spikes.fp, 'shuffle') && ~iscell(spikes.fp.shuffle)
    spikes.fp = rmfield(spikes.fp, 'shuffle');
end
if ~isfield(spikes.fp, 'shuffle')
    spikes.fp.shuffle = {}; 
end

for iS = progress(1:Opt.nShuffle, 'Title', 'Shuffles')
    % get
    Shuf = {animal, day};
    [item, Shuf] = units.shuffle.get(Shuf, iS, 'debug', false);
    % run
    spikes.fp.shuffle{iS} = coding.futurepast.main(item, beh, props, ...
        futurepastKws); 
end
spikes.fp.shuffle = cat(1, spikes.fp.shuffle{:});
util.notify.pushover('FuturePast','Finished computing shifted shuffle effects')

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
disp("Starting main FP calculations")
spikes.fp.main = coding.futurepast.main(spikes, beh, props, futurepastKws);
util.notify.pushover('FuturePast','Computed main data effects');

% ---------------------------------------------
% De-biasing the main effect using shuftle data
% ---------------------------------------------
spikes.fp.main = coding.futurepast.shuffcorrect(spikes.fp);
