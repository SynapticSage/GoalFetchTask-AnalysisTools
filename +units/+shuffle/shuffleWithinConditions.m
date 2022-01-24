function [out, groups] = shuffleWithinConditions(beh, spikes, groupby, varargin)
% Circularly shuffles within times per some conditioned group
%
%
% General shuffle method that caches the shuffles for later use....
%
%
% Used in...
% (1) Sarel 2017 and Dotson 2021 each requires circular spike time shuffles
% within trials. Doable by setting groupby to your relevent trial var name in
% the behavior table.
%
% -----
% Notes
% -----

% --------------------------------------------------------------
%                                              
%                     ,---.     |    o          
%                    |   |,---.|--- .,---.,---.
%                    |   ||   ||    ||   ||   |
%                    `---'|---'`---'``---'`   '
%                         |                    
% --------------------------------------------------------------
Opt = units.shuffle.optargs(varargin{:});
Opt.kws_atBehavior.maxNeuron = numel(spikes.spikeTimes);
Opt.shufflename = sprintf('shuffleWithinConditions=[%s]', ...
    join(groupby, "-"));
Opt.shufflename = string(Opt.shufflename);
    

if ~isempty(Opt.props)
    Opt.props = union("time", string(Opt.props));
    Opt.props = union(groupby, string(Opt.props));
    beh = beh(:, Opt.props);
end

                                    
%----------------------------------------
%                     ,-.-.     o     
%                     | | |,---..,---.
%                     | | |,---|||   |
%                     ` ' '`---^``   '
%----------------------------------------
% Create set of shuffle times
%----------------------------------------
% S x G x 1 if uniform or S x G x N if unit-based
if isempty(Opt.groups)
    disp("Finding groups")
    groups = util.table.findgroups(beh, groupby);
else
    groups = Opt.groups;
end
if Opt.dropGroupby
    beh(:, groupby) = [];
end

checksum = all(ismember(unique(diff(groups.time.groups(groups.time.groups~=-1))), [0,1]));
assert(checksum, 'Non-contiguous regions in your conditionals!')

measure = units.shuffle.helper.measuretimeperiods(beh, groups);
shifts = units.shuffle.helper.calculateshifts(spikes, groups, measure, Opt);

% Prepare the output structure and document our settings
[out, Opt] = units.shuffle.helper.prepareOuts(spikes, groupby, shifts, groups,  Opt);

% -------
% Shuffle
% -------
switch Opt.shiftWhat
    case 'behavior'
        % if matfile cache, then input(out) gets filed under out.shuffle,
        % i.e., out.shuffle = out, on the output end
        %
        % potentially confusing, I must admit
        out = units.shuffle.helper.behaviorbasedshuffle(out, shifts, groups, spikes, beh, Opt);
    case 'spikes'
end

if isfield(out, 'beh') && iscell(out.beh) && istable(out.beh{1})
    out.beh = util.table.icat(out.beh);
end
