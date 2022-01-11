function [out, groups] = conditional_time(beh, spikes, groupby, varargin)
% Circularly shuffles within times per some conditioned group
%
%
% Uses:
% (1) Sarel 2017 and Dotson 2021 each requires circular spike time shuffles
% within trials. Doable by setting groupby to your relevent trial var name in
% the behavior table.
%
%
% -----
% Notes
% -----

%                                              
%                     ,---.     |    o          
%                    |   |,---.|--- .,---.,---.
%                    |   ||   ||    ||   ||   |
%                    `---'|---'`---'``---'`   '
%                         |                    
ip = inputParser;

ip.KeepUnmatched = true; % Any unmatched go to the called method units.atBehavior.m
% Usual params for that: 

% HOW do we shift?
ip.addParameter('nShuffle', 100, @isnumeric);
ip.addParameter('startShuffle', 1, @isnumeric); % the first shuffle index
ip.addParameter('endShuffle', [], @isnumeric); % the first shuffle index
ip.addParameter('shuffleunits', 'unitwise'); % shuffle neurons so that {unitwise}|uniform
ip.addParameter('shiftstatistic', 'uniform'); % what statistic of shift? {uniform}|normal
ip.addParameter('shifttype', 'circshift');
ip.addParameter('width', 'whole'); % draw the 'whole' period of time, or some specified amount or standard deviation
ip.addParameter('shiftWhat', 'behavior'); % It's equibalent to shift behavior times repeatedly per cell or spike times per cell, but my estimate is that it's less memory intense for behavior

% SAVE space?
ip.addParameter('throwOutNonGroup', false);
ip.addParameter('preallocationSize', 2); % number of shuffles to run at a time
ip.addParameter('props',[]);
ip.addParameter('dropGroupby', true);

% CACHE specific
ip.addParameter('skipShuffled', false); % requires a cacheToDisk scenario: if true, if it detects an existing shuffle at that index on the disk, it skips
ip.addParameter('cacheMethod', 'matfile'); % {matfile} | parquet | RAM, Speed: RAM>matfile>parquet, RAMspace: matfile>parquet>RAM
ip.addParameter('parquetfile', @(shuff) shuff + ".parquet"); % lambda defining parquet file name
ip.addParameter('cacheToDisk', {}); % used to define folder or cache file; if out to disk, this takes the parameters for coding.file.shufflefilename or coding.file.parquetfoldername
ip.addParameter('groups', []); % pass in computed groups rather than labels from which we have to compute them? for saving computational time,
ip.addParameter('outfolder', []); % pass this to set the outfolder, instead of deriving it from cache-specific methods called on cacheToDisk parameters
%ip.addParameter('lazy', false); % instead save the parameters needed to yield shuffles on the fly

% UNMATCHED PARAMS --> go to units.atBehavior_singleCell.m

ip.parse(varargin{:})
Opt = ip.Results;
Opt.shuffleunits   = string(Opt.shuffleunits);
Opt.shiftstatistic = char(Opt.shiftstatistic);
Opt.cacheMethod = char(lower(Opt.cacheMethod));

Opt.kws_atBehavior = ip.Unmatched;
Opt.kws_atBehavior.maxNeuron = numel(spikes.spikeTimes);
if isempty(Opt.endShuffle)
    Opt.endShuffle = Opt.nShuffle;
end

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
[out, Opt] = units.shuffle.helper.prepareOuts(spikes, groupby, shifts, groups,  Opt);

% -------
% Shuffle
% -------
switch Opt.shiftWhat
    case 'behavior'
        out = units.shuffle.helper.behaviorbasedshuffle(out, shifts, groups, spikes, beh, Opt);
    case 'spikes'
end

if isfield(out, 'beh') && iscell(out.beh) && istable(out.beh{1})
    out.beh = util.table.icat(out.beh);
end
