function [out, groups] = shuffle_within_condition(beh, spikes, groupby, varargin)
% Circularly shuffles within times per some conditioned group
%
%
% Uses:
% (1) Sarel 2017 and Dotson 2021 each requires circular spike time shuffles
% within trials. Doable by setting groupby to your relevent trial var name in
% the behavior table.
%
% -------
% Inputs
% -------
% 
% beh : table-type
%   This is a time length table, where columns are properties of the animal's
%   behavior
%
% spikes : struct
%   This is the spikes struct emitted by the units.getRate() method. It contains
%   all of the information that we will use about spiking
%
%
% groupby : list[string]
%   This is a list of propertiees to shuffle within. We shuffle within the
%   logical and of unique elements of these properties.
%
%
% ---------
% Optionals
% ---------
%
% Descriptions of the optional inputs can be found in the optional section 
% below.
%
% -------
% Outputs
% -------
%
% out : 
%
%   contains the ram or matfile based cache if those methods are seelceteed
%
% groups : 
%
%   the groups struct from util.findgroups(beh, props)
%
% -----
% Notes
% -----
%
%
% This is an extremely flexible method that can be used for overall shuffles
% or within property
%
% If you're using a shifts property in your shuffles, be aware that the set of
% tables can be very large. 100 shuffles of just x/y/time data for all neurons
% takes up 70GB on my machine. if you have those shuffles resampled for time-shifts
% (not shuffle-time shifts, but time-shifts as in the dotson and yartsev paper), 
% then you will take up 1400GB of space.
% 
% In order to scale to that level, there's an option for the atBehavior method
% that allows returning the indices of the beh table to sample rather
% than the actual behavior/time at those indices. This is a much more compressed
% way of doing things.

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
ip.addParameter('shuffleunits', 'unitwise'); % shuffle neurons so that {unitwise}|uniform
ip.addParameter('shiftstatistic', 'uniform'); % what statistic of shift? {uniform}|normal
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

ip.parse(varargin{:})
Opt = ip.Results;
Opt.shuffleunits   = string(Opt.shuffleunits);
Opt.shiftstatistic = char(Opt.shiftstatistic);
Opt.cacheMethod    = char(lower(Opt.cacheMethod));
Opt.kws_atBehavior = ip.Unmatched;

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
    groups = util.table.findgroups(beh, groupby);
else
    groups = Opt.groups;
end
if Opt.dropGroupby
    beh(:, groupby) = [];
end

checksum = all(ismember(unique(diff(groups.time.groups(groups.time.groups~=-1))), [0,1]));
assert(checksum, 'Non-contiguous regions in your conditionals!')

[measure] = measuretimeperiods(beh, groups);
[shifts] = calculateshifts(spikes, groups, measure, Opt);
[out, Opt] = prepareOuts(spikes, groupby, shifts, groups,  Opt);

% -------
% Shuffle
% -------
switch Opt.shiftWhat

    case 'behavior'
        out = behaviorbasedshuffle(out, shifts, groups, spikes, beh, Opt);

    case 'spikes'

end

if isfield(out, 'beh') && iscell(out.beh) && istable(out.beh{1})
    out.beh = util.table.icat(out.beh);
end

%---------------------------------------------------------------------- 
%--------------- HELPER FUNCTIONS ------------------------------------- 
%---------------------------------------------------------------------- 
function measure = measuretimeperiods(beh, groups)
% Obtains the start and stop periods

    nG = groups.nGroups;
    measure = table(nan(nG,1), nan(nG,1), nan(nG,1),...
        'VariableNames', ["start","stop","len"]);

    for g = groups.uGroups'
        measure.start(g) = beh.time(find(groups.time.groups == g, 1, 'first'));
        measure.stop(g)  = beh.time(find(groups.time.groups == g, 1, 'last'));
    end
    measure.len = measure.stop - measure.start;

function [shifts] = calculateshifts(spikes, groups, measure, Opt)

    % Calculate shifts
    if strcmp(Opt.shuffleunits, 'uniform')
        shifts = nan(Opt.nShuffle, 1, groups.nGroups, 'single');
    elseif strcmp(Opt.shuffleunits, 'unitwise')
        shifts = nan(Opt.nShuffle, height(spikes.cellTable), groups.nGroups, 'single');
    end

    [~, ~, G] = util.ndgrid.coord(shifts);
    if Opt.width == "whole"
        W = measure.len(G);
    else
        error("Not implemented")
    end

    switch Opt.shiftstatistic
        case 'uniform'
            shifts = W .* (rand(size(shifts)) - 0.5);
        case 'normal'
            shifts = W .* randn(size(shifts));
    end

function [out, Opt] = prepareOuts(spikes, groupby, shifts, groups, Opt)
% SET UP OUTPUTS
    
    % Setup outfolder
    % ---------------
    if isempty(Opt.outfolder)
        switch Opt.cacheMethod
        case 'parquet'
            if isempty(Opt.cacheToDisk)
                error('Must provide arguments to coding.file.parquetfolder as a cell array passed to Opt.cacheToDisk')
            end
            Opt.outfolder = coding.file.parquetfolder(Opt.cacheToDisk{:});
        end
    end
    if strcmp(Opt.cacheMethod, 'parquet') && ~Opt.skipShuffled
        disp("Deleting previous shuffles")
        % Remove any former shuffle data
        for file = dir(fullfile(Opt.outfolder), '*.parquet')'
            delete(fullfile(file.folder, file.name));
        end
    end
    % Setup struct
    % ---------------
    switch Opt.cacheMethod
    case 'matfile'
        out = coding.file.shufflematfile(Opt.cacheToDisk{:});
    case 'parquet'
        out = struct();
    otherwise
        out = struct();
    end
    out.Groups = groups;
    out.nShuffle = Opt.nShuffle;
    out.nNeurons = height(spikes.cellTable);
    out.nGroups  = groups.nGroups;
    out.groupby = groupby;
    out.shifts = shifts;
    out.Opt = Opt;

function out = behaviorbasedshuffle(out, shifts, groups, spikes, beh, Opt)

        % Ready input and output variables
        beh = util.table.castefficient(beh, 'negativeNan', true);
        nNeurons = height(spikes.cellTable);
        original_time = beh.time;
        Opt.preallocationSize = min(Opt.nShuffle, Opt.preallocationSize);
        if ismember('shuffle', fieldnames(out))
            util.matfile.rmsinglevar(out.Properties.Source, 'shuffle');
        end

        % Iterate through blocks of shufflse
        prevShuff = -1;
        singleShuffle = cell(1, nNeurons);
        for s = progress(1:Opt.preallocationSize:Opt.nShuffle,...
                'Title', 'Shuffling behavior')

            endPoint = min(s+Opt.preallocationSize-1, Opt.nShuffle);
            piece.shifts = shifts(s:endPoint,:,:);
            beh.time = original_time;

            %% THE MAGIC : where group-shuffle happens 
            %% -------------------------------------   
            piece.newtimes = ...
                units.shuffle.helper.preallocateBehaviorShifts(piece.shifts,...
                beh, groups);
            %% ------------------------------------- %%
            clear shuffle

            % QUERY shuffle times rom neural data
            % ------------------------------------
            SN = util.indicesMatrixForm([endPoint-s+1, nNeurons]);
            shuffle = cell(1, Opt.nShuffle);
            count = 0;
            pfExist = @(folder, file) exist(fullfile(folder, file), 'file')>0;
            for sn = progress(SN', 'Title', 'Shifting cells')
                count = count + 1;
                iShuff  = sn(1) + s - 1;
                iPreallocShuff  = sn(1);
                iNeuron = sn(2);

                % If we change to a new shuffle, store what we computed per cell
                notSkipProcessingPreviousShuff = ~Opt.skipShuffled || ~pfExist(Opt.outfolder, Opt.parquetfile(prevShuff));
                if (prevShuff ~= -1) && iShuff ~= prevShuff ...
                    && notSkipProcessingPreviousShuff
                    tmp = util.cell.icat(singleShuffle, 2);
                    assert(numel(tmp) == 1);
                    clear singleShuffle
                    singleShuffle = cell(1, nNeurons);
                    [out, shuffle] = cache(out, tmp{1}, shuffle, prevShuff, Opt);
                end

                skipProcessingCurrentShuffle = Opt.skipShuffled && pfExist(Opt.outfolder, Opt.parquetfile(iShuff));
                if skipProcessingCurrentShuffle
                    disp("skipping")
                    continue
                end

                beh.time = squeeze(piece.newtimes(iPreallocShuff, iNeuron, :));
                tmp = ...
                    units.atBehavior_singleCell(spikes.spikeTimes{iNeuron}, ...
                    beh, Opt.kws_atBehavior);
                if iscell(tmp) && istable(tmp{1})
                    tc = @util.type.castefficient;
                    for i = 1:numel(tmp)
                        tmp{i}.neuron  = tc(iNeuron * ones(height(tmp{i}), 1));
                        tmp{i}.shuffle = tc(iShuff * ones(height(tmp{i}), 1));
                        tmp{i}.shift   = tc(i * ones(height(tmp{i}), 1));
                    end
                    tmp = cat(1,tmp{:});
                end
                singleShuffle{iNeuron} = tmp;
                prevShuff = iShuff;
            end
        end % END of shuffle loop

        % Final expunge of shuffle contents
        if ~isempty(singleShuffle) %happens when processing last is skipped
            tmp = util.cell.icat(singleShuffle, 2);
            clear singleShuffle 
            [out, shuffle] = cache(out, tmp{1}, shuffle, prevShuff, Opt);
        end

        % ------------------------------------
        % Concatonate and label and cache
        % ------------------------------------
        if ~isempty(shuffle{1})
            shuffle = util.cell.icat(shuffle, 2);
            shuffle = shuffle{1};
        end
        out.beh = shuffle;

function [out, shuffle] = cache(out, singleShuffle, shuffle, prevShuff, Opt)
    switch Opt.cacheMethod
    case 'matfile'
        thing = table2struct(singleShuffle, 'ToScalar', true);
        out.beh(shuff, 1) = thing;
        shuffle{prevShuff} =  singleShuffle;
    case 'parquet'
        parquetwrite(fullfile(Opt.outfolder, Opt.parquetfile(prevShuff)), singleShuffle);
    case 'ram'
        shuffle{prevShuff} = singleShuffle;
    end
