function [shuffle, data_source] = get(cache_args_OR_cache_obj, iShuff, varargin)
% Acquires a shuffle from a cache object or from cache arguments


ip = inputParser;
ip.addParameter('cacheMethod', 'matfile'); % {matfile} | parquet | ''
ip.addParameter('debug', false); 
ip.addParameter('shiftless', false); % when true, if the shuffles have a shift property, this selects the shift=0seconds portion (unshifted)
ip.addParameter('shift', []); % when not empty, this shift index is looked up, instead of returning all
ip.parse(varargin{:})
Opt = ip.Results;


if istable(cache_args_OR_cache_obj)
    % -----
    % TABLE
    % -----

    if Opt.debug
        disp("Indexing shuffle out of table");
        keyboard
    end
    if any(strcmp('shuffle', cache_args_OR_cache_obj))
        shuffle = cache_args_OR_cache_obj(cache_args_OR_cache_obj.shuffle == iShuff, :);
    else
        shuffle = cache_args_OR_cache_obj;
    end
    data_source = [];

    % lookup a specific shift
    if ~isempty(Opt.shift)
        shuffle = shuffle(shuffle.shift == Opt.shift, :);
    end

elseif isa(cache_args_OR_cache_obj, 'matlab.io.MatFile')
    % -------
    % MATFILE
    % -------

    if Opt.debug
        disp("Pulling struct out of matfile and converting to table");
        keyboard
    end
    tmp = cache_args_OR_cache_obj.shuffle(iShuff,1);
    tmp = struct2table(tmp);

    if Opt.shiftless
        if ~ismember(fieldnames(cache_args_OR_cache_obj), 'zeroshiftIndex')
            zeroshiftIndex = median(unique(tmp.shift));
            cache_args_OR_cache_obj.Properties.Writable = true;
            cache_args_OR_cache_obj.zeroshiftIndex = median(unique(tmp.shift));
            cache_args_OR_cache_obj.Properties.Writable = false;
        else
            zeroshiftIndex = cache_args_OR_cache_obj.zeroshiftIndex;
        end
        V = util.struct.varargin2struct(varargin);
        V.shift = zeroshiftIndex;
        shuffle = units.shuffle.get(tmp, iShuff, V);
    else
        shuffle = units.shuffle.get(tmp, iShuff, varargin{:});
    end

    % Format the result into an output struct that partially matches
    % the spikes struct used in all my other functions
    shuffle = struct('beh', shuffle);
    if ismember('indices', shuffle.beh.Properties.VariableNames)
        shuffle.behtype = 'indices';
        if ismember(fieldnames(cache_args_OR_cache_obj), 'uShift')
            shuffle.uShift = cache_args_OR_cache_obj.uShift;
        else
            uShift = unique(shuffle.beh.shift);
            cache_args_OR_cache_obj.Properties.Writable = true;
            cache_args_OR_cache_obj.uShift = uShift;
            cache_args_OR_cache_obj.Properties.Writable = false;
            shuffle.uShift = uShift;
        end
        if ismember(fieldnames(cache_args_OR_cache_obj), 'uNeuron')
            shuffle.uNeuron = cache_args_OR_cache_obj.uNeuron;
        else
            uNeuron = unique(shuffle.beh.neuron);
            cache_args_OR_cache_obj.Properties.Writable = true;
            cache_args_OR_cache_obj.uNeuron = uNeuron;
            cache_args_OR_cache_obj.Properties.Writable = false;
            shuffle.uNeuron = uNeuron;
        end
    end
    
elseif iscell(cache_args_OR_cache_obj)
    % --------------
    % FILE ARGUMENTS
    % --------------

    we_have_cache_args = @(x) isstring(x) || ischar(x);
    if we_have_cache_args(cache_args_OR_cache_obj{1})
        switch lower(Opt.cacheMethod)
            case 'matfile'
                if Opt.debug
                    disp("Opening matfile");
                    keyboard
                end
                data_source = coding.file.shufflematfile(cache_args_OR_cache_obj{:});
                shuffle = units.shuffle.get(data_source, iShuff, varargin{:});
            case 'parquet'
                folder = coding.file.shuffleparquetfolder(cache_args_OR_cache_obj{:});
                file = fullfile(folder, iShuff + ".parquet");
                shuffle = parquetread(file);
            otherwise
                error("Not implemented")
        end
    else
        error("Please provide cache arguments... see coding.file.* method of interest")
    end
end
