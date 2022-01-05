function [shuffle, data_source] = get(cache_args_OR_cache_obj, iShuff, varargin)
% Acquires a shuffle from a cache object or from cache arguments


ip = inputParser;
ip.addParameter('cacheMethod', 'matfile'); % {matfile} | parquet | ''
ip.addParameter('debug', false); 
ip.parse(varargin{:})
Opt = ip.Results;


if istable(cache_args_OR_cache_obj)

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

elseif isa(cache_args_OR_cache_obj, 'matlab.io.MatFile')

    if Opt.debug
        disp("Pulling struct out of matfile and converting to table");
        keyboard
    end
    tmp = cache_args_OR_cache_obj.shuffle(iShuff,1);
    tmp = struct2table(tmp);
    shuffle = units.shuffle.get(tmp, iShuff, varargin{:});
    
elseif iscell(cache_args_OR_cache_obj)

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
