function filename_full = shufflematfilename(animal, index, varargin)
% function filename_full = filename(animal, index)
ip = inputParser;
ip.addParameter('shifttype', []);
ip.addParameter('groupby', []);
ip.addParameter('kws_atBehavior', []);
ip.addParameter('query', []);
ip.KeepUnmatched = true;
ip.parse(varargin{:})
Opt = ip.Results;

if ~isempty(Opt.query)
    filtportion = "_filt=" + ...
        replace(replace(strip(Opt.query), '$', ''), ' ', '');
else
    filtportion = "";
end

filename_full = coding.file.filename(animal, index, varargin{:});
replace_string = sprintf('_shuffle_groupby=%s_shifttype=%s%s.mat', ...
    join(Opt.groupby, "-"),...
    Opt.shifttype,...
    filtportion);
filename_full = replace(filename_full, '.mat', replace_string);
