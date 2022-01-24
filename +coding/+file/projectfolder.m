function directory = projectfolder(varargin)

if numel(varargin) == 2 && strcmp(varargin{1}, 'set')
    Opt.set = varargin{2};
    varargin(1:2) = [];
else
    Opt.set = [];
end

if ~isempty(Opt.set)
    directory = string(Opt.set);
    setpref(mfilename(), 'directory', directory);
else
    directory = getpref(mfilename(), 'directory');
    directory = fullfile(directory, varargin{:});
end
