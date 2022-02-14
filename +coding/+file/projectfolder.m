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

% And make sure the referenced thing exists, else create the directory path
[subdirectory,~,ext] = fileparts(directory);
if ~isempty(ext) % addressing a file
    ensurePathExist(subdirectory);
else  % addressing a directory
    ensurePathExist(directory);
end

function ensurePathExist(directory)

    if ~exist(directory, 'dir')
        [subdirectory, ~, ~] = fileparts(directory);
        ensurePathExist(subdirectory);
        disp("Making " + directory);
        mkdir(directory)
    else
        return
    end
    if isequal(directory, '/')
        return
    end
