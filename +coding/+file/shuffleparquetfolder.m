function folder = shuffleparquetfolder(animal, index, varargin)
% Returns folder to use for a set of parquet shuffles

[~, Props] = coding.file.filename(animal, index, varargin{:});
name = Props.headerdesc;
name = name + "_parquet";
folder = fullfile(Props.folder, name);
if ~exist(folder, 'dir')
    mkdir(folder);
end
