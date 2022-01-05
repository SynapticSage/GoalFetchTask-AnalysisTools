function m = matfile(animal, index, varargin)
filename_full = coding.file.filename(animal, index, varargin{:});
m = matfile(filename_full, varargin{:});
