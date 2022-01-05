function m = shufflematfile(animal, index, varargin)

filename_full = coding.file.shufflematfilename(animal, index, varargin{:});
m = matfile(filename_full, varargin{:});

