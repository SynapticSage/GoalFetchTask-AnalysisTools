function filename_full = shufflematfilename(animal, index, varargin)
% function filename_full = filename(animal, index)
 
filename_full = coding.file.filename(animal, index, varargin{:});
filename_full = replace(filename_full, '.mat', '_shuffle.mat');
