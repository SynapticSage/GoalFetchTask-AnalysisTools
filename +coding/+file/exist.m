function answer = exist(animal, index, varargin)
% function answer = exist(animal, index)

ip = inputParser;
ip.KeepUnmatched = true;
ip.parse(varargin{:})
Opt = ip.Results;

filename_full = coding.file.filename(animal, index, varargin{:});
answer        = exist(filename_full, 'file') > 0;
