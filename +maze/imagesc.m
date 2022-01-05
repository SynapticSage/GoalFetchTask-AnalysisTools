function imagesc(animal, day, varargin{:});
% Maze specific version of imagesc

ip = inputParser;
ip.KeepUnmatched = true;
ip.addParameter('barriers', false);
ip.parse(varargin{:})
Opt = ip.Results;

unmatched = ip.Unmatched;


