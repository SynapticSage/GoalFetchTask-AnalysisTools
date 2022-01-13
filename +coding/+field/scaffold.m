function scaff = scaffold(tab, varargin)
% Uses Opt.scaff instructions to provide

ip = inputParser;
ip.KeepUnmatched = true;
ip.addParameter('grid', 25);
ip.addParameter('scaff', []);
ip.addParameter('props', []);
ip.parse(varargin{:})
Opt = ip.Results;
Opt.props = string(Opt.props);


%  CALCULATE CENTERS
%  -----------------
scaff = struct('center', nan,...
    'edge', nan, 'col', nan, 'field', string(nan));

if ~isempty(Opt.scaff)
    scaff = util.struct.update(scaff, Opt.scaff);
end

if ~isempty(Opt.grid) 
    if ~isstruct(Opt.grid)
        scaff = util.struct.update(scaff, struct('grid', {Opt.grid}));
        useStruct = false;
    else
        useStruct = true;
    end
end

% ITERATE EACH FIELD
% ------------------
% column combination and acertain the centers and edges
scaffcount = 0;
for fcount = 1:numel(Opt.props)

    field = Opt.props(fcount);

    for col = 1:size(tab.(field),2)
        scaffcount = scaffcount + 1;
        scaff(scaffcount)      = scaff(1);
        scaff(scaffcount).field = field;
        scaff(scaffcount).col   = col;
        scaff(scaffcount).min  = min(tab.(field));
        scaff(scaffcount).max  = max(tab.(field));
        if useStruct
            scaff(scaffcount).grid = Opt.grid.(field);
        end
        scaff(scaffcount).edge = linspace(...
            scaff(scaffcount).min, ...
            scaff(scaffcount).max, ...
            scaff(scaffcount).grid+1);
        scaff(scaffcount).center = mean([scaff(scaffcount).edge(1:end-1);...
                                         scaff(scaffcount).edge(2:end)], 1);
    end
end

