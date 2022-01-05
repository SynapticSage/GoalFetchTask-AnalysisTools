function place(ratetime, rate, postime, position, varargin)
% Calculates a field over properties in position

ip = inputParser;
ip.addParameter('grid', 100);
ip.parse(varargin{:})
Opt = ip.Results;

bounds = [min(position); max(position)];
edges   = arrayfun(@(x) linedge(bounds(x,1),bounds(x,2), Opt.grid+1), 1:size(bounds<2));
centers = arrayfun(@(x) linedge(bounds(x,1),bounds(x,2), Opt.grid+1), 1:size(bounds<2));

position = num2cell(position, 2);

% OCCUPANCY
clear bin
if numel(edges) == 2 % 2D
    [N,~,~,bin{1},bin{2}] = histcounts2(position{:}, edges{:});
elseif numel(edges) == 1 % 1D
    [N,~,bin{1}] = histcounts(position{:}, edges{:});
else % ND
    gridOfSpace = cell(size(edges));
    [gridOfSpace{:}] = ndgrid(edges{:});
    error("ND code not completed yet");
end

% RATE PER POSITION
pos2rate = interp1(ratetime, 1:numel(ratetime), postime, 'nearest');
binrate  = cellfun(@(x) x(rate2pos,:), bin);

% Occupancy normalize the rate


function y = distance(x, varargin)
% Distance metric used to match elements

n = numel(x);
metric = arrayfun(@(i) shift(-1, (x(i) - varargin{i}).^2), 1:n, 'UniformOutput', false);
y = cat(1, metric{:});
y = sqrt(sum(y, 1));



