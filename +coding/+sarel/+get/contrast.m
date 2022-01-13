function contrast = contrast(vals, dim, varargin)

ip = inputParser;
ip.addParameter('nullIndex', 1);
ip.KeepUnmatched = true;
ip.parse(varargin{:});
Opt = ip.Results;

% Obtain null distributions
select = num2cell(repmat(':', ndims(vals)));

select_null = select;
select_null{dim} = Opt.nullIndex;

measurement_inds = setdiff(1:size(vals, dim), Opt.nullIndex);
select_measurement = select;
select_measurement{dim} = measurement_inds;

nulvals = vals(select_null{:});
measurements = vals(select_measurement{:});

% Compute index
contrast = (measurements - nulvals)./(measurements + nulvals);
