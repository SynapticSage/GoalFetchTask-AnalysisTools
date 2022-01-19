function comparison = compare(sarel, varargin)
% compares shuffle to main data in the sarel struct
%
% Input
% Struct details
%   | Levels     | Purpose                                     |
%   |------------+---------------------------------------------+
%   | splits     | ways of splitting time length               |
%   |            | to measure tuning curves                    |
%   |            |                                             |
%   | datafields | the main tuning curve types                 |
%   |            |                                             |
%   | metrics    | *measurements* done on tuning curves        |
%   |            |                                             |
%   | stats      | further measures done either on metrics     |
%   |            | or main/shuffle of tuning curves or metrics |

ip = coding.sarel.helper.mainStructInputParser();
ip.addParameter('method', @minus);
ip.parse(varargin{:})
Opt = ip.Results;

first = struct('type', '()', 'subs', {{1}});
tuning_curveFields = @(x) util.struct.matchingfields(x,...
    'logic', 'or',...
    'keyConditions', {@(y) subsref(isstrprop(y, 'lower'), first)},...
    'valueConditions', {@(y) ~isstruct(y)});

metricFields = @(x) util.struct.matchingfields(x,...
    'logic', 'or',...
    'keyConditions', {@(y) all(isstrprop(y, 'lower'))},...
    'valueConditions', {@(y) isstruct(y)});

comparison = struct();

% Iterate the splits
splits = tuning_curves(sarel);
for split = setdiff(splits,"shuffle")

    S = sarel.(split);
    metrics = intersect(metricFields(S), Opt.metrics);

    for tuning_curve = tuning_curveFields(S)
        comparison.(tuning_curve) = coding.sarel.shuffle.stat(sarel, [split, tuning_curve]);
        for metric = metrics
            if ~isfield(S.(metric), tuning_curve); continue; end
            comparison.(metric).(tuning_curve) = coding.sarel.shuffle.stat(sarel, [split, metric, tuning_curve]);
        end
    end
end
