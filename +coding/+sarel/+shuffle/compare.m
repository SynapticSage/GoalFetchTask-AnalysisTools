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
ip.addParameter('debug', false);
ip.KeepUnmatched = true;
ip.parse(varargin{:})
Opt = ip.Results;

first = struct('type', '()', 'subs', {{1}});
tuning_curveFields = @(x) util.struct.matchingfields(x,...
    'logic', 'and',...
    'keyConditions', {@(x) ismember(x, coding.sarel.table.field.tuningCurves)});
splitFields = @(x) util.struct.matchingfields(x,...
    'logic', 'or',...
    'keyConditions', {@(y) subsref(isstrprop(y, 'lower'), first)},...
    'valueConditions', {@(y) ~isstruct(y)});

metricFields = @(x) util.struct.matchingfields(x,...
    'logic', 'or',...
    'keyConditions', {@(y) all(isstrprop(y, 'lower'))},...
    'valueConditions', {@(y) isstruct(y)});

comparison = struct();

% Iterate the splits
splits = splitFields(sarel);
disp(splits);
for split = setdiff(splits,"shuffle")

    S = sarel.(split);
    metrics = intersect(metricFields(S), Opt.metrics);

    % Figure out shuffle dim across all metrics
    D = structfun(@ndims, coding.sarel.shuffle.stat(sarel, split, 'onlyShuffle', true));
    D = max(D);

    for tuning_curve = tuning_curveFields(S)
        C.Dimensions = nd.nestedFieldCat(sarel, [split, "Dimensions"]);
        assert(isstring(C.Dimensions));
        for metric = metrics
            if ~isfield(S.(metric), tuning_curve); continue; end
            %if split == "stops"; keyboard; end
            OptStat = util.struct.update(struct(), ip.Unmatched);
            OptStat.shuffleDim = D;
            C.(metric).(tuning_curve) = coding.sarel.shuffle.stat(sarel, ...
                [split, metric, tuning_curve], OptStat);
            %assert(isfield(C,metric) && isfield(C.(metric), tuning_curve));
            if util.struct.nestedisfield(C, {'vm', 'currentAngle_thetahat_varDifferences'})
                keyboard
            end
        end
        if util.struct.nestedisfield(C, {'vm', 'currentAngle_thetahat_varDifferences'})
            keyboard
        end
        Q = coding.sarel.shuffle.stat(sarel, [split, tuning_curve], 'shuffleDim', D);
        C.(tuning_curve) = Q;
        if util.struct.nestedisfield(C, {'vm', 'currentAngle_thetahat_varDifferences'})
            keyboard
        end
    end

    comparison.(split) = C;
    if Opt.debug
        util.struct.printstruct(comparison);
    end

end

comparison.Binning = sarel.Binning;
