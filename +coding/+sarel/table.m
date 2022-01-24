function Out = table(sarel, varargin)
% Creates a tidy data table about goal vector analyses
%
% needed for ggplot-spirited visuals
%
%
%  ---------------------
%  Input: *Sarel struct*
%  ---------------------
%  | Levels     | Purpose                                     |   |
%  |------------+---------------------------------------------+---|
%  | splits     | ways of splitting time length               |   |
%  |            | to measure tuning curves                    |   |
%  |            | (splits contain a sub-dimension of          |   |
%  |            | dimension, which all fields below have)     |   |
%  |            |                                             |   |
%  | datafields | the main tuning curve types                 |   |
%  |            |                                             |   |
%  | metrics    | *measurements* done on tuning curves        |   |
%  |            |                                             |   |
%  | stats      | further measures done either on metrics     |   |
%  |            | or main/shuffle of tuning curves or metrics |   |
%  
%
%  ---------------------
%  Output: *TIDY DATA*
%  ---------------------
%  For the above, the decision is how to lay these out into tidy data
%  structure.
%  |--------------- INDEX -------------|----------------DATA COLS---------------|
%  v                                   v                                        v
%  | dim1  | ... | dimD | tuning_curve | metric1_stat1   | ... | metricM_statN  |
%
%
%  And for each split, we output one table. These are placed into a variable
%  as Out.split
%
%  -----
%  Notes
%  -----
%  This method shares much of its logic with sarel.shuffle.compare

ip = coding.sarel.helper.mainStructInputParser();
ip.parse(varargin{:});
Opt = ip.Results;

first = struct('type', '()', 'subs', {{1}});
%tuning_curveFields = @(x) util.struct.matchingfields(x,...
%    'logic', 'and',...
%    'keyConditions', {@(y) subsref(isstrprop(y, 'lower'), first)},...
%    'valueConditions', {@(y) ~isstruct(y)});
tuning_curveFields = @(x) util.struct.matchingfields(x,...
    'logic', 'and',...
    'keyConditions', {@(x) ismember(x, coding.sarel.table.field.tuningCurves)});
metricFields = @(x) util.struct.matchingfields(x,...
    'logic', 'and',...
    'keyConditions', {@(y) subsref(isstrprop(y, 'lower'), first)},...
    'valueConditions', {@(y) isstruct(y)});
splitFields = metricFields;

Out = struct();
splits = splitFields(sarel);
for split = setdiff(splits, "shuffle")

    S = sarel.(split);
    metrics = intersect(metricFields(S), Opt.metrics);

    for tuning_curve = tuning_curveFields(S)
        clear Metrics
        Out.(split).(tuning_curve) = table();

        % Add any raw level data
        Metrics.raw = coding.sarel.table.stat(sarel, [split, tuning_curve]);
        descriptors = struct('tuning_curve', tuning_curve);

        % Add metric level data
        for metric = metrics
            if ~isfield(S.(metric), tuning_curve); continue; end
            Metrics.(metric) = coding.sarel.table.stat(sarel,...
                [split, metric, tuning_curve]);

        end

        %if numel(S.Dimensions) > 2; keyboard; end
        descriptors.Dimensions = coding.sarel.table.dimension(Metrics.raw, ...
            sarel.Binning, S.Dimensions, tuning_curve);

        Out.(split).(tuning_curve) = coding.sarel.table.append(Out.(split).(tuning_curve), Metrics, descriptors);
    end
end
