function Out = table(sarel)
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

first = struct('type', '()', 'subs', {{1}});
tuning_curveFields = @(x) util.struct.matchingfields(x,...
    'logic', 'or',...
    'keyConditions', {@(y) subsref(isstrprop(y, 'lower'), first)},...
    'valueConditions', {@(y) ~isstruct(y)});
metricFields = @(x) util.struct.matchingfields(x,...
    'logic', 'or',...
    'keyConditions', {@(y) all(isstrprop(y, 'lower'))},...
    'valueConditions', {@(y) isstruct(y)});

splits = tuning_curves(sarel);
for split = setdiff(splits, "shuffle")

    S = sarel.(split);
    metrics = intersect(metricFields(S), Opt.metrics);

    Out.(split) = table();
    for tuning_curve = tuning_curveFields(S)

        % Add any raw level data
        Metrics.raw = coding.sarel.shuffle.stat(sarel, [split, tuning_curve]);
        descriptors = struct('metric', 'raw', 'tuning_curve', tuning_curve);

        % Add metric level data
        for metric = metrics
            if ~isfield(S.(metric), tuning_curve); continue; end
            Metrics.(metric) = coding.sarel.table.stat(sarel,...
                [split, metric, tuning_curve]);
        end

        descriptors.Dimensions = util.table.dimension(Metrics.raw, ...
            sarel.Binning, S.Dimensions);

        Out.(split) = coding.sarel.table.append(Out.(split), Metrics, descriptors);
    end
end
