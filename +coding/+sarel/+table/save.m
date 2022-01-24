function save(tabStruct, varargin)
% Saves tabular goal-vector analysis

ip = inputParser;
ip.addParameter('target', []); % {matlab} | csv | {all}
ip.addParameter('tag', []); % {matlab} | csv | {all}
ip.parse(varargin{:});
Opt = ip.Results;

datafolder = coding.file.datafolder('exp_pro');

if isempty(Opt.target) || strcmp(Opt.target, 'all')
    Opt.target = {'matlab', 'csv'};
end
if ~iscell(Opt.target)
    Opt.target = {Opt.target};
end
if ~isempty(Opt.tag)
    Opt.tag = ['=e' Opt.tag];
end

for target = Opt.target
    switch Opt.target{1}
        case 'matlab'
            save(fullfile(datafolder, 'goal-vector_matlab%s.mat', Opt.tag), '-struct', 'tabStruct');
        case 'csv'
            for split = string(fieldnames(tabStruct))'
                for derivative = string(fieldnames(tabStruct.(split)))'
                    writetable(tabStruct.(split).(derivative), ...
                        fullfile(datafolder, sprintf('goal-vector_%s_%s%s.csv', split, derivative, Opt.tag)));
                end
            end
    end
end
