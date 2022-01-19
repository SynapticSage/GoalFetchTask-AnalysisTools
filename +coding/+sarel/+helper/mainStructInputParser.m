function ip = mainStructInputParser()
% Contains the optional arguments used in functions that process details of the
% sarel structs (structs containing sarel paper style results)
%
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


ip = inputParser;
ip.addParameter('metrics', ["vm", "occNorm", "rayleigh", "maxmean_indices"]);
