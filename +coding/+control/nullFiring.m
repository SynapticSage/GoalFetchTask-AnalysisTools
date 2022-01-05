aunction nullFR = nullFiring(FR)
% Creates a nullfiring distribution.

ip = inputParser;
ip.addParameter('animal', []); % Can be a string for some of the simpler null distributions or can be a struct for some fo the more complex ones
ip.addParameter('method', "ones"); % Can be a string for some of the simpler null distributions or can be a struct for some fo the more complex ones
ip.parse(varargin{:})
Opt = ip.Results;

%Aliases 
NEURON = 2;
TIME = 1;

% Create Null Data
if isstring(method) || ischar(method)
    switch char(method)
        case 'ones'
            nullFR = ones(size(FR));
        case 'uniform'
            sample = rand(size(FR));
        case 'normal'
            sample = randn(size(FR));
        case 'shuffle_time'
            reorder = randperm(1:size(FR,TIME));
            nullFR = FR(reorder, :);
        case 'shuffle_neuron'
            reorder = randperm(1:size(FR,NEURON));
            nullFR = FR(:, reorder);
        case 'shuffle'
            reorder = randperm(1:numel(FR));
            nullFR  = reshape(FR(reorder), size(FR));
        otherwise
            error("Bad method");
    end
elseif isstruct(method)
    switch char(method)
        case 'placefield'
        case 'headdirfield'
        otherwise
            error("Bad method");
    end
else
    error("Unrecognized type!");
end
