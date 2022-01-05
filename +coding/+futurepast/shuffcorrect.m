function D = shuffcorrect(D)
% 
% ------
% Inputs
% ------
%
% D : struct
%
% D.MI
% D.shuffles
%
% --------------------
% Output modifications
% --------------------
%
% D.MIcorrected
% D.MIcorrectedDistribution

ip = inputParser;
ip.addParameter('normalization', 'divide');
ip.addParameter('shufflereduce', @median);


% Actual mutual information
MI = D.MI;

% Shuffle mutual information
MIshuffles = cellfun(@(x) x.MI, D.shuffles, "UniformOutput", false);
MIshuffles = cat(1, MIshuffles{:});
%MIshuffles = Opt.shufflereduce(MIshuffles, 1);

switch lower(char(Opt.normalization))
case 'divide'
    D.MIcorrected = MI ./ Opt.shufflereduce(MIshuffles);
    D.MIcorrectedDistribution = MI ./ MIshuffles;
    D.
case 'subtract'
    D.MIcorrected = MI - Opt.shufflereduce(MIshuffles);
    D.MIcorrectedDistribution = MI - MIshuffles;
end

