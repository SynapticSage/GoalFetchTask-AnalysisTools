function D = estOptTau(D)
% Estimates optimall tau time shift
% 
% Inputs
% ---------
%
% D : struct
%   D.MI
%   D.MIcorrected
%   D.MIcorrectedDistribution
% 
% Output Modifications
% --------------------
%
% D.bestTau, best time shift
% D.bestMI, best time shift
% D.bestTauInd, index of best tau
% D.tauImprovement, how much improvement over zero shift
% D.tauImprovemntCI
% D.optimalTable, stores all of values per cell in a table


% Determine the best time shift
[D.bestMI, D.bestTauInd] = max(D.MIcorrected);
D.bestTau = D.shifts(D.bestTauInd);

% How much does this tau shift improve the MI above the null assumption of time_shift = 0?
D.tauImprovement = D.bestMI / D.MIcorrected(D.shifts==0);

% We want to get a handle into how much this tauImproment is. Is it
% significant? Ways to get at this: 
%
% (1) Compare to deviation sizes - standard
% deviation of MI vals
%
% (2) Distribution approach - 
%  Distribution of (bestTau, bestMI) 
%  in MIcorrectedDistribution

[D.distMI, D.distTauInd] = max(D.MIcorrectedDistribution);
D.distTau = D.shifts(D.distTauInd);


