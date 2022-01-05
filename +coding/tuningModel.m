function estFiringRates = tuningModel(propertyTuningCurve, binsTuningCurve, propertyTimeSeries)
% Estimate firing rates

isin = propertyTimeSeries >= binsTuningCurve(:,1) & propertyTimeSeries < binsTuningCurve(:,2);
[~, locs] = max(isin); % Find 1s
estFiringRates = propertyTuningCurve(locs);
