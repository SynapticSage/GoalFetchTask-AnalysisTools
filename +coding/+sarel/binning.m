function [Binning, Opt]  = binning(beh, varargin)
% Bins for the sarel functions

ip = inputParser;
ip.KeepUnmatched = true;
ip.addParameter('anglebins', 18);
ip.addParameter('distbins',  20);
ip.addParameter('xy_grid', 22);
ip.addParameter('anglefield', 'currentAngle');
ip.addParameter('distancefield', 'currentEucDist');
ip.parse(varargin{:});
Opt = ip.Results;

% Create the grids of angle and distance to measure
% -------------------------------------------------
% Binning scheme
Binning.angleEdges      = linspace(min(beh.(Opt.anglefield)),  max(beh.(Opt.anglefield)),  Opt.anglebins+1);
Binning.distEdges       = linspace(min(beh.(Opt.distancefield)), max(beh.(Opt.distancefield)), Opt.distbins+1);
Binning.angleCenters    = mean([Binning.angleEdges(1:end-1)', Binning.angleEdges(2:end)'], 2)';
Binning.distCenters     = mean([Binning.distEdges(1:end-1)', Binning.distEdges(2:end)'] ,2)';
Binning.possibleStops   = unique(beh.stopWell);
Binning.possibleStarts  = unique(beh.stopWell);
Binning.possibleCuemem  = unique(beh.cuemem);


% Other anlyses that require place binning
% ----------------------------------------
scaffold = coding.field.scaffold(beh, 'props', ["x","y"], 'grid', 22);
Binning.xCenters = scaffold(1).center;
Binning.yCenters = scaffold(2).center;
Binning.xEdges = scaffold(1).edge;
Binning.yEdges = scaffold(2).edge;
