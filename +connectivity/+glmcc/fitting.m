function out = fitting(spikes, varargin)
% Wrapper for glmcc_fitting.py

ip = inputParser;
ip.addParameter('resultLoc', []);
ip.parse(varargin{:})
Opt = ip.Results;

nCells = numel(spikes.spikeTimes);
WFile  = '/tmp/wfile';


currdir = pwd;
destruct = onCleanup(@() cd(currdir));
cd('/tmp/');

% Setup the inputs
connectivity.glmcc.makecelltmpfiles(spikes);

% Run GLMcc
outfolder = '/tmp/glmcc';
if ~exist(outfolder,'dir')
    mkdir(outfolder)
end
if ~exist(fullfile(outfolder,'weights'),'dir')
    mkdir(fullfile(outfolder,'weights'))
end
Wfile = '/tmp/glmcc/sorted_w.csv'
command = join(["!python -m pdb ~/Code/analysis/+connectivity/+glmcc/GLMCC/glmcc_fitting.py", nCells, "/tmp/cellFiring", "exp", Wfile, "all", "LR"]," ");
eval(command)

% Process the weight file



if Opt.resultLoc
end
