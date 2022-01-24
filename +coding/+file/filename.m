function [filename_full, Props] = filename(animal, index, varargin)
% function filename_full = filename(animal, index)
% 

ip = inputParser;
ip.KeepUnmatched = true;
ip.addParameter('projectFolder', []);
ip.parse(varargin{:})
Opt = ip.Results;

if isempty(Opt.projectfolder)
    Opt.projectfolder = coding.file.projectfolder();
end

folder = fullfile(Opt.projectfolder);
if ~exist(folder, 'dir')
    mkdir(folder);
end

header = string(animal); 
desc = sprintf("days=%s", join(string(index),","));
ext = ".mat";
headerdesc = join([header, desc], "_");
filename = headerdesc + ext;
filename_full = fullfile(folder, filename);

Props = struct();
Props.animal = animal;
Props.index  = index;
Props.folder =  folder;
Props.header = header;
Props.desc = desc;
Props.headerdesc = headerdesc;
Props.ext = ext;
