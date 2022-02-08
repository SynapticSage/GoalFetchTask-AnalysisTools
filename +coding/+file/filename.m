function [filename_full, Props] = filename(animal, index, varargin)
% function filename_full = filename(animal, index)

ip = inputParser;
ip.KeepUnmatched = true;
ip.addParameter('datafolder', {'exp_raw'});
ip.parse(varargin{:})
Opt = ip.Results;

if iscell(Opt.datafolder)
    Opt.datafolder = coding.file.datafolder(Opt.datafolder{:});
elseif ~isempty(Opt.datafolder)
    Opt.datafolder = Opt.datafolder;
else
    Opt.datafolder = coding.file.datafolder();
end

folder = fullfile(Opt.datafolder);
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
