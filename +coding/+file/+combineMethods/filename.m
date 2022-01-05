function filename_full = filename(animal, index)
% function filename_full = filename(animal, index)
% 

filename = sprintf("checkpoint_days=%s", join(string(index),","));
filename_full = fullfile(folder, filename);

