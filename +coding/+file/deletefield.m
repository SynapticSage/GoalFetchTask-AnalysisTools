function deletefield(animal, index, fields, varargin)

filename = coding.file.filename(animal, index);
fields = string(fields);
for field = fields(:)'
    util.matfile.rmvar(filename, field);
end
