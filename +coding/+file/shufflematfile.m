function m = shufflematfile(animal, index, varargin)

filename_full = coding.file.shufflematfilename(animal, index, varargin{:});
if ~(numel(varargin) == 1 && isstruct(varargin{1}))
    S = util.struct.varargin2struct(varargin);
else
    S = varargin{1};
end
if ~isfield(S,'writable')
    m = matfile(filename_full);
else
    m = matfile(filename_full, 'Writable', S.writable);
end

