function tab = info(animal, datatype, varargin)
% function tab = info(animal, datatype, varargin)
% Creates an easier to work with table of cellinfo

ip = inputParser;
ip.addParameter('maximizeFields', false); % only accept entries with all fields? (prioritizes days that have complete preprocessing)
ip.addParameter('filterByFields', @(x) ismember("area", string(fieldnames(x))));
ip.parse(varargin{:})
Opt = ip.Results;

if nargin < 2
    datatype = "cellinfo";
end

disp("Loading");
info = ndb.load(animal, datatype);
originalIndices = ndb.indicesMatrixForm(info);

% -----
% Label
% -----
disp("Labeling dimensions");
inds = ndb.indicesMatrixForm(info);
topDim = size(inds, 2);
dims = 1:topDim;
dimLabels = ["day","epoch","tetrode","cell"];
dimLabels = dimLabels(1:topDim);
info = ndb.dimLabel(info, dims, dimLabels);

% --------------
% Filter results
% --------------
if Opt.maximizeFields
    fieldcount = ndb.apply(info, @(x) numel(fieldnames(x)), "outputList");
    fieldcount = cat(1, fieldcount{:});
    m = max(fieldcount(:));
    filt = fieldcount == m;
    info = ndb.query(info, filt);
end

if isempty(Opt.filterByFields)
    fieldcount = ndb.apply(info, Opt.filterByFields, "outputList");
    fieldcount = cat(1, fieldcount{:});
    m = max(fieldcount(:));
    filt = fieldcount == m;
    info = ndb.query(info, filt);
end

% Convert to tidy
disp("Converting to tidy data");
tab = tidyData.fromNdb(info);
assert(all(ismember(["area"], string(fieldnames(tab)))))

% Save
% ----
datadir = animaldef(animal);
tableName = string(animal) + datatype + "Table";
assign(tableName, tab);
disp("Saving");
save(fullfile(datadir{2}, tableName), tableName);
