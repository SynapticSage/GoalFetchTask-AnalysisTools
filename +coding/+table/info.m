function tab = info(animal, datatype, varargin)
% function tab = info(animal, datatype, varargin)
% Creates an easier to work with table of cellinfo

ip = inputParser;
ip.addParameter('maximizeFields', false); % only accept entries with all fields? (prioritizes days that have complete preprocessing)
ip.addParameter('filterByFields', @(x) ismember("area", string(fieldnames(x))));
ip.addParameter('inds', []); % if inds, just compute what we need and add it to the table
ip.parse(varargin{:})
Opt = ip.Results;

if nargin < 2
    datatype = "cellinfo";
end

disp("Loading");
info = ndb.load(animal, datatype);
%originalIndices = ndb.indicesMatrixForm(info);

% -----
% Label
% -----
disp("Labeling dimensions");
dimLabels = ["day","epoch","tetrode","cell"];
if ~isempty(Opt.inds)
    infoTable = ndb.load(animal, string(datatype) + "Table");

    pullTheseLabels = dimLabels(1:numel(Opt.inds));
    addTheseLater{:,1} = pullTheseLabels;
    addTheseLater{:,2} = Opt.inds;
    info = ndb.get(info, Opt.inds);

    infoQuery = "$" + pullTheseLabels + " == " + string(Opt.inds);
    removeTheseRows = util.table.query_logical(infoTable, infoQuery);
    infoTable(find(removeTheseRows), :) = [];

    topDim = size(ndb.indicesMatrixForm(info), 2);
    remainingDimLabels = dimLabels(numel(Opt.inds)+1:numel(Opt.inds)+topDim);
    remainingDims = 1:numel(remainingDimLabels);
    [dims, dimLabels] = deal(remainingDims, remainingDimLabels);
else
    topDim = size(ndb.indicesMatrixForm(info), 2);
    dims = 1:topDim;
    dimLabels = dimLabels(1:topDim);
    infoTable = [];
end
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

% ---------------
% Convert to tidy
% ---------------
disp("Converting to tidy data");
tab = tidyData.fromNdb(info);
assert(all(ismember(["area"], string(fieldnames(tab)))), "Area should be a field")

% ---------------
% Append to existing or is this it?
% ---------------
if isempty(infoTable)
    infoTable = tab;
else
    for row = addTheseLater'
        tab.(row{1}) = repmat(row{2}, height(tab), 1);
    end
    labels    = intersect(["day", "epoch", "tetrode", "cell"],...
                         string(infoTable.Properties.VariableNames));
    infoTable = util.cell.icat({infoTable, tab}, 'fieldCombine', 'union');
    infoTable = sortrows(infoTable, labels);
end

% Save
% ----
datadir = animaldef(animal);
tableName = datatype + "Table";
assign(tableName, infoTable);
disp("Saving");
save(fullfile(datadir{2}, string(animal) + tableName), tableName);
