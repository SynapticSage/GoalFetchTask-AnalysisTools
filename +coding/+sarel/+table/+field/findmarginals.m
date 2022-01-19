function marginals = findmarginals(fieldname)

marginal_names = ["starts", "stops", "cuemem", "stop", "start"];

locs = arrayfun(@(x) strfind(fieldname, x), marginal_names, 'UniformOutput', false);
locs = cellfun(@min, locs, 'UniformOutput', false);
marginal_names(cellfun(@isempty, locs)) = [];
locs = cat(1, locs{:});
marginal_names = marginal_names(:);

tab = table(locs, marginal_names);
tab = sortrows(tab, 'locs');

marginals = tab.marginal_names;
