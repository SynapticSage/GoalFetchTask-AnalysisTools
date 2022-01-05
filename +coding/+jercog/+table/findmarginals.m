function marginals = findmarginals(S)

common_names  = intersect(coding.jercog.table.field.standard,       fieldnames(S));
non_marginals = intersect(coding.jercog.table.field.nonmarginal(S), fieldnames(S));
marginals = setdiff(setdiff(fieldnames(S), common_names), non_marginals);
marginals = marginals(:)';
