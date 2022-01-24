function Stat = stat(X, mainFields, higherFields)

if nargin < 3
    higherFields = string([]);
end

main    = nd.nestedFieldCat(X, [mainFields, higherFields]);
Stat = struct();
if ~isstruct(main) && isnumeric(main)
    Stat = main;
elseif isstruct(main)
    for field = string(fieldnames(main))'
        newfield = join([higherFields,field], "_");
        Stat.(newfield) =...
            coding.sarel.table.stat(X, mainFields, [higherFields, field]);
        %if isstruct(Stat.(newfield))
        %    error('Stat.(%s) is a freaking struct', newfield);
        %end
    end
else
    warning('Unprocessable %s', join(mainFields, '.'))
end

if isstruct(Stat)
    Stat = nd.fullyUnnest(Stat);
end

