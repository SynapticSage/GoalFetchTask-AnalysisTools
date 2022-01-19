function Stat = stat(X, mainFields, higherFields)

if nargin < 3
    higherFields = string();
end

main    = nd.nestedFieldCat(X, [mainFields, higherFields]);
Stat = struct();
if ~isstruct(main) && isnumeric(main)
    Stat = main;
elseif isstruct(main)
    for field = string(fieldnames(main))'
        Stat.(join([higherFields,fields], "_")) =...
            coding.sarel.shuffle.stat(X, mainFields, [higherFields, field]);
    end
else
    warning('Unprocessable %s', join(mainFields, '.'))
end

Stat = nd.fullyUnnest(Stat);
