parquetFolder = '/Volumes/FastData/goal/RY16_days=36_parquet/';
matfileLocation = '/Volumes/FastData/goal/RY16_days=36_shuffle.mat';

M = matfile(matfileLocation, 'Writable', true);
util.matfile.rmvar(matfileLocation, 'shuffle');

iCount = 0;
for file = progress(dir(fullfile(parquetFolder, '*.parquet'))', 'Title', 'Bringing parquet files into matfile')

    iCount = iCount + 1;

    path = fullfile(file.folder, file.name);
    P = parquetread(path);
    P = util.type.castefficient(P);

    M.shuffle(iCount, 1) = table2struct(P, 'ToScalar', true);

end
