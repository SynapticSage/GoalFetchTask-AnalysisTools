animal = "RY16";

dayepochs = ndb.indicesMatrixForm(behavior);
days = unique(dayepochs(:,1));
for day = days'

    codingStruct = coding.coordcentric(animal, day); 

end

