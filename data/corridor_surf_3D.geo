SetFactory("OpenCASCADE");
h=0.399;
Mesh.CharacteristicLengthMin = h;
Mesh.CharacteristicLengthMax = h;
Block(1) = {-2,-2,-2,4,4,4};
Block(2) = {-1,-1,-1,2,2,2};
BooleanUnion{Surface{2}; Delete;}{Surface{1}; Delete;}
