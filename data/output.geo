SetFactory("OpenCASCADE");
h=0.05;
Mesh.CharacteristicLengthMin = h;
Mesh.CharacteristicLengthMax = h;
Disk(1) = {0,0,0,1.2};
Disk(2) = {0,0,0,4};
BooleanDifference{Surface{2}; Delete;}{Surface{1}; Delete;}
