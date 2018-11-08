SetFactory("OpenCASCADE");
h=0.05;
Mesh.CharacteristicLengthMin = h;
Mesh.CharacteristicLengthMax = h;
Rectangle(2) = {-1.9,-1.9,0,3.8,3.8};
Rectangle(1) = {-1.1,-1.1,0,2.2,2.2};
BooleanDifference{Surface{2}; Delete;}{Surface{1}; Delete;}
