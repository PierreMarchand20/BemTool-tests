SetFactory("OpenCASCADE");
Mesh.MshFileVersion = 2.2;
h=0.005;
Mesh.CharacteristicLengthMin = h;
Mesh.CharacteristicLengthMax = h;
Rectangle(1) = {-2,-2,0,4,4};
Rectangle(2) = {-1,-1,0,2,2};
BooleanUnion{Line{2}; Delete;}{Line{1}; Delete;}
