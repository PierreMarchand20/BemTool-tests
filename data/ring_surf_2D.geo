SetFactory("OpenCASCADE");
h=0.0001;
Mesh.CharacteristicLengthMin = h;
Mesh.CharacteristicLengthMax = h;
Circle(1) = {0,0,0,1};
Circle(2) = {0,0,0,2};
BooleanUnion{Line{2}; Delete;}{Line{1}; Delete;}
