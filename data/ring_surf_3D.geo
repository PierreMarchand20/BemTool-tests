SetFactory("OpenCASCADE");
h=0.2;
Mesh.CharacteristicLengthMin = h;
Mesh.CharacteristicLengthMax = h;
Sphere(1) = {0,0,0,1};
Sphere(2) = {0,0,0,2};
BooleanUnion{Surface{2}; Delete;}{Surface{1}; Delete;}
