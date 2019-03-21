Mesh.MshFileVersion = 2.2;
SetFactory("OpenCASCADE");
h=0.01;
Mesh.CharacteristicLengthMin = h;
Mesh.CharacteristicLengthMax = h;
Circle(1) = {0,0,0,1};
