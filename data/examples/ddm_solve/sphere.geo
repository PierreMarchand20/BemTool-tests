Mesh.MshFileVersion = 2.2;
SetFactory("OpenCASCADE");
h=0.07;
Mesh.CharacteristicLengthMin = h;
Mesh.CharacteristicLengthMax = h;
Sphere(1) = {0,0,0,1};
