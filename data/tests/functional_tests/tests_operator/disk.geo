SetFactory("OpenCASCADE");
h=0.05;
Mesh.CharacteristicLengthMin = h;
Mesh.CharacteristicLengthMax = h;
Disk(1) = {0,0,0,1.9};
