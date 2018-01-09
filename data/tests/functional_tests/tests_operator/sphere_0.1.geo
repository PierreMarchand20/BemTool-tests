SetFactory("OpenCASCADE");
h = 0.1;
Mesh.CharacteristicLengthMin = h;
Mesh.CharacteristicLengthMax = h;
Sphere(1) = {0,0,0,0.5};
