SetFactory("OpenCASCADE");
h = 0.04;
Mesh.CharacteristicLengthMin = h;
Mesh.CharacteristicLengthMax = h;
Sphere(1) = {0,0,0,0.5};