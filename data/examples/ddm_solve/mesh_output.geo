Mesh.MshFileVersion = 2.2;
SetFactory("OpenCASCADE");
h=0.05;
Mesh.CharacteristicLengthMin = h;
Mesh.CharacteristicLengthMax = h;
Disk(1) = {0,0,0,1.1};
Disk(2) = {0,0,0,3};
BooleanDifference{Surface{2};Delete;}{Surface{1};Delete;}