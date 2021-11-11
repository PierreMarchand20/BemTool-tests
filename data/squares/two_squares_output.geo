SetFactory("OpenCASCADE");
Mesh.MshFileVersion = 2.2;

DefineConstant[
    npplo = {10, Name "npplo"},
    k = {5, Name "k"}
];
Printf("npplo : %f",npplo);
Printf("k : %f",k);

meshsize=2*Pi/(npplo*k);
Mesh.CharacteristicLengthMin = meshsize;
Mesh.CharacteristicLengthMax = meshsize;
Printf("Meshsize : %f",meshsize);

Rectangle(1)={-1.1,-0.1,0,1.2,1.2};
Rectangle(2)={Pi-0.1,-0.1,0,1.2,1.2};
Rectangle(3)={-2,-1,0,4+Pi,3};

BooleanDifference{Surface{3};Delete;}{Surface{1};Delete;}
BooleanDifference{Surface{3};Delete;}{Surface{2};Delete;}