Mesh.MshFileVersion = 2.2;

lc = 0.005;
L= 1;
R = 1;

Point(1) = {0, 0, 0, lc};
Point(2) = {L, 0,  0, lc} ;
Point(3) = {L, L/2, 0, lc} ;
Point(4) = {L/2,  L/2, 0, lc/10.} ;
Point(5) = {L/2,  L, 0, lc} ;
Point(6) = {0,  L, 0, lc} ;

Line(1) = {1,2} ;
Line(2) = {2,3} ;
Line(3) = {3,4} ;
Line(4) = {4,5} ;
Line(5) = {5,6} ;
Line(6) = {6,1} ;

Curve Loop(1) = {6,1,2,3,4,5} ;

Plane Surface(1) = {1} ;