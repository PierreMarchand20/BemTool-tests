lc = 0.00075;
Point(0) = { 0 , 0 , 0 , lc};
Point(1) = { 0.5 , 0 , 0 , lc};
Point(2) = { 0 , 0.5 , 0 , lc};
Point(3) = { -0.5 , 0 , 0 , lc};
Point(4) = { 0 , -0.5 , 0 , lc};
Circle(1) = { 1 , 0 , 2};
Circle(2) = { 2 , 0 , 3};
Circle(3) = { 3 , 0 , 4};
Circle(4) = { 4 , 0 , 1};
Physical Line(0) ={ 1 , 2 , 3 , 4};
