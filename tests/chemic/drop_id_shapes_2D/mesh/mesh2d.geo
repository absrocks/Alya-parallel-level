//+
Point(1) = {-1, -1, 0, 0.04};
//+
Point(2) = {1, -1, 0, 0.04};
//+
Point(3) = {1, 1, 0, 0.04};
//+
Point(4) = {-1, 1, 0, 0.04};
//+
Line(1) = {4, 3};
//+
Line(2) = {3, 2};
//+
Line(3) = {2, 1};
//+
Line(4) = {1, 4};
//+
Line Loop(1) = {1, 2, 3, 4};
//+
Plane Surface(1) = {1};
//+ Switch mesh orientation:
Reverse Surface(1);
//+
Physical Line("Inlet") = {4};
//+
Physical Line("Outlet") = {2};
//+
Physical Line("Walls") = {1, 3};
//+
Physical Surface("Fluid") = {1};
