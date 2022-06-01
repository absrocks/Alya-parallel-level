L = 2.0;
H = 1.0;
E = 0.0;
nx = 100;
h = L/nx;

Point(1) = {-L/2,-H/2,-E/2,h};
Point(2) = { L/2,-H/2,-E/2,h};
Point(3) = { L/2, H/2,-E/2,h};
Point(4) = {-L/2, H/2,-E/2,h};

Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 1};

Line Loop(1) = {1,2,3,4};
Surface(1) = {1};

Transfinite Surface {1} Alternated;
Recombine Surface {1:8};

Physical Curve("wall") = {3, 1};
Physical Surface("fluid") = {1};

Mesh.MshFileVersion = 2.2;

