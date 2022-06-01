pi = 3.14159265;
h = 2*pi/16;

Point(1) = {-pi,-pi,-pi,h};
Point(2) = { pi,-pi,-pi,h};
Point(3) = { pi, pi,-pi,h};
Point(4) = {-pi, pi,-pi,h};

Line(1) = {1,2};
Line(2) = {2,3};
Line(3) = {3,4};
Line(4) = {4,1};

Line Loop(1) = {1,2,3,4};
Surface(1) = {1};

Transfinite Surface {1};
Recombine Surface {1};

Extrude {0, 0, 2*pi} {
  Surface{1}; Layers{16}; Recombine;
}

Mesh.MshFileVersion = 2.2;
//+
Physical Volume("fluid") = {1};
