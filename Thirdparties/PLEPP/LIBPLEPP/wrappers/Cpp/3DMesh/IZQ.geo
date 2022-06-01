Point(1) = {-1, 0, 0, 1.0};
Point(2) = {0, 0, 0, 1.0};
Point(3) = {0, 1, 0, 1.0};
Point(4) = {-1, 1, 0, 1.0};
Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 1};
Line Loop(5) = {4, 1, 2, 3};
Plane Surface(6) = {5};
Extrude {0, 0, 1} {
  Surface{6};
}
Mesh.CharacteristicLengthFactor = 0.3;
Surface Loop(29) = {28, 15, 6, 19, 23, 27};
Volume(30) = {29};
Physical Volume("BULK") = {1};
