Point(1) = {0, 0, 0, 1.0};
Point(2) = {1, 0, 0, 1.0};
Point(3) = {1, 1, 0, 1.0};
Point(4) = {0, 1, 0, 1.0};

Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 1};

Line Loop(1) = {1, 2, 3, 4};
Surface(1) = {1};

Extrude {0, 0, 1} {
  Surface{1};
  Layers{4};
  Recombine; 
}

Transfinite Line {1, 2, 3, 4} = 5 Using Progression 1;
Transfinite Surface {1};

//+
Physical Surface("top") = {21};
//+
Physical Surface("bottom") = {13};
//+
Physical Surface("front") = {26};
//+
Physical Surface("back") = {1};
//+
Physical Surface("left") = {17};
//+
Physical Surface("right") = {25};
//+
Physical Volume("fluid") = {1};
