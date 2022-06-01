Point(1) = {0.0, 0.0, 0.0, 1.0};
Point(2) = {1.0, 0.0, 0.0, 1.0};
Point(3) = {1.0, 1.0, 0.0, 1.0};
Point(4) = {0.0, 1.0, 0.0, 1.0};

Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 1};

Line Loop(1) = {1,2,3,4};
Surface(1) = {1};

Transfinite Curve {1, 2, 3, 4} = 5 Using Progression 1;
Transfinite Surface {1} Alternated;

//Field[1] = BoundaryLayer;
//Field[1].EdgesList = {1,2,3,4};
////Field[1].ExcludedFaceList = {2};
//Field[1].AnisoMax = 1.0;
//Field[1].hfar = 0.25;
//Field[1].hwall_n = 0.25/4;
//Field[1].thickness = 0.25;
//Field[1].ratio = 1.1;
//Field[1].Quads = 0;
//Field[1].IntersectMetrics = 1;
//BoundaryLayer Field = 1;

Extrude {0, 0, 1} {
  Surface{1};
  Layers{4}; 
}

Physical Surface("shear") = {21};
Physical Surface("wall") = {26, 1, 13, 17, 25};
Physical Volume("fluid") = {1};

Mesh.MshFileVersion = 2.2;
