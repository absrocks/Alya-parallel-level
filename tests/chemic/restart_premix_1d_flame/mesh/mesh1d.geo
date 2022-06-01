


//                                            
//                                            
//    4          3             3               
//    o========================o              
//    |                        |                
//   4|                        |2               
//    |                        |                
//    o========================o              
//    1          1             2                
//                                            
//                                            
//                                            
//                                            
//                                            




Lx=0.0025;
Ly=0.00005;
nx=51;
ny=2;

Point(1) = {-Lx/2, -Ly/2, 0.0, 1.0};
Point(2) = {Lx/2,  -Ly/2, 0.0, 1.0};
Point(3) = {Lx/2,   Ly/2, 0.0, 1.0};
Point(4) = {-Lx/2,  Ly/2, 0.0, 1.0};

Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 1};
Transfinite Line{1,3} = nx;
Transfinite Line{2,4} = ny;
Line Loop(5) = {3, 4, 1, 2};
Plane Surface(6) = {5};
Transfinite Surface{6};
Recombine Surface {6};

Physical Line("LEFT") = {4};
Physical Line("RIGHT") = {2};
Physical Line("TB") = {3, 1};
Physical Surface("DOMAIN") = {6};
