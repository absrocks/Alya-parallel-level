// -------------------------------------------------------------------------
// Minimal flame case
// -------------------------------------------------------------------------
//
//         4     3      3
//          0=========0
//          |         |
//          |         |
//         4|         |2
//          |         |
//          |         |
//          0=========0
//         1     1     2
//


gx  = 1.0;        
gy  = 1.0;       
gz  = 1.0;       

nGapx = 3; // DAMN  COARSE
nGapy = 3;
nGapz = 3;

z1 = 0.0;

// Points
Point(1)  = { -gx/2  , -gy/2  , -gz/2, z1}; 
Point(2)  = {  gx/2  , -gy/2  , -gz/2, z1}; 
Point(3)  = {  gx/2  ,  gy/2  , -gz/2, z1}; 
Point(4)  = { -gx/2  ,  gy/2  , -gz/2, z1}; 
                   
Point(5)  = { -gx/2  , -gy/2  ,  gz/2, z1}; 
Point(6)  = {  gx/2  , -gy/2  ,  gz/2, z1}; 
Point(7)  = {  gx/2  ,  gy/2  ,  gz/2, z1}; 
Point(8)  = { -gx/2  ,  gy/2  ,  gz/2, z1}; 

// Lines
Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 1};

Line(5) = {5, 6};
Line(6) = {6, 7};
Line(7) = {7, 8};
Line(8) = {8, 5};

Line(9)  = {1, 5};
Line(10) = {2, 6};
Line(11) = {3, 7};
Line(12) = {4, 8};

// Line loops
Line Loop(1) = {1,2,3,4};
Line Loop(2) = {-5,-8,-7,-6};
Line Loop(3) = {-1,9,5,-10};
Line Loop(4) = {-2,10,6,-11};
Line Loop(5) = {-3,11,7,-12};
Line Loop(6) = {-4,12,8,-9};

// Surface definition
Plane Surface(1) = {1};
Plane Surface(2) = {2};
Plane Surface(3) = {3};
Plane Surface(4) = {4};
Plane Surface(5) = {5};
Plane Surface(6) = {6};

// Surface loop
Surface Loop(666) = {1,2,3,-4,-5,-6};

// Volume
Volume(1) = {666};


//  Transfinite lines and surfaces for structured grids 
Transfinite Line{1,3,5,7}     =  nGapx;
Transfinite Line{2,6,8,4}     =  nGapy;
Transfinite Line{9,10,11,12}  =  nGapz;
Transfinite Surface{1,2,3,4,5,6};
Transfinite Volume{1};

// Recombination to get QUADS instead of TRIA
Recombine Surface {1,2,3,4,5,6};
Recombine Volume {1};

// Boundary conditions
Physical Surface("Left") = {6};
Physical Surface("Right") = {4};
Physical Surface("Wall") = {1,2,3,5};
Physical Volume("Fluid") = {1};

