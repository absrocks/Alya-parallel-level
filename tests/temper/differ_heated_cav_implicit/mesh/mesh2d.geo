
l  = 1;  
Ne = 10;
Np = Ne+1;


z1 = 0.0;

// Corners:
Point(1)  = {-l/2   , -l/2   , z1}; 
Point(2)  = {l/2    , -l/2   , z1}; 
Point(3)  = {l/2    , l/2    , z1}; 
Point(4)  = {-l/2   , l/2    , z1}; 


// Lines
Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 1};


// Line loops
Line Loop(1) = {1,2,3,4};

// Surface definition
Plane Surface(51) = {1};

//  Transfinite lines and surfaces for structured grids
Transfinite Line{1,2,3,4} = Np;
Transfinite Surface{51};

// Recombination to get QUADS instead of TRIA
Recombine Surface {51};

// Boundary conditions
Physical Line("Bottom")  = {1};
Physical Line("Right")   = {2};
Physical Line("Top")     = {3};
Physical Line("Left")    = {4};
Physical Surface("Surf") = {51};

