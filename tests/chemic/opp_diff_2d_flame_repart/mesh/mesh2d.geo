// -------------------------------------------------------------------------
// Counterflow diffusion flame case
// -------------------------------------------------------------------------

ymin = -10e-3;
ymax = 20e-3;

g  = ymax-ymin;      // Gap distance
d1 = g*1.0;      // Streams diameter
l  = 5.0*g;      // Domain length

nGap = 16;
nIn  = nGap*d1/(1.5*g);
nLen = nGap*(l-d1)/2/(2.0*g);


z1 = 0.0;

// Bottom side

Point(1)  = {-l/2    , ymin  , z1}; 
Point(2)  = {-d1/2   , ymin  , z1}; 
Point(3)  = {d1/2    , ymin  , z1}; 
Point(4)  = {l/2     , ymin  , z1}; 

// Top side

Point(5)  = {l/2      , ymax   , z1}; 
Point(6)  = {d1/2     , ymax   , z1}; 
Point(7)  = {-d1/2    , ymax   , z1}; 
Point(8)  = {-l/2     , ymax   , z1};


// Lines

Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 5};
Line(5) = {5, 6};
Line(6) = {6, 7};
Line(7) = {7, 8};
Line(8) = {8, 1};

// Auxiliary lines

Line(9) = {2, 7};
Line(10) = {3, 6};

// Line loops

Line Loop(1) = {7,8,1,9};
Line Loop(2) = {6, -9, 2, 10};
Line Loop(3) = {5, -10, 3, 4};

// Surface definition

Plane Surface(51) = {1};
Plane Surface(52) = {2};
Plane Surface(53) = {3};

//  Transfinite lines and surfaces for structured grids
 
Transfinite Line{2,6}       =  nIn;    // Nozzle inlet -> Reactants
Transfinite Line{4,8,9,10}  =  nGap;    // Gap
Transfinite Line{1,3,5,7}   =  nLen;    // Walls

Transfinite Surface{51};
Transfinite Surface{52};
Transfinite Surface{53};

// Recombination to get QUADS instead of TRIA

Recombine Surface {51};
Recombine Surface {52};
Recombine Surface {53};

// Boundary conditions

Physical Line("Fuel") = {2};
Physical Line("Oxidizer") = {6};
Physical Line("Wall") = {1, 3, 5, 7};
Physical Line("Outlet") = {4,8};
Physical Surface("Surf") = {51,52,53};

