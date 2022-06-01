// -------------------------------------------------------------------------
// Counterflow diffusion flame case
// -------------------------------------------------------------------------

ymin = -2e-3;
ymax = 5e-3;

g  = ymax-ymin;      // Gap distance
d1 = g*1.0;      // Streams diameter
l  = 5.0*g;      // Domain length

nGap = 26;
nIn  = nGap*d1/(1.5*g);
nLen = nGap*(l-d1)/2/(2.0*g);


z1 = 0.0;

// Bottom side

Point(1)  = {-l/2    , ymin  , z1}; 
Point(2)  = {-d1/2   , ymin  , z1}; 
Point(3)  = {d1/2    , ymin  , z1}; 

// Top side

Point(4)  = {d1/2     , ymax  , z1}; 
Point(5)  = {-d1/2    , ymax  , z1}; 
Point(6)  = {-l/2     , ymax  , z1};


// Lines

Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 5};
Line(5) = {5, 6};
Line(6) = {6, 1};

// Auxiliary lines

Line(7) = {2, 5};

// Line loops

Line Loop(1) = {1,7,5,6};
Line Loop(2) = {2,3,4,-7};

// Surface definition

Plane Surface(51) = {1};
Plane Surface(52) = {2};

//  Transfinite lines and surfaces for structured grids
 
Transfinite Line{2,4}       =  nIn;     // Nozzle inlet -> Reactants
Transfinite Line{3,6,7}     =  nGap;    // Gap
Transfinite Line{1,5}       =  nLen;    // Walls

Transfinite Surface{51};
Transfinite Surface{52};

// Recombination to get QUADS instead of TRIA

Recombine Surface {51};
Recombine Surface {52};

// Boundary conditions

Physical Line("Fuel") = {2};
Physical Line("Oxidizer") = {4};
Physical Line("Wall") = {1, 5};
Physical Line("Sepatator") = {3};
Physical Line("Outlet") = {6};
Physical Surface("Surf") = {51,52};

