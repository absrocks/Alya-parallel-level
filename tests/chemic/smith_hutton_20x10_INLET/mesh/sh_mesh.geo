// -------------------------------------------------------------------------
// Smith-Hitton problem
//   R.M.Smith and A.G.Hutton (1982)
//
//    ^ |- - - - - - - - - - - - |
//    | |                        |
//    | |                        |
//  L | |                        |
//    | |          (0,0)         |
//    | |- - - - - - - - - - - - |
//    ^              |
//             L     |     L
//       <---------->|<---------->
//                  
//   Daniel Mira
//   28/06/2017
// -------------------------------------------------------------------------

L  = 1.0;      // Characteristic length channel L

z1 =  0.0;
c  =  0.0;

R1 = 4.0;        // Ratio length to jet
R2 = 25.0;        // Ratio length from jet to end
R3 = 2.0;        // Length pipe

// Top plane
Point(1)  = { c    , c    , z1};   // Domain center
Point(2)  = {-L + c, c    , z1};
Point(3)  = {-L + c, c + L, z1};
Point(4)  = { L + c, c + L, z1};
Point(5)  = { L + c, c    , z1};

// Axiliary points
Point(6)  = { c    , c + L, z1};


// Inlet
Line(1) = {1, 2};

// Left boundary
Line(2) = {2, 3};

// Top boundary (left)
Line(3) = {3, 6};

// Top boundary (right)
Line(6) = {6, 4};

// Right boundary
Line(4) = {4, 5};

// Outlet
Line(5) = {5, 1};

// Symmetry plane
Line(7) = {6, 1};


// Line loops (better use GUI to define it)

Line Loop(8) = {2, 3, 7, 1};
Plane Surface(9) = {8};

Line Loop(10) = {7, -5, -4, -6};
Plane Surface(11) = {10};

//
// Transfinite lines and surfaces for structured grids

// Recommendation: define length scale lc

 lc = 0.1; 
 
 Transfinite Line{1,2,3,4,5,6,7}  =  L/lc ;

// 
 Transfinite Surface{9};
 Transfinite Surface{11};


// Recombination to get QUADS instead of TRIA

 Recombine Surface {9};
 Recombine Surface {11};
 Mesh.Smoothing = 1;

// Boundary conditions

 Physical Line("INFLOW") = {1};
 Physical Line("OUTFLOW") = {5};
 Physical Line("LEFT") = {2};
 Physical Line("RIGHT") = {4};
 Physical Line("TOPL") = {3};
 Physical Line("TOPR") = {6};
 Physical Surface("Surf") = {9,11};


