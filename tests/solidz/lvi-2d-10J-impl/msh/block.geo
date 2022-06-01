/*

block - 2d

Inputs:
L	- Length
T	- Thickness

Mesh:
lc_# 	- Characteristic Length of Edge/Point
nnZ     - Number of Nodes in Z or width direction

*/

// - - - - - - - - - - - - - - - - - - - - - - -
// Inputs
// - - - - - - - - - - - - - - - - - - - - - - -

// Geometry
a    =  8.0;  // mm
T    =  1.0;  // mm

// Mesh
lc_S   =  0.75;  // element size

// - - - - - - - - - - - - - - - - - - - - - - -
// Start of Code
// - - - - - - - - - - - - - - - - - - - - - - -

lc=lc_S;

// Create 4 points to define the base
//                       { x     y  z mesh}
p01 = newp; Point(p01) = { 0,    0, 0,  lc};
p02 = newp; Point(p02) = { a,    0, 0,  lc};
p03 = newp; Point(p03) = { a,    T, 0,  lc};
p04 = newp; Point(p04) = { 0,    T, 0,  lc};

// Define lines
l01 = newl; Line(l01) = {p01,p02};
l02 = newl; Line(l02) = {p02,p03};
l03 = newl; Line(l03) = {p03,p04};
l04 = newl; Line(l04) = {p04,p01};

// Line loop and define surface (base)
ll01=newll; Line Loop(ll01) = {l01:l04};

// Base surface
s01=news; Plane Surface(s01) = {ll01};

// - - - - - - - - - - - - - - - - - - - - - - -
// Material Codes
// - - - - - - - - - - - - - - - - - - - - - - -

Physical Surface("MATER_BULK") = {s01};

// - - - - - - - - - - - - - - - - - - - - - - -
// Boundary Conditions
// - - - - - - - - - - - - - - - - - - - - - - -

Physical Line("BCBOU_BOT")   = {l01};
Physical Line("BCBOU_TOP")   = {l03};

Mesh.CharacteristicLengthFactor = lc_S;

