/*

Piezoelectric - 3d

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
b    =  8.0;  // mm
T    =  1.0;  // mm

// Mesh
lc_S   =  0.75;      // element size
nnL  = Ceil(a/lc_S); // Number of Nodes lines
nnT  = Ceil(T/lc_S); // Number of Nodes lines
tetra = 1;
structured =1;

// - - - - - - - - - - - - - - - - - - - - - - -
// Start of Code
// - - - - - - - - - - - - - - - - - - - - - - -
/* numv[]      contains in the following order:
   numv[0]   = surfaceNumber of opposing Plane
   numv[1]   = volumeNumber of extruded Volume
   numv[...] = surfaceNumber of surrounding surfaces in the same order as in Line Loop */

lc=lc_S;

// Create 4 points to define the base
//                       { x     y  z mesh}
p01 = newp; Point(p01) = { 0,    0, 0,  lc};
p02 = newp; Point(p02) = { a,    0, 0,  lc};
p03 = newp; Point(p03) = { a,    b, 0,  lc};
p04 = newp; Point(p04) = { 0,    b, 0,  lc};

// Define lines
l01 = newl; Line(l01) = {p01,p02};
l02 = newl; Line(l02) = {p02,p03};
l03 = newl; Line(l03) = {p03,p04};
l04 = newl; Line(l04) = {p04,p01};

// Line loop and define surface (base)
ll01=newll; Line Loop(ll01) = {l01:l04};

// Base surface
s01=news; Plane Surface(s01) = {ll01};

// Extrude bottom arm
//If (structured == 1)
//  If (tetra == 1)
    numv1[] = Extrude {0,0,T} {Surface{s01}; Layers{nnT};};
//  Else
//    numv1[] = Extrude {0,0,T} {Surface{s01}; Layers{nnT}; Recombine;};
 // EndIf
//Else
//  numv1[] = Extrude {0,0,T} {Surface{s01}; Layers{nnT};};
//EndIf
// - - - - - - - - - - - - - - - - - - - - - - -
// Material Codes
// - - - - - - - - - - - - - - - - - - - - - - -

Physical Volume("MATER_PZT") = {numv1[1]};

// - - - - - - - - - - - - - - - - - - - - - - -
// Boundary Conditions
// - - - - - - - - - - - - - - - - - - - - - - -

Physical Surface("BCBOU_BOT")   = {s01};
Physical Surface("BCBOU_TOP")   = {numv1[0]};

//Mesh.CharacteristicLengthFactor = lc_S;

//Transfinite Line {l01,l02,l03,l04} = nnL;
//Transfinite Surface {s01};

//If (structured == 1 && tetra == 0)
//  Recombine Surface "*";
//EndIf