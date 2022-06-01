// -------------------------------------------------------------------------
// Minimal flame case
// -------------------------------------------------------------------------
// POINTS                   23            33           43  
//            13o------------o-------------o------------o              
//             /|           /|            /|           /|               
//            / |        22/ |         32/ |        42/ |              
//         12o------------o-------------o------------o  |              
//           |  |         |  |          |  |         |  |              
//           |  |         |  |24        |  |34       |  |44            
//           |14o---------|--o----------|--o---------|--o              
//           | /          | /           | /          | /               
//           |/           |/            |/           |/                
//         11o------------o-------------o------------o                 
//                       21            31           41             
//                                                         
// LINES                107                 207                 307        
//              o   ------------   o   -------------   o   ------------   o              
//          102/|              202/|               302/|              402/|               
//            / |     106        / |     206         / |     306        / |              
//           o   ------------   o   -------------   o   ------------   o  |              
//           |  |103            |  |203             |  |303            |  |403           
//           |  |      108      |  |       208      |  |       308     |  |              
//        101|  o   ---------   |--o   ----------   |--o   ---------   |--o              
//           | /104          201| /204           301| /304          401| /404            
//           |/                 |/                  |/                 |/                
//           o   ------------   o   -------------   o   ------------   o                 
//                  105                 205                 305                  
//                                                         
//                                                         
//                                                         
//                                                        


gx  = 0.33;        
gy  = 1.0;        
gz  = 1.0;       

nGapx = 4; // DAMN  COARSE
nGapy = 4;
nGapz = 4;

z1 = 0.0;


// Points
Point(11)  = {   0.0  , -gy/2  , -gz/2, z1}; 
Point(12)  = {   0.0  , -gy/2  ,  gz/2, z1}; 
Point(13)  = {   0.0  ,  gy/2  ,  gz/2, z1}; 
Point(14)  = {   0.0  ,  gy/2  , -gz/2, z1}; 

Point(21)  = {   gx   , -gy/2  , -gz/2, z1}; 
Point(22)  = {   gx   , -gy/2  ,  gz/2, z1}; 
Point(23)  = {   gx   ,  gy/2  ,  gz/2, z1}; 
Point(24)  = {   gx   ,  gy/2  , -gz/2, z1}; 

Point(31)  = { 2.0*gx , -gy/2  , -gz/2, z1}; 
Point(32)  = { 2.0*gx , -gy/2  ,  gz/2, z1}; 
Point(33)  = { 2.0*gx ,  gy/2  ,  gz/2, z1}; 
Point(34)  = { 2.0*gx ,  gy/2  , -gz/2, z1}; 

Point(41)  = { 3.0*gx , -gy/2  , -gz/2, z1}; 
Point(42)  = { 3.0*gx , -gy/2  ,  gz/2, z1}; 
Point(43)  = { 3.0*gx ,  gy/2  ,  gz/2, z1}; 
Point(44)  = { 3.0*gx ,  gy/2  , -gz/2, z1}; 


// Lines
Line(101) = {11, 12};
Line(102) = {12, 13};
Line(103) = {13, 14};
Line(104) = {14, 11};
Line(105) = {11, 21};
Line(106) = {12, 22};
Line(107) = {13, 23};
Line(108) = {14, 24};

Line(201) = {21, 22};
Line(202) = {22, 23};
Line(203) = {23, 24};
Line(204) = {24, 21};
Line(205) = {21, 31};
Line(206) = {22, 32};
Line(207) = {23, 33};
Line(208) = {24, 34};

Line(301) = {31, 32};
Line(302) = {32, 33};
Line(303) = {33, 34};
Line(304) = {34, 31};
Line(305) = {31, 41};
Line(306) = {32, 42};
Line(307) = {33, 43};
Line(308) = {34, 44};

Line(401) = {41, 42};
Line(402) = {42, 43};
Line(403) = {43, 44};
Line(404) = {44, 41};

// Line loops
Line Loop(611) = { 101, 102, 103, 104};
Line Loop(612) = {-101, 105, 201,-106};
Line Loop(613) = {-102, 106, 202,-107};
Line Loop(614) = {-103, 107, 203,-108};
Line Loop(615) = {-104, 108, 204,-105};

Line Loop(621) = { 201, 202, 203, 204};
Line Loop(622) = {-201, 205, 301,-206};
Line Loop(623) = {-202, 206, 302,-207};
Line Loop(624) = {-203, 207, 303,-208};
Line Loop(625) = {-204, 208, 304,-205};

Line Loop(631) = { 301, 302, 303, 304};
Line Loop(632) = {-301, 305, 401,-306};
Line Loop(633) = {-302, 306, 402,-307};
Line Loop(634) = {-303, 307, 403,-308};
Line Loop(635) = {-304, 308, 404,-305};

Line Loop(641) = {-404,-403,-402,-401};

// Surface definition
Plane Surface(701) = {611};
Plane Surface(702) = {612};
Plane Surface(703) = {613};
Plane Surface(704) = {614};
Plane Surface(705) = {615};
Plane Surface(706) = {621};

Plane Surface(802) = {622};
Plane Surface(803) = {623};
Plane Surface(804) = {624};
Plane Surface(805) = {625};

Plane Surface(901) = {631};
Plane Surface(902) = {632};
Plane Surface(903) = {633};
Plane Surface(904) = {634};
Plane Surface(905) = {635};
Plane Surface(906) = {641};

// Surface loop
Surface Loop(666) = {701,702,703,704,705,706};
Surface Loop(667) = {706,802,803,804,805,901};
Surface Loop(668) = {901,902,903,904,905,906};

// Volume
Volume(1) = {666};
Volume(2) = {667};
Volume(3) = {668};


//  Transfinite lines and surfaces for structured grids 
Transfinite Line{105,106,107,108,205,206,207,208,305,306,307,308}     =  nGapx;
Transfinite Line{102,104,202,204,302,304,402,404}                     =  nGapy;
Transfinite Line{101,103,201,203,301,303,401,403}                     =  nGapz;
//Transfinite Surface{701,702,703,704,705,706,802,803,804,805,901,902,903,904,905,906};
Transfinite Surface{701,702,703,704,705,706,901,902,903,904,905,906};
Transfinite Volume{1,3};

// Recombination to get QUADS instead of TRIA
//Recombine Surface{701,702,703,704,705,706,802,803,804,805,901,902,903,904,905,906};
Recombine Surface{701,702,703,704,705,706,901,902,903,904,905,906};
Recombine Volume {1};

// Boundary conditions
Physical Surface("Left") = {701};
Physical Surface("Right")= {901};
Physical Surface("Wall") = {702,703,704,705,802,803,804,805,902,903,904,905};
Physical Volume("Fluid") = {1,2,3};

