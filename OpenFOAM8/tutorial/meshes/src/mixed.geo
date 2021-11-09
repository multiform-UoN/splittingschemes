
l1 = 10.0;
dl = 0.2;
l2 = 1.8;
h = 1.0;
hbl = 0.15;


Point(1)  = {       0, -0.5,     0, 1.0};
Point(2)  = {l1      , -0.5,     0, 1.0};
Point(3)  = {l1+dl   , -0.5,     0, 1.0};
Point(4)  = {l1+dl+l2, -0.5,     0, 1.0};
Point(5)  = {       0, -0.5,     h, 1.0};
Point(6)  = {l1      , -0.5,     h, 1.0};
Point(7)  = {l1+dl   , -0.5,     h, 1.0};
Point(8)  = {l1+dl+l2, -0.5,     h, 1.0};
Point(9)  = {l1+dl   , -0.5,   hbl, 1.0};
Point(10) = {l1+dl   , -0.5, h-hbl, 1.0};
Point(11) = {l1+dl+l2, -0.5,   hbl, 1.0};
Point(12) = {l1+dl+l2, -0.5, h-hbl, 1.0};
//+
Line(1) = {2, 1};
Line(2) = {6, 5};
Line(3) = {7, 8};
Line(4) = {10, 12};
Line(5) = {9, 11};
Line(6) = {3, 4};
Line(7) = {7, 6};
Line(8) = {3, 2};
Line(9) = {1, 5};
Line(10) = {2, 6};
Line(11) = {3, 9};
Line(12) = {7, 10};
Line(13) = {10, 9};
Line(14) = {8, 12};
Line(15) = {4, 11};
Line(16) = {11, 12};
//+
Curve Loop(1) = {1, 9, -2, -10}; Plane Surface(1) = {1};
Curve Loop(2) = {8, 10, -7, 12, 13, -11}; Plane Surface(2) = {2};
Curve Loop(3) = {4, -14, -3, 12}; Plane Surface(3) = {3};
Curve Loop(4) = {6, 15, -5, -11}; Plane Surface(4) = {4};
Curve Loop(5) = {5, 16, -4, 13}; Plane Surface(5) = {5};

nbl  = 60;
nzcnt = 100;
nxl1 = 200;
nzl1 = 50;
nxdl = 50;
ratio = 2.0;
//+
Transfinite Curve {15, 11, 12, 14} = nbl Using Progression 1.01;
Transfinite Curve {3, 4, 5, 6} = ((l2*nbl)/(hbl*ratio)) Using Progression 1;
Transfinite Curve {16, 13} = nzcnt Using Bump 0.3;
Transfinite Curve {2, 1} = nxl1 Using Progression 1;
Transfinite Curve {10, 9} = nzl1 Using Bump 0.2;
Transfinite Curve {7, 8} = nxdl Using Progression 1.025;
//+
Transfinite Surface {1};
Transfinite Surface {3};
Transfinite Surface {5};
Transfinite Surface {4};
Recombine Surface {1, 4, 5, 3};
//+
Extrude {0, 1, 0} {
  Surface{3}; Surface{5}; Surface{4}; Surface{2}; Surface{1}; Layers{1}; Recombine;
}
//+
Physical Surface("frontAndBack") = {38, 60, 82, 114, 136, 1, 2, 4, 5, 25, 3};
//+
Physical Surface("inlet") = {127};
//+
Physical Surface("outlet") = {29, 51, 73};
//+
Physical Surface("topInlet") = {131, 101};
//+
Physical Surface("topOutlet") = {33};
//+
Physical Surface("groundInlet") = {123, 93};
//+
Physical Surface("groundOutlet") = {69};
//+
Physical Volume("internal") = {1, 2, 3, 4, 5};
