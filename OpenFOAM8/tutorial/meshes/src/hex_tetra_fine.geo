// this is a gmsh file
L = 16.0;
H = 1.0;
L_bed = 6.2;
X_bed = 6.0;
Z1 = 0.2;
Z2 = 0.95;
//+
Point(1)  = {0,           -0.5, 0, 1.0};
Point(2)  = {X_bed,       -0.5, 0, 1.0};
Point(3)  = {X_bed+L_bed, -0.5, 0, 1.0};
Point(4)  = {L,           -0.5, 0, 1.0};

Point(5)  = {0,           -0.5, Z1, 1.0};
Point(6)  = {X_bed,       -0.5, Z1, 1.0};
Point(7)  = {X_bed+L_bed, -0.5, Z1, 1.0};
Point(8)  = {L,           -0.5, Z1, 1.0};

Point(9)  = {0,           -0.5, Z2, 1.0};
Point(10) = {X_bed,       -0.5, Z2, 1.0};
Point(11) = {X_bed+L_bed, -0.5, Z2, 1.0};
Point(12) = {L,           -0.5, Z2, 1.0};

Point(13) = {0,           -0.5, H/2, 1.0};
Point(14) = {L,           -0.5, H/2, 1.0};
Point(15) = {X_bed,       -0.5, H/2, 1.0};
Point(16) = {X_bed+L_bed, -0.5, H/2, 1.0};

Point(17) = {0,           -0.5, H, 1.0};
Point(18) = {L,           -0.5, H, 1.0};
Point(19) = {X_bed,       -0.5, H, 1.0};
Point(20) = {X_bed+L_bed, -0.5, H, 1.0};

Line(1) = {1, 5};
Line(2) = {2, 6};
Line(3) = {3, 7};
Line(4) = {4, 8};
Line(5) = {17, 9};
Line(6) = {19, 10};
Line(7) = {20, 11};
Line(8) = {18, 12};

Line(9) = {5, 13};
Line(10) = {6, 15};
Line(11) = {7, 16};
Line(12) = {8, 14};
Line(13) = {9, 13};
Line(14) = {10, 15};
Line(15) = {11, 16};
Line(16) = {12, 14};

Line(17) = {2, 1};
Line(18) = {6, 5};
Line(19) = {15, 13};
Line(20) = {10, 9};
Line(21) = {19, 17};
Line(22) = {3, 4};
Line(23) = {7, 8};
Line(24) = {16, 14};
Line(25) = {11, 12};
Line(26) = {20, 18};
Line(27) = {2, 3};
Line(28) = {6, 7};
Line(29) = {15, 16};
Line(30) = {10, 11};
Line(31) = {19, 20};

//+
Curve Loop(1) = {17, 1, -18, -2}; Plane Surface(1) = {1};
Curve Loop(2) = {18, 9, -19, -10}; Plane Surface(2) = {2};
Curve Loop(3) = {19, -13, -20, 14}; Plane Surface(3) = {3};
Curve Loop(4) = {21, 5, -20, -6}; Plane Surface(4) = {4};
Curve Loop(5) = {27, 3, -28, -2}; Plane Surface(5) = {5};
Curve Loop(6) = {28, 11, -29, -10}; Plane Surface(6) = {6};
Curve Loop(7) = {29, -15, -30, 14}; Plane Surface(7) = {7};
Curve Loop(8) = {30, -7, -31, 6}; Plane Surface(8) = {8};
Curve Loop(9) = {22, 4, -23, -3}; Plane Surface(9) = {9};
Curve Loop(10) = {23, 12, -24, -11}; Plane Surface(10) = {10};
Curve Loop(11) = {24, -16, -25, 15}; Plane Surface(11) = {11};
Curve Loop(12) = {25, -8, -26, 7}; Plane Surface(12) = {12};

nGround = 100; ratioInGr = 3.0; ratioCntGr = 1.5; ratioOutGr = 2.5;
nTop = 25; ratioInTop = 3.0; ratioCntTop = 3.0; ratioOutTop = 3.0;
nzCntGr = 30; nzCntTop = 30;
nxCntIn = 150; nxCntCnt = 180; nxCntOut = 100;

Transfinite Curve {1, 2, 3, 4} = nGround Using Progression 1.01;
Transfinite Curve {17, 18} = (X_bed*nGround)/(Z1*ratioInGr) Using Progression 1.001;
Transfinite Curve {27, 28} = (L_bed*nGround)/(Z1*ratioCntGr) Using Progression 1;
Transfinite Curve {22, 23} = ((L-(X_bed+L_bed))*nGround)/(Z1*ratioOutGr) Using Progression 1.001;

Transfinite Curve {5, 6, 7, 8} = nTop Using Progression 1;
Transfinite Curve {21, 20} = (X_bed*nTop)/((H-Z2)*ratioInTop) Using Progression 1;
Transfinite Curve {31, 30} = (L_bed*nTop)/((H-Z2)*ratioCntTop) Using Progression 1;
Transfinite Curve {26, 25} = ((L-(X_bed+L_bed))*nTop)/((H-Z2)*ratioOutTop) Using Progression 1;

Transfinite Curve {9, 10, 11, 12} = nzCntGr Using Progression 1.05;
Transfinite Curve {13, 14, 15, 16} = nzCntTop Using Progression 1.08;
Transfinite Curve {19} = nxCntIn Using Progression 1;
Transfinite Curve {29} = nxCntCnt Using Bump 1;
Transfinite Curve {24} = nxCntOut Using Progression 1;

Transfinite Surface {1};
Transfinite Surface {5};
Transfinite Surface {9};
Transfinite Surface {4};
Transfinite Surface {8};
Transfinite Surface {12};
Recombine Surface {1, 5, 9, 4, 8, 12};//+
Extrude {0, 1, 0} {
  Surface{1}; Surface{2}; Surface{3}; Surface{4}; Surface{5}; Surface{6}; Surface{7}; Surface{8}; Surface{9}; Surface{10}; Surface{11}; Surface{12}; Layers{1}; Recombine;
}
//+
Physical Surface("inlet") = {44, 66, 88, 110};
Physical Surface("outlet") = {220, 242, 264, 286};
Physical Surface("top") = {290, 202, 106};
Physical Surface("ground_inlet") = {40};
Physical Surface("ground_sand_bed") = {128};
Physical Surface("ground_outlet") = {216};
Physical Surface("front") = {9, 10, 11, 12, 229, 251, 273, 295, 5, 6, 7, 8};
Physical Surface("back") = {141, 163, 185, 207, 1, 2, 3, 4, 53, 75, 97, 119};
//+
Physical Volume("internal") = {12, 10, 9, 11, 5, 7, 8, 6, 1, 2, 4, 3};
