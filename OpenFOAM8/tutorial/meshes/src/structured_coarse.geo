// this is a gmsh file
L = 20;
H = 1.0;
X_bed = 10.0;
L_bed = 4.0;
O_bed = L - (X_bed+L_bed);
h_blayer = 0.2;

//+
Point(1)  = {0,           -0.1, 0, 1.0};
Point(2)  = {X_bed,       -0.1, 0, 1.0};
Point(3)  = {X_bed+L_bed, -0.1, 0, 1.0};
Point(4)  = {L,           -0.1, 0, 1.0};
Point(5)  = {0,           -0.1, h_blayer, 1.0};
Point(6)  = {X_bed,       -0.1, h_blayer, 1.0};
Point(7)  = {X_bed+L_bed, -0.1, h_blayer, 1.0};
Point(8)  = {L,           -0.1, h_blayer, 1.0};
Point(9)  = {0,           -0.1, H-h_blayer, 1.0};
Point(10) = {X_bed,       -0.1, H-h_blayer, 1.0};
Point(11) = {X_bed+L_bed, -0.1, H-h_blayer, 1.0};
Point(12) = {L,           -0.1, H-h_blayer, 1.0};
Point(13) = {0,           -0.1, H, 1.0};
Point(14) = {L,           -0.1, H, 1.0};
Point(15) = {X_bed,       -0.1, H, 1.0};
Point(26) = {X_bed+L_bed, -0.1, H, 1.0};

//+
Line(1) = {2, 1};
Line(2) = {6, 5};
Line(3) = {10, 9};
Line(4) = {15, 13};
Line(5) = {2, 3};
Line(6) = {6, 7};
Line(7) = {10, 11};
Line(8) = {15, 26};
Line(9) = {3, 4};
Line(10) = {7, 8};
Line(11) = {11, 12};
Line(12) = {26, 14};
//+
Line(13) = {1, 5};
Line(14) = {2, 6};
Line(15) = {3, 7};
Line(16) = {4, 8};
//+
Line(17) = {13, 9};
Line(18) = {15, 10};
Line(19) = {26, 11};
Line(20) = {14, 12};
//+
Line(21) = {9, 5};
Line(22) = {10, 6};
Line(23) = {11, 7};
Line(24) = {12, 8};
//+
Curve Loop(1) = {1, 13, -2, -14};   Plane Surface(1) = {1};
Curve Loop(2) = {2, -21, -3, 22};   Plane Surface(2) = {2};
Curve Loop(3) = {4, 17, -3, -18};   Plane Surface(3) = {3};
Curve Loop(4) = {7, -19, -8, 18};   Plane Surface(4) = {4};
Curve Loop(5) = {5, 15, -6, -14};   Plane Surface(5) = {5};
Curve Loop(6) = {6, -23, -7, 22};   Plane Surface(6) = {6};
Curve Loop(7) = {9, 16, -10, -15};  Plane Surface(7) = {7};
Curve Loop(8) = {10, -24, -11, 23}; Plane Surface(8) = {8};
Curve Loop(9) = {11, -20, -12, 19}; Plane Surface(9) = {9};
//+
//+
//+

nz_blayer = 60;
nz_cnt = 60;

ratio_inlet  = 16.0;
ratio_cnt    = 16.0;
ratio_outlet = 16.0;

bumpFactor = 0.2;
progressionFactor_blayer = 1.0;

//+
Transfinite Curve {13, 14, 15, 16} = nz_blayer Using Progression progressionFactor_blayer;
Transfinite Curve {17, 18, 19, 20} = nz_blayer Using Progression progressionFactor_blayer;
Transfinite Curve {21, 22, 23, 24} = nz_cnt Using Bump bumpFactor;
//+
Transfinite Curve {1, 2, 3, 4}    = (X_bed/h_blayer)*(nz_blayer/ratio_inlet)Using Progression 1;
Transfinite Curve {5, 6, 7, 8}    = (L_bed/h_blayer)*(nz_blayer/ratio_cnt)Using Progression 1;
Transfinite Curve {9, 10, 11, 12} = (O_bed/h_blayer)*(nz_blayer/ratio_outlet)Using Progression 1;
//+
Transfinite Surface {1};
Transfinite Surface {2};
Transfinite Surface {3};
Transfinite Surface {5};
Transfinite Surface {6};
Transfinite Surface {4};
Transfinite Surface {7};
Transfinite Surface {8};
Transfinite Surface {9};
Recombine Surface {1, 2, 3, 6, 5, 4, 7, 8, 9};
//+
Extrude {0, 0.2, 0} {
  Surface{1}; Surface{2}; Surface{3}; Surface{4}; Surface{6}; Surface{5}; Surface{7}; Surface{8}; Surface{9}; Layers{1}; Recombine;
}
//+
Physical Surface("inlet") = {81, 59, 37};
Physical Surface("outlet") = {213, 191, 169};
Physical Surface("top") = {217, 107, 77};
Physical Surface("ground_inlet") = {33};
Physical Surface("ground_sand_bed") = {143};
Physical Surface("ground_outlet") = {165};
Physical Surface("frontAndBack") = {222, 200, 178, 134, 112, 156, 90, 68, 46, 3, 2, 1, 5, 6, 4, 9, 8, 7};
Physical Volume("internal") = {7, 8, 9, 4, 5, 6, 3, 2, 1};
//+
