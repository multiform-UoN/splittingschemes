L = 6.0;
H = 0.8;
W = 3;
h = 0.02;

Point(1)  = {0,         0,   0, 1.0};
Point(2)  = {0.5*(L+H), 0,   0, 1.0};
Point(3)  = {0.5*(L-H), 0,   0, 1.0};
Point(4)  = {L,         0,   0, 1.0};
Point(5)  = {L,         0,   H, 1.0};
Point(6)  = {0.5*(L+H), 0,   H, 1.0};
Point(7)  = {0.5*(L+H), 0,   W, 1.0};
Point(8)  = {0.5*(L-H), 0,   W, 1.0};
Point(9)  = {0.5*(L-H), 0,   H, 1.0};
Point(10) = {0,         0,   H, 1.0};
Point(11) = {0.5*(L+H), 0, W-H, 1.0};
Point(12) = {0.5*(L-H), 0, W-H, 1.0};
Point(13) = {L,         0, W-H, 1.0};
Point(14) = {L,         0,   W, 1.0};

Line(1) = {1, 3};
Line(2) = {3, 2};
Line(3) = {2, 4};
Line(4) = {10, 9};
Line(5) = {9, 6};
Line(6) = {6, 5};
Line(7) = {3, 9};
Line(8) = {9, 12};
Line(9) = {12, 8};
Line(10) = {8, 7};
Line(11) = {7, 14};
Line(12) = {2, 6};
Line(13) = {6, 11};
Line(14) = {11, 13};
Line(15) = {11, 12};
Line(16) = {11, 7};
Line(17) = {13, 14};
Line(18) = {4, 5};
Line(19) = {1, 10};

Curve Loop(1) = {18, -6, -12, 3};
Plane Surface(1) = {1};
Curve Loop(2) = {12, -5, -7, 2};
Plane Surface(2) = {2};
Curve Loop(3) = {7, -4, -19, 1};
Plane Surface(3) = {3};
Curve Loop(4) = {13, 15, -8, 5};
Plane Surface(4) = {4};
Curve Loop(5) = {9, 10, -16, 15};
Plane Surface(5) = {5};
Curve Loop(6) = {11, -17, -14, 16};
Plane Surface(6) = {6};

Transfinite Curve {11, 14, 6, 3, 4, 1} = (0.5*(L-H)/h)+1 Using Progression 1;
Transfinite Curve {13, 8} = ((W-2*H)/h)+1 Using Progression 1;
Transfinite Curve {17, 16, 10, 9, 15, 18, 5, 12, 2, 7, 19} = (H/h)+1 Using Progression 1;

Recombine Surface {6, 5, 4, 1, 2, 3};
Transfinite Surface {6};
Transfinite Surface {5};
Transfinite Surface {4};
Transfinite Surface {1};
Transfinite Surface {2};
Transfinite Surface {3};
Extrude {0, -0.2, 0} {
  Surface{6}; Surface{5}; Surface{4}; Surface{1}; Surface{2}; Surface{3}; Layers{1}; Recombine;
}
Physical Surface("anode") = {146};
Physical Surface("cathode2") = {32};
Physical Surface("cathode1") = {94};
Physical Surface("frontAndBack") = {6, 41, 5, 63, 4, 85, 2, 129, 1, 107, 3, 151};
Physical Surface("wall") = {28, 54, 50, 80, 142, 72, 98, 36, 106, 128, 150};
Physical Volume("internal") = {2, 1, 3, 5, 4, 6};
