w = 1;
h = 1.5;

Point(1)  = {  -h/2, -w/2,   0, 1.0};
Point(2)  = {   h/2, -w/2,   0, 1.0};
Point(3)  = {   h/2, -w/2,   h, 1.0};
Point(4)  = {  -h/2, -w/2,   h, 1.0};
Point(5)  = {-0.5, -w/2,   0, 1.0};
Point(6)  = {  -0, -w/2, 0.5, 1.0};
Point(7)  = { 0.5, -w/2, 0.5, 1.0};
Point(8)  = { 0.5, -w/2, 0.4, 1.0};
Point(9)  = { 0.2, -w/2, 0.4, 1.0};
Point(10) = { 0.2, -w/2, 0.2, 1.0};
Point(11) = { 0.5, -w/2,   0, 1.0};
//+
Line(1) = {5, 1};
Line(2) = {5, 6};
Line(3) = {6, 7};
Line(4) = {7, 8};
Line(5) = {8, 9};
Line(6) = {9, 10};
Line(7) = {10, 11};
Line(8) = {11, 2};
Line(9) = {1, 4};
Line(10) = {2, 3};
Line(11) = {4, 3};
//+
Curve Loop(1) = {10, -11, -9, -1, 2, 3, 4, 5, 6, 7, 8};
Plane Surface(1) = {1};
//+
dl = 0.03;
Transfinite Curve {9, 11, 10} = (h/dl)+1 Using Progression 1;
Transfinite Curve {3} = (0.5/dl)+1 Using Progression 1;
Transfinite Curve {8, 1} = (((h/2)-0.5)/dl)+1 Using Progression 1;
Transfinite Curve {2} = (0.707/dl)+1 Using Progression 1;
Transfinite Curve {7} = (0.361/dl)+1 Using Progression 1;
Transfinite Curve {6} = (0.2/dl)+1 Using Progression 1;
Transfinite Curve {5} = (0.3/dl)+1 Using Progression 1;
Transfinite Curve {4} = (0.1/dl)+1 Using Progression 1;
//+
Extrude {0, w, 0} {
  Surface{1}; Layers{1}; Recombine;
}
//+
Physical Surface("back") = {68};
Physical Surface("front") = {1};
Physical Surface("top") = {31};
Physical Surface("left") = {35};
Physical Surface("right") = {27};
Physical Surface("ground") = {67, 39};
Physical Surface("wall") = {63, 59, 55, 51, 47, 43};
Physical Volume("internal") = {1};
