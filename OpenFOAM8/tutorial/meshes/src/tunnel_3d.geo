
Point(1) = {0, 0,-0.5, 1.0};
Point(2) = {0,-0, 0.5, 1.0};
Point(3) = {0, 1, 0.5, 1.0};
Point(4) = {0, 1,-0.5, 1.0};
//+
Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 1};
//+
Curve Loop(1) = {4, 1, 2, 3};
Plane Surface(1) = {1};
//+
Transfinite Curve {1, 2, 3, 4} = 41 Using Bump 0.2;
//+
Recombine Surface {1};
//+
Transfinite Surface {1};
//+
Extrude {16, 0, 0} {
  Surface{1}; Layers{500}; Recombine;
}
//+
Extrude {2.5, 0, 0} {
  Surface{26}; Layers{130}; Recombine;
}
//+
Extrude {1.5, 0, 0} {
  Surface{48}; Layers{80}; Recombine;
}
//+
Physical Surface("inlet") = {1};
//+
Physical Surface("outlet") = {70};
//+
Physical Surface("top") = {25, 47, 69};
//+
Physical Surface("ground_inlet") = {17};
//+
Physical Surface("ground_sand_bed") = {39};
//+
Physical Surface("ground_outlet") = {61};
//+
Physical Surface("front") = {21, 43, 65};
//+
Physical Surface("back") = {13, 35, 57};
//+
Physical Volume("internal") = {1, 2, 3};
