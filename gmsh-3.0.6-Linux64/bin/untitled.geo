//+
Point(1) = {0, 0, 0, 1.0};
//+
Point(2) = {0, 1, 0, 1.0};
//+
Point(3) = {0.7, 0.4, -0, 1.0};
//+
Point(4) = {1, 0, -0, 1.0};
//+
Line(1) = {2, 1};
//+
Line(2) = {1, 4};
//+
Line(3) = {2, 4};
//+
Recursive Delete {
  Line{3}; Point{3}; 
}
//+
Point(5) = {1, 1, -0, 1.0};
//+
Line(3) = {2, 5};
//+
Line(4) = {5, 4};
//+
Line Loop(1) = {1, 2, -4, -3};
//+
Plane Surface(1) = {1};
//+
Physical Surface("e") = {1};
