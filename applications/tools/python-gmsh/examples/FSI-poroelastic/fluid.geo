// Turn on labels
SetFactory("OpenCASCADE");

//Inputs
gridsize = 0.5;

// paper dimensions (all in cm)
//b = .3;
//a = 0.06;
//h = 0.09;
//scale = 80;
//l = 1.3;

// square cavity
//b = 2.5;
//a = 0.06;
//h = 0.09;
//scale = 60;
//l = 8;

a = 1;
h = 1.5;
b = 2.5;
l = 10;
scale = 3;
// Now make the box and flap

Point(1) = {0, -(h+b+a), 0, gridsize};
Point(2) = {0, -(b+a), 0, gridsize};
Point(3) = {0, -a, 0, gridsize};
Point(4) = {0, a, 0, gridsize};
Point(5) = {0, (b+a), 0, gridsize};
Point(6) = {0, (h+b+a), 0, gridsize};

Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 5};
Line(5) = {5, 6};
Transfinite Curve {2, 4} = scale*b+1 Using Progression 1;
Transfinite Curve {1, 5} = scale*h+1 Using Progression 1;
Transfinite Curve {3} = scale*2*a+1 Using Progression 1;

surfaces[] = Extrude { l, 0,0 }{Line{1, 2, 3, 4, 5}; Layers{4*scale*(a+h+b)}; Recombine;};
//BooleanUnion{Surface{1};Delete;}{Surface{2};Delete;}
//BooleanUnion{Surface{1};Delete;}{Surface{3};Delete;}



//Now make 3D by extrusion.
newEntities[] = Extrude { 0,0,1 }{Surface{1,2, 3, 4, 5}; Layers{1}; Recombine;};


Physical Surface("frontAndBack") = {1, 2, 3, 4, 5, 26, 22, 18, 14, 10};
Physical Surface("bottom") = {6};
Physical Surface("right") = {25, 21, 17, 13, 9};
Physical Surface("top") = {23};
Physical Surface("left") = {12, 20};
Physical Surface("outlet") = {24, 8};
Physical Surface("inlet") = {16};

Physical Volume("fluid") = {1, 2, 3, 4, 5};


