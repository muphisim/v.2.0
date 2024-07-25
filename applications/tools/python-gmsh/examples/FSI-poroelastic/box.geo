// Turn on labels
SetFactory("OpenCASCADE");

//Inputs
gridsize = 0.5;

// paper dimensions
//b = .3;
//a = 0.6;
//h = 0.9;
//scale = 4;
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

Point(1) = {l, -(h+b+a), 0, gridsize};
Point(2) = {l, (h+b+a), 0, gridsize};

Line(1) = {1, 2};

//Transfinite Curve {1} = scale*(2*(a+b+h))+5 Using Progression 1;
Transfinite Curve {1} = 29 Using Progression 1;

surfaces[] = Extrude { l, 0,0 }{Line{1}; Layers{4*scale*(a+h+b)}; Recombine;};
//BooleanUnion{Surface{1};Delete;}{Surface{2};Delete;}
//BooleanUnion{Surface{1};Delete;}{Surface{3};Delete;}



//Now make 3D by extrusion.
//newEntities[] = Extrude { 0,0,.1 }{Surface{1}; Layers{1}; Recombine;};


Physical Surface(1) = {1};
Physical Curve(1) = {1};
Physical Curve(3) = {3};
Physical Curve(2) = {2};
Physical Curve(4) = {4};




