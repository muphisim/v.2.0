// Turn on labels
SetFactory("OpenCASCADE");

//Inputs
gridsize = 0.1;
smallGrid = 0.1;

// Make the box and flap
Rectangle(1) = {0, 0, 0, 100, 1};

// Fine tune mesh
box[] = PointsOf{Surface{1};};
Characteristic Length{box[]} = gridsize;


Physical Curve(1) = {4};
Physical Curve(2) = {2};
Physical Curve(3) = {3};
Physical Curve(4) = {1};
Physical Surface(1) = {1};

