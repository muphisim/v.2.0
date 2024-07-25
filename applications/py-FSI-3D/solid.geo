// Turn on labels
SetFactory("OpenCASCADE");


// Mesh parameters:
meshSize = 0.5;
// Box parameters
Lx = 8; 
Ly = 6;
Lz = 30;
flapWidth = 0.2;
flapHeight = 3.5;
flapBase = flapWidth;
flapStart = 0;
flapLength = Lz;


// Make flap
baseR = newp; Point(baseR) = {flapBase/2, 0, flapStart, meshSize}; 
p = newp; Point(p) = {3/4*flapWidth, 1/4*flapHeight, flapStart, meshSize}; 
p = newp; Point(p) = {3/4*flapWidth, 3/4*flapHeight, flapStart, meshSize}; 
p = newp; Point(p) = {flapWidth/2, flapHeight, flapStart, meshSize}; 
p = newp; Point(p) = {-flapWidth/2, flapHeight, flapStart, meshSize}; 
p = newp; Point(p) = {-3/4*flapWidth, 3/4*flapHeight, flapStart, meshSize}; 
p = newp; Point(p) = {-3/4*flapWidth, 1/4*flapHeight, flapStart, meshSize}; 
baseL = newp; Point(baseL) = {-flapBase/2, 0, flapStart, meshSize}; 
ls1 = newl; Spline(ls1)={baseR:baseL:1}; 
base = newl; Line(base) = {baseL, baseR};
cl = newcl; Curve Loop(cl) = {ls1, base};
flapS = news; Plane Surface(flapS) = {cl};

If (flapStart + flapLength > Lz)
    flapLength = Lz - flapStart;
EndIf
flapFull[] = Extrude { 0, 0, flapLength }{Surface{flapS};};
flap = flapFull[1];

Physical Surface("moving") = {5};
Physical Surface("slip") = {};
Physical Surface("fixed") = {6};

flapEnd = flapStart + flapLength;
If (flapEnd < Lz && flapStart > 0)
    Physical Surface("moving", 0) += {4, 7}; 
ElseIf (flapEnd < Lz && flapStart == 0)
    Physical Surface("moving", 0) += {7}; 
    Physical Surface("slip", 2) += {4}; 
ElseIf (flapEnd == Lz && flapStart > 0)
    Physical Surface("moving", 0) += {4}; 
    Physical Surface("slip", 2) += {7}; 
Else
    Physical Surface("slip", 2) += {4, 7}; 
EndIf

Physical Volume("solid") = {flap};

Characteristic Length{Point{:}} = meshSize;
