// Turn on labels
SetFactory("OpenCASCADE");
Mesh.Algorithm = 8;



// Initialise parameters - user coded. 
// Mesh parameters:
meshSize = 0.15;
// Box parameters
Lx = 8; 
Ly = 6;
Lz = 30;
flapWidth = 0.2;
flapHeight = 3.5;
flapBase = flapWidth;
flapStart = 0;
flapLength = Lz;
// Catheter parameters
xc = 2;
yc = 2.5;
r = 1.25;
ri = 0.75;
catheterLength = Lz;

// Hole parameters updated in simulation
theta = 10;
rotation = 0;
t[] = {theta, theta, theta, theta}; 
phi[] = {rotation, rotation+90, rotation+180, rotation+270};
offset[] = {0, 90,0, 90,0, 90,0, 90,0, 90};
startDist = {0.0};
endDist = {Lz};

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
flapFull[] = Extrude { 0, 0, flapLength }{Surface{flapS};};
flap = flapFull[1];
// End solid


// Make box 
boxR = newp; Point(boxR) = {flapBase/2, 0, 0, meshSize}; 
p = newp; Point(p) = {Lx/2, Ly/2, 0, meshSize}; 
p = newp; Point(p) = {flapBase/2, Ly, 0, meshSize}; 
p = newp; Point(p) = {-flapBase/2, Ly, 0, meshSize}; 
p = newp; Point(p) = {-Lx/2, Ly/2, 0, meshSize}; 
boxL = newp; Point(boxL) = {-flapBase/2, 0, 0, meshSize}; 
spline = newl; Spline(spline)={boxR:boxL:1}; 
base = newl; Line(base) = {boxL, boxR};
cl = newcl; Curve Loop(cl) = {spline, base};
boxS = news; Plane Surface(boxS) = {cl}; 

boxFull[] = Extrude { 0,0, Lz}{Surface{boxS};};
box = boxFull[1];


//Make catheter 
shuntOuter = newv;
Cylinder(shuntOuter) = {xc, yc, 0, 0, 0, catheterLength, r, 2*Pi};
shuntTip = newv;
Sphere(shuntTip) = {xc, yc, catheterLength, r};
BooleanUnion{Volume{shuntOuter};Delete;}{Volume{shuntTip};Delete;}
shuntInner = newv;
Cylinder(shuntInner) = {xc, yc, 0, 0, 0, catheterLength, ri, 2*Pi};
shuntOuter = box+1;
BooleanDifference{Volume{shuntOuter};Delete;}{Volume{shuntInner};Delete;}


numHolesLength = #offset[];
numHolesAxial = #t[]; 
holesEnd = catheterLength - endDist;
holesLength = endDist - startDist;
//Printf("start = %g, length = %g", holesEnd, holesLength);


// Set up each of the hole wedges    
For i In {0:numHolesAxial-1:1}
    For j In {1:numHolesLength:1}
    	zPos = holesEnd + j*holesLength/(numHolesLength+1);
        xPos =  xc + r * Cos(Pi/180*(phi[i] + offset[j-1])); yPos = yc + r * Sin(Pi/180*(phi[i] + offset[j-1])); rHole = r * Tan(Pi/180* t[i]);
        If (rHole > 1/2 * holesLength / (numHolesLength+1))
        	rHole = 1/2 * holesLength /(numHolesLength+1);
            Printf("rHole modified");
        EndIf
        xAxis = xPos - xc; yAxis = yPos - yc;
        v = newv; 
        If (t[i] < 10)
	        Cylinder(v) = {xPos, yPos, zPos, -xAxis, -yAxis, 0, rHole};
        Else
	        Cone(v) = {xPos, yPos, zPos, -xAxis, -yAxis, 0, rHole, 0};
        EndIf
	    BooleanDifference{Volume{shuntOuter};Delete;}{Volume{v}; Delete;}
    EndFor
EndFor
BooleanDifference{Volume{box}; Delete;}{Volume{shuntOuter}; Delete;}
Coherence;
BooleanDifference{Volume{box};Delete;}{Volume{flap}; Delete;}




// Identify groups of surfaces with phyiscal entities
surfs[]=Surface "*";
Physical Surface("walls") = {};
Physical Surface("flap") = {};
Physical Surface("shunt") = {};
Physical Surface("shuntHoles") = {};
// Identify back and front walls
For i In {0:#surfs[]-1:1}
	s = surfs(i);
	Printf("i %g, s[i] %g", i, s);
	bdyline() = Boundary{Surface{s};};
	k=#bdyline();
	Printf("len bdyline %g", k);
	bdy[] = Boundary{Line{bdyline[0]};};
	coords = Point{bdy[0]};
	If (k > 3) // Selects the back and front walls
		Printf("found k > 3, k= %g", k);
		If ((coords[0]-xc)*(coords[0]-xc) + (coords[1]-yc)*(coords[1]-yc) < r*r + .1)
			Physical Surface("shunt",3) += {s};
		Else
			Physical Surface("walls",1) += {s};
		EndIf
	ElseIf (k == 1)
		Printf("found k ==1 k=%g", k);
		If (coords[2]==0)
		    Physical Surface("shuntHoles", 4) += {s};
		Else
		    Physical Surface("shunt", 3) += {s};
		EndIf
	ElseIf (k == 3) //(coords[0] < xc - r - 0.01)
		Physical Surface("flap",1) += {s};
	EndIf
		
EndFor

Physical Volume("fluid") = {box};
Characteristic Length{Point{:}} = meshSize;
//Characteristic Length{PointsOf{Physical Surface{1};}} = 3*meshSize;

