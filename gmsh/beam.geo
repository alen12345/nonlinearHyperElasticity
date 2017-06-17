// MESH PARAMS
a = 0.25;
lc = a / 3;

// POINT GEOMETRY
Point(1) = {-a ,-a , 0, lc};
Point(2) = {a, -a, 0, lc};
Point(3) = {a, a, 0, lc};
Point(4) = {-a, a, 0, lc};

// LINE GEOMETRY
Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 1};

// LINE LOOPS
Line Loop(1) = {1, 2, 3, 4};

// SURFACES
Plane Surface(1) = {1};

// EXTRUSION
beam[] = Extrude {0, 0, a*40} {Surface{1}; };


// PHYSICAL ENTITIES
Physical Surface("allSurface", 1) = {beam[]};
Physical Surface("clamp", 2) = {1};
Physical Surface("topFace", 3) = {beam[2]};
Physical Surface("rightFace", 4) = {beam[0]};

Physical Volume("allVolume", 1) = {beam[1]};
