// ---------------------------
// Balken 1x1x6 mm, fein genug für axiales Locking
// ---------------------------
Lx = 1;
Ly = 1;
Lz = 3;
Nx = 4;  // Anzahl Elemente in x-Richtung
Ny = 4;  // Anzahl Elemente in y-Richtung
Nz = 12; // Anzahl Elemente in z-Richtung

// Punkte
Point(1) = {0,0,0,1};
Point(2) = {Lx,0,0,1};
Point(3) = {Lx,Ly,0,1};
Point(4) = {0,Ly,0,1};
Point(5) = {0,0,Lz,1};
Point(6) = {Lx,0,Lz,1};
Point(7) = {Lx,Ly,Lz,1};
Point(8) = {0,Ly,Lz,1};

// Linien
Line(1) = {1,2}; Line(2) = {2,3}; Line(3) = {3,4}; Line(4) = {4,1};
Line(5) = {5,6}; Line(6) = {6,7}; Line(7) = {7,8}; Line(8) = {8,5};
Line(9) = {1,5}; Line(10) = {2,6}; Line(11) = {3,7}; Line(12) = {4,8};

// Flächen
Line Loop(1) = {1,2,3,4}; Plane Surface(1) = {1};
Line Loop(2) = {5,6,7,8}; Plane Surface(2) = {2};
Line Loop(3) = {1,10,-5,-9}; Plane Surface(3) = {3};
Line Loop(4) = {2,11,-6,-10}; Plane Surface(4) = {4};
Line Loop(5) = {3,12,-7,-11}; Plane Surface(5) = {5};
Line Loop(6) = {4,9,-8,-12}; Plane Surface(6) = {6};

// Volume
Surface Loop(1) = {1,2,3,4,5,6};
Volume(1) = {1};

// ---------------------------
// Transfinite Lines
Transfinite Line {1,2,3,4} = Nx;  // Breite x-Richtung
Transfinite Line {5,6,7,8} = Nx;
Transfinite Line {9,10,11,12} = Nz; // Länge z-Richtung

// Querschnitt y-Richtung (Höhe)
Transfinite Line {1,2,3,4,5,6,7,8} = Ny;

// Transfinite Surfaces & Volume
Transfinite Surface {1,2,3,4,5,6};
Transfinite Volume {1};

// ---------------------------
// Kein Recombine
// Recombine Volume {1}; // NICHT verwenden

// Mesh 3D
Mesh 3;