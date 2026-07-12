// Gmsh project created on Tue Nov 18 22:52:35 2025
SetFactory("OpenCASCADE");

Merge "teststab.stp";

// Mesh densities
nx = 40;
ny = 2;
nz = 2;

// Structured mesh generation
Transfinite Line "*", nx+1;
Transfinite Surface "*";
Transfinite Volume "*";

Recombine Surface "*";
Recombine Volume "*";

Mesh.Algorithm3D = 8;   // Transfinite hex
Mesh 3;

