# myFEM-
This is a private project where I combined my skills in Python, structural mechanics and numerical analysis to develope a 3D finite element programm for structural analysis based on geometries meshed in gmsh in gmsh 4.1. dataformat.

The data in environment_fem.yml includes every necessary package as well as there versions to insure compatability throughout all libaries. Download this as a requirement to be able to run the code myFEM correctly.
myFEM is only supporting following settings so far:

        mesh:
          - hex8 meshed geometries with gmsh 4.1 data format

        boundary conditions:
          - force boundary conditions
          - fixed supports 

        parts:
          - only single parts, no assemblies supportes yet

        post processing:
          - only displacements in x,y,z direction can be shown so far


        Solver:
          - cholesky decomposition


Next updates will include increased number of supported boundary conditions as well as element types. Furthermore the post processing part will have expanded to calculate distortions and stresses as wells.
        
