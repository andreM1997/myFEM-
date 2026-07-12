import pyvista as pv
import numpy as np
import meshio



meshfile = "C:/Users/andre/OneDrive/Desktop/myFEM/validierungsbalken.msh"

#with open(meshfile, "rb") as f:
#    print(f.read(3))


#exit()

a = meshio.read(meshfile)
print(a)

























#COO_stiffness_matrix = [[],[],[]]
COO_stiffness_matrix = []


#COO_stiffness_matrix = np.empty((3,0))

print(np.shape(COO_stiffness_matrix))

a=1
b=2
c=3

COO_stiffness_matrix.append([a,b,c])
print(np.shape(COO_stiffness_matrix))

COO_stiffness_matrix.append([b,b,c*b])

print(COO_stiffness_matrix)
print(np.shape(COO_stiffness_matrix))

print(COO_stiffness_matrix[0][2])



#np.append(COO_stiffness_matrix,[[a],[b],[c]])

print(COO_stiffness_matrix)





exit()







mesh = meshio.read("C:/Users/andre/OneDrive/Desktop/myFEM/teststab2.msh")

node2cells = np.asarray([[0,3,4,5,5,2,7,9],[9,2,5,3,9,9,9,9]])

lord = np.asarray([[3,3,3,3,3,3,3,3,3,3,3],[3,3,3,3,3,3,3,3,3,3,3]])

print(mesh.cells)
print(mesh.cells[26].data)


exit()

for i in range (0,len(mesh.cells)): 
    if mesh.cells[i].type =='hexahedron':
       node2cells = np.append(node2cells,mesh.cells[i].data, axis=0)


if np.shape(lord)[1]> np.shape(node2cells)[1]:
    print("lord [1] existiert")
    node2cells =  np.hstack((node2cells,np.zeros((np.shape(node2cells)[0],np.shape(lord)[1]-np.shape(node2cells)[1]))))

elif np.shape(lord)[1]==np.shape(node2cells)[1]:
    print("blabla")




node2cells = np.append(node2cells,lord,axis=0)
print(node2cells)