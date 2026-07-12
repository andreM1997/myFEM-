import numpy as np
from datapreparation import def_data
import meshio
import sympy as sp
import scipy
from sksparse.cholmod import cholesky
from scipy.sparse import csr_matrix, csc_matrix
from scipy.sparse.linalg import spsolve, eigsh




class solveFEM:

    def __init__(self,lineare_FEM):  #am besten in 3 klassen getrennt    meshdata, FEM_aufstellen, GUI einstellungen (z.b welcher solver + welche einstellungen für solver.... oä.) 
        
        self.row             = np.array(lineare_FEM.row)
        self.column          = np.array(lineare_FEM.column)
        self.value           = np.array(lineare_FEM.value)

        self.ForceVector     = np.array(lineare_FEM.globalForceVector)

        self.nodedata        = lineare_FEM
        self.nodecoordinate  = self.nodedata.nodecoordinate
        self.node2dof        = self.nodedata.node2dof
        self.node2cell       = self.nodedata.node2cell


        #steifigkeitsmatrix
        #self.FEMdata = FEMdata    # in FEMdata werden zb forcevector, stiffnessmatrix und später auch massenmatrix dämpfungsmatrix gespeichert
        #self.row     = self.FEMdata.row
        #self.column  = self.FEMdata.column
        #self.value   = self.FEMdata.value

        # hier stehen grundlegende Daten des meshes drinnen (irgendwie noch welche elementtypen oder so das sind?? evtl für visulalisierung)
        #self.nodedata        = nodedata
        #self.nodecoordinate  = self.nodedata.nodecoordinate
        #self.node2dof        = self.nodedata.node2dof
        #self.node2cell       = self.nodedata.node2cell

        #self.solverdata           = solverdata

        

    def solve_static_structural(self):


        dim_K=int(np.shape(self.node2dof)[0]*np.shape(self.node2dof)[1])
        K_csr = csr_matrix(((self.value, self.column, self.row)), shape=(dim_K, dim_K))
        K_csc = K_csr.tocsc()





        # === DEBUG: globale Steifigkeitsmatrix nach Penalty ===
        K_test = csr_matrix((self.value, self.column, self.row))

        ndof = K_test.shape[0]

        K_csr = K_csr.astype(np.float64)
        self.ForceVector = self.ForceVector.astype(np.float64).reshape(-1,1)

        displacements = spsolve(K_csr, self.ForceVector)



        # Cholesky-Zerlegung
        factor = cholesky(K_csc)

        print(np.shape(self.ForceVector))

        # Lösen K u = f
        displacements = factor(self.ForceVector)

        #print("shape des lösungsvektors",np.shape(u))
        return displacements

        
        


    
        
