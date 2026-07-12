import numpy as np
from GUI import FEMGUI
import meshio
import sympy as sp
import inspect
import scipy

from scipy.sparse.linalg import eigsh
from scipy.sparse import csr_matrix, coo_matrix



class def_data:

    xi, eta, zeta = sp.symbols('xi eta zeta')

    # Daten für trilineares Mesh mit hexaeder (formfunktion, ableitung formfunktion, Koordinaten Gaußpunkte, Gewichte Gaußpunkte)
    shape_hex8      =       np.array([sp.lambdify((xi, eta, zeta),  1/8*(1-xi)*(1-eta)*(1-zeta)),  
                                      sp.lambdify((xi, eta, zeta),  1/8*(1+xi)*(1-eta)*(1-zeta)),
                                      sp.lambdify((xi, eta, zeta),  1/8*(1+xi)*(1+eta)*(1-zeta)),
                                      sp.lambdify((xi, eta, zeta),  1/8*(1-xi)*(1+eta)*(1-zeta)),
                                      sp.lambdify((xi, eta, zeta),  1/8*(1-xi)*(1-eta)*(1+zeta)),
                                      sp.lambdify((xi, eta, zeta),  1/8*(1+xi)*(1-eta)*(1+zeta)),
                                      sp.lambdify((xi, eta, zeta),  1/8*(1+xi)*(1+eta)*(1+zeta)),
                                      sp.lambdify((xi, eta, zeta),  1/8*(1-xi)*(1+eta)*(1+zeta))], dtype="object")  # noch korrigieren stimmt nicht mehr mit den ableitungen von der reihenfolge überein




    # ableitungen nach xi, eta , zeta                                                                                                                   #koordinate 3 überarbeiten erste normale shapefunktion aktualisieren und darauf basierend spalte 3 anpassen
    d_shape_hex_8 =         np.array([[sp.lambdify((xi, eta, zeta), -1/8*(1-eta)*(1-zeta)),     sp.lambdify((xi, eta, zeta), -1/8*(1-xi)*(1-zeta)),       sp.lambdify((xi, eta, zeta), -1/8*(1-xi)*(1-eta))],
                                      [sp.lambdify((xi, eta, zeta),  1/8*(1-eta)*(1-zeta)),     sp.lambdify((xi, eta, zeta), -1/8*(1+xi)*(1-zeta)),       sp.lambdify((xi, eta, zeta), -1/8*(1+xi)*(1-eta))],
                                      [sp.lambdify((xi, eta, zeta),  1/8*(1+eta)*(1-zeta)),     sp.lambdify((xi, eta, zeta),  1/8*(1+xi)*(1-zeta)),       sp.lambdify((xi, eta, zeta), -1/8*(1+xi)*(1+eta))],
                                      [sp.lambdify((xi, eta, zeta), -1/8*(1+eta)*(1-zeta)),     sp.lambdify((xi, eta, zeta),  1/8*(1-xi)*(1-zeta)),       sp.lambdify((xi, eta, zeta), -1/8*(1-xi)*(1+eta))],
                                      [sp.lambdify((xi, eta, zeta), -1/8*(1-eta)*(1+zeta)),     sp.lambdify((xi, eta, zeta), -1/8*(1-xi)*(1+zeta)),       sp.lambdify((xi, eta, zeta),  1/8*(1-xi)*(1-eta))],
                                      [sp.lambdify((xi, eta, zeta),  1/8*(1-eta)*(1+zeta)),     sp.lambdify((xi, eta, zeta), -1/8*(1+xi)*(1+zeta)),       sp.lambdify((xi, eta, zeta),  1/8*(1+xi)*(1-eta))],
                                      [sp.lambdify((xi, eta, zeta),  1/8*(1+eta)*(1+zeta)),     sp.lambdify((xi, eta, zeta),  1/8*(1+xi)*(1+zeta)),       sp.lambdify((xi, eta, zeta),  1/8*(1+xi)*(1+eta))],
                                      [sp.lambdify((xi, eta, zeta), -1/8*(1+eta)*(1+zeta)),     sp.lambdify((xi, eta, zeta),  1/8*(1-xi)*(1+zeta)),       sp.lambdify((xi, eta, zeta),  1/8*(1-xi)*(1+eta))]], dtype="object")


    gaus_pos_hex8     =     np.zeros(3)   # gauspunkt bei koordinate 0,0,0 in referenzelement
    
    # Koordinaten Gauspunkte xi, eta, zeta
    gaus_pos_hex8_2   =     np.array([[-0.57735026919, -0.57735026919,  -0.57735026919],
                                      [ 0.57735026919, -0.57735026919,  -0.57735026919],
                                      [ 0.57735026919,  0.57735026919,  -0.57735026919],
                                      [-0.57735026919,  0.57735026919,  -0.57735026919],
                                      [-0.57735026919, -0.57735026919,   0.57735026919],
                                      [ 0.57735026919, -0.57735026919,   0.57735026919],
                                      [ 0.57735026919,  0.57735026919,   0.57735026919],
                                      [-0.57735026919,  0.57735026919,   0.57735026919]])

    # Gewichte Gauspunkte 
    gaus_weight_hex8_2 =     np.array([1,1,1,1,1,1,1,1]) #bereitsmehrdimesnional verrechnet als wi*wj*wk  (alle halt 1)




    def __init__(self, fem_window):

        print("bin in datapreparation angelangt")
        
        self.gui                = fem_window
        self.nodecoordinate     = None
        self.node2cell          = None
        self.node2dof           = None
        self.mat_values         = None
        self.mat_combo          = self.gui.mat_combo
        self.mat_tensor         = None
        self.globalForceVector  = []
        self.jacob_gl           = np.zeros((len(self.gui.mesh_points[:,0])*3,3))
        self.el_type            = self.gui.el_type
        self.el_type_num        = self.gui.el_type_num
        self.n_SF               = self.gui.n_SF
        self.n_elements         = self.gui.n_elements
        self.n_el_node          = self.gui.n_el_node
        self.n_dim              = 3     # 3d FEM
        self.n_GP               = [8]   # für das erste beispiel eigentlich müsste man das aus der Usereingabe erhalten 

        #Stefigkeitsmatrix
        self.COO_stiffness_matrix = []  # row, colum, value
        
        

    
    def nodematrices(self):
        self.nodecoordinate = self.gui.mesh_points                                                                        # Knotenkoordinatenmatrix
        self.node2cell      = self.gui.mesh_cells                                                                         # Knotenkoinzidenzmatrix 
        self.node2dof       = np.arange(len(self.nodecoordinate[:,0])*3).reshape((len(self.nodecoordinate[:,0]),3))       # ElementDOFmatrix



        # HEX8 meshio -> FEM Referenzordnung

        for e in range(len(self.node2cell)):

             nodes = self.node2cell[e]

             coords = self.nodecoordinate[nodes]

             xmin = np.min(coords[:,0])
             xmax = np.max(coords[:,0])

             ymin = np.min(coords[:,1])
             ymax = np.max(coords[:,1])

             zmin = np.min(coords[:,2])
             zmax = np.max(coords[:,2])


             new_order = [
                 np.argmin(np.sum((coords-[xmin,ymin,zmin])**2,axis=1)),
                 np.argmin(np.sum((coords-[xmax,ymin,zmin])**2,axis=1)),
                 np.argmin(np.sum((coords-[xmax,ymax,zmin])**2,axis=1)),
                 np.argmin(np.sum((coords-[xmin,ymax,zmin])**2,axis=1)),

                 np.argmin(np.sum((coords-[xmin,ymin,zmax])**2,axis=1)),
                 np.argmin(np.sum((coords-[xmax,ymin,zmax])**2,axis=1)),
                 np.argmin(np.sum((coords-[xmax,ymax,zmax])**2,axis=1)),
                 np.argmin(np.sum((coords-[xmin,ymax,zmax])**2,axis=1)),
             ]


             self.node2cell[e] = nodes[new_order]





        if self.nodecoordinate is None or self.node2cell is None:
            print("Kein Mesh geladen!")                         # wie kann man ab hier den ablauf aus lineare FEM abbrechen??
            return None
        nodedata = (self.nodecoordinate, self.node2cell,self.node2dof)






        return nodedata
    



    def mat_model(self):

        self.mat_tensor = np.zeros((6,6))
        if self.gui.mat_combo.currentText()=="isotrop" and self.gui.mat_combo2.currentText()=="linear elastisch":
            self.mat_values = self.gui.mat_values 
            E       =   float(self.mat_values["E"][0].text())   # E-Modul
            v       =   float(self.mat_values["v"][0].text())   # Querkontraktionszahl
            lambdaa =   E*v/((1+v)*(1-2*v))
            mu      =   E/(2+2*v)
            self.mat_tensor[0][0] = self.mat_tensor[1][1] = self.mat_tensor[2][2] = lambdaa + 2*mu
            self.mat_tensor[3][3] = self.mat_tensor[4][4] = self.mat_tensor[5][5] = mu
            self.mat_tensor[0][1] = self.mat_tensor[0][2] = self.mat_tensor[1][0] = self.mat_tensor[2][0]  = lambdaa
            self.mat_tensor[1][2] = self.mat_tensor[2][1] = lambdaa
        self.mat_tensor = np.asarray(self.mat_tensor)


        


        



    def jacobian(self):

        #nachher aus eingabven automatisch entnehmen
        #self.n_dim  = 3
        #self.n_GP   = 8
        #self.n_SF   = 8
        #self.n_elements = 40
        #self.n_el_node = 8 # wie viele Knoten ein bestimmtes element hat zb. hex 8 enthält 8 Knoten 
        k=0

        

        #diese drei attribute werden durch calc jacobian geändert
        self.Jacobi     = np.zeros((3,3,self.n_elements,self.n_GP[0]))   #korrigieren hier muss noch bessere lösung gefunden werden das self.n_GP[0] nicht richtig ist
        self.det_Jacobi = np.zeros((self.n_elements,self.n_GP[0]))
        self.inv_Jacobi = np.zeros((3,3,self.n_elements,self.n_GP[0]))

        vector = np.zeros((8))
        d_shape_fun_eval = np.zeros((192))  # auswertung der 8 formfunktionen für die 3 ableitungsrichtungen an den 8 Gauspunkten (später flexible für verschiedene elemente und anzahl an Gauspunkten)

        
        for i in range(0,len(self.n_SF)):
            for n_shape_fun in range(self.n_SF[i]):   
                for d_xyz in range(0,self.n_dim):
                    d_shape_fun = self.d_shape_hex_8[:,d_xyz]
                    for GP in range(0,self.n_GP[i]):
                        d_shape_fun_eval[GP+d_xyz*8+n_shape_fun*24] = d_shape_fun[n_shape_fun](self.gaus_pos_hex8_2[GP][0],self.gaus_pos_hex8_2[GP][1],self.gaus_pos_hex8_2[GP][2])

        for element in range(0,self.n_elements):                                    #   abhängig von der anzahl an Elementen (die hexader sind, am ende schleife um ganze def_el_stiffness für die verschiedenen elementtypen)
            for GP in range(0,self.n_GP[0]):                                        #   muss noch angepasst werden schleife darüber darf nicht direkt über alle element gehen sondern muss blockweise elements abarbeiten um auch unterschiedliche elemente innerhalb eines meshes zu unterstützen
                for x_y_z in range(0,self.n_dim):
                    for d_xyz in range(0,self.n_dim):
                        for el_node in range(0,self.n_el_node[0]):
                            gl_node = self.node2cell[element,el_node]
                            coordinate = self.nodecoordinate[gl_node,x_y_z]         # el_node: Elementknoten(jedes hex element hat die Knoten 0-8) gl_node (globale elementnnummer)
                            vector[el_node] = coordinate*d_shape_fun_eval[GP+d_xyz*8+el_node*24]
                        self.Jacobi[x_y_z,d_xyz, element, GP] = sum(vector)

            for GP in range(0,self.n_GP[0]):
                Jacobi_2d             = self.Jacobi[:,:,element,GP]
                self.det_Jacobi[element,GP]   = np.linalg.det(Jacobi_2d)
                self.inv_Jacobi[:,:,element,GP] = np.linalg.inv(Jacobi_2d)
                k+=1

            



    def deformation_tensor(self,model):

        #wie erstelle ich die matrix des defomrationstensors am anfang mit nur 0en? zuvor noch schauen welcher elementtyp  meisten zeilen braucht und das als struktur für deformationstensor wählen
        self.deform_tensor = np.zeros((2*self.n_dim,self.n_SF[0]*self.n_dim,self.n_GP[0],self.n_elements))

        # noch überlegen wie die logik flexibel werden kann so dass sie auch für meshes mit unterschiedlichen Elementtypen funktioniert --> unterschiedliche Anzahl an einträgen in Verzerrungstensor
        # bis jetzt nur für hexaeder --> wie soll formdes deform tensors aussehen wenn dieser felxibel für alle elemente anzuwenden ist??


        if model=="linear":   # lineares Verzerrungsmaß später auch erweiterbar sodass nichtlineare Verzerrungsmase eingefügt werden können         
            for j in range(0,len(self.el_type)):
                for element in range(0,self.n_elements):
                    if self.el_type[j]=="hex8":
                        for GP in range(0,self.n_GP[0]):
                            for shape_function in range (0,self.n_SF[0]):
                                #zahl 0 gibt an nach welcher variable die shapefunktion abgeleitet wurde                                                 #zahl 0 gibt an um welche ableitung es sich handelt
                                dN_dxi    = self.d_shape_hex_8[shape_function][0](self.gaus_pos_hex8_2[GP][0],self.gaus_pos_hex8_2[GP][1],self.gaus_pos_hex8_2[GP][2])
                                dN_deta   = self.d_shape_hex_8[shape_function][1](self.gaus_pos_hex8_2[GP][0],self.gaus_pos_hex8_2[GP][1],self.gaus_pos_hex8_2[GP][2])
                                dN_dzeta  = self.d_shape_hex_8[shape_function][2](self.gaus_pos_hex8_2[GP][0],self.gaus_pos_hex8_2[GP][1],self.gaus_pos_hex8_2[GP][2])
                                
                                inv_Jacobi = self.inv_Jacobi[:,:,element, GP]

                                grad_nat = np.array([dN_dxi, dN_deta, dN_dzeta])
                                grad_global = inv_Jacobi @ grad_nat

                                dN_dx = grad_global[0] 
                                dN_dy = grad_global[1]
                                dN_dz = grad_global[2]

                                col = shape_function * self.n_dim

                                # ε_xx
                                self.deform_tensor[0, col + 0, GP, element] = dN_dx

                                # ε_yy
                                self.deform_tensor[1, col + 1, GP, element] = dN_dy

                                # ε_zz
                                self.deform_tensor[2, col + 2, GP, element] = dN_dz

                                # γ_xy
                                self.deform_tensor[3, col + 0, GP, element] = dN_dy
                                self.deform_tensor[3, col + 1, GP, element] = dN_dx

                                # γ_yz
                                self.deform_tensor[4, col + 1, GP, element] = dN_dz
                                self.deform_tensor[4, col + 2, GP, element] = dN_dy

                                # γ_xz
                                self.deform_tensor[5, col + 0, GP, element] = dN_dz
                                self.deform_tensor[5, col + 2, GP, element] = dN_dx


                    else: # hier dass die gleichen schleifen nur formfunktionen werden entsprechen von anderem element herangezogen und größe des verzerrungstensors ist anders
                        print("fail")
        else:
            print("FEHLER leider hat die if bedingung in calculate deformation tensor nicht funktioniert")



        # nicht irritieren lassen durch viele null einträge das liegt an der struktur des defomrationstensors
        # print("Deformationstensor für alle element wurde erstellt")
        # print(self.deform_tensor[:,:,0,0])


 
    def def_stiffness(self):

        # schauen wieso self.n_SF oder self.n_elements hier nicht vorhanden ist obwohl es in Jacobian definiert wurde und deformation tensor die variable auch kennt
        el_stiffness_GP      = np.zeros((24,24))
        self.el_stiffness     = np.zeros((24,24))                   #   jetzt für lineare Hexaeder (8 Knoten) und 3 dimensionen (3d FEM)

        self.jacobian()
        self.deformation_tensor("linear")                           #   später wenn mehrere verzerrungstensoren vorhanden sind als variable aus der GUI ziehen
        
        for element in range(0,self.n_elements):
            self.el_stiffness     = np.zeros((24,24))
            for GP in range(0,self.n_GP[0]):
                deform                    = self.deform_tensor[:,:,GP,element]
                transpose_deform          = deform.T
                el_stiffness_GP[:,:]      = transpose_deform@self.mat_tensor@deform*self.det_Jacobi[element,GP]*self.gaus_weight_hex8_2[GP]        
                self.el_stiffness[:,:]   += el_stiffness_GP


            self.COO_stiffness_matrix = self.matrices_to_COO(element,self.COO_stiffness_matrix, self.el_stiffness)    #übergibt die matrix damit bereits befüllte element schon vorhanden sind und neue elemente nur hinzugefügt werden         

        self.row, self.column, self.value  = self.COO_to_CSR(self.COO_stiffness_matrix)     # für direkte Solver (CSR besser für iterative Solver)




    def matrices_to_COO(self,element,COO_matrix, element_matrix):
        

        # sollte jetzt flexibel auch für massen und dämpfungsmatrix einsetzbar sein
        #24x24 matrix (len(self.el_stiffness)=24!)  8 Knoten mit je 3 Freiheitsgraden
        #1-3 1. Knoten 4.-6 2. Knoten ....
        n_nodes=int((len(element_matrix)/self.n_dim))

        for i in range(0,n_nodes):
            for j in range(0,self.n_dim):
                for k in range(0,n_nodes):                        # jetzt für lineare hexaeder (woher kann man für jedes element die anzahl an Knoten herbekommen??)
                    for l in range(0,self.n_dim):
                        if abs(element_matrix[i*self.n_dim+j, k*self.n_dim+l])>0:        #  1e-9:  # numerische artefakte die nahe null liegen werden nicht berücksichtigt
                            if self.node2dof[self.node2cell[element,i],j]!=self.node2dof[self.node2cell[element,k],l]:
                                COO_row    = self.node2dof[self.node2cell[element,i],j]          # elementnode i -> globalnode          # Reihe 
                                COO_column = self.node2dof[self.node2cell[element,k],l]          # elementnode i -> globalnode          # Spalte
                                value      = element_matrix[i*self.n_dim+j, k*self.n_dim+l]                                       # Wert                       
                                COO_matrix.append([COO_row,COO_column,value])

                            else:
                                COO_row    = self.node2dof[self.node2cell[element,k],l]          # elementnode i -> globalnode          # Reihe 
                                COO_column = self.node2dof[self.node2cell[element,i],j]
                                value      = element_matrix[i*self.n_dim+j, k*self.n_dim+l]                                       # Wert                       
                                COO_matrix.append([COO_row,COO_column,value])
                        
        return COO_matrix  
    

    
    def COO_to_CSR(self, COO_matrix):

        zwischenspeicher = []   #[[entity_of_COO_matrix],[value]]
        row    = []
        column = []
        value  = []
        COO_row= [] 
        j=0
        
        
        n_dof = self.node2dof.max() + 1

        rows, cols, data = zip(*COO_matrix)
        K_csr = coo_matrix((data, (rows, cols)), shape=(n_dof, n_dof)).tocsr()

        row = K_csr.indptr.tolist()    # Zeiger auf Zeilenanfang
        column = K_csr.indices.tolist()  # Spaltenindizes der Werte
        value = K_csr.data.tolist()      # Werte


        # soll nur COO zu CSR Format umschreiben (egal ob matrix oder vektor damit für Kraftvektor auch verwendet werden kann)
        #while j < len(COO_matrix):
            #k=0
            #print("j:",j)
        #for i in range(0,len(COO_matrix)): #alle einträge der COO_matrix werden durchsucht
        #    m=0
        #    zwischenspeicher.append([COO_matrix[i][0],COO_matrix[i][1],COO_matrix[i][2]])   # Spalte und Wert werden angehängt
        #    m+=1
                      
        #zwischenspeicher=sorted(zwischenspeicher, key=lambda row: row[0])  
      
        #idx = np.lexsort((np.asarray(zwischenspeicher)[:,1], np.asarray(zwischenspeicher)[:,0]))
        #filtered_zwischenspeicher = np.array(zwischenspeicher)[idx]

        #ALTE VERSION DIE NUR NACH ROW SORTIERT HAT
        #filtered_zwischenspeicher = list([row for row in zwischenspeicher if row[2] != 0]) # entfernt alle einträge bei denen value (also zwischenspeicher[][]==0)


        #COO_row.append(filtered_zwischenspeicher[0][0])
        #column.append(filtered_zwischenspeicher[0][1])
        #value.append(filtered_zwischenspeicher[0][2])
      
        
        #for i in range(1,len(filtered_zwischenspeicher)):
            #if filtered_zwischenspeicher[i][0]==row[len(row)-1] and filtered_zwischenspeicher[i][1]==column[len(column)-1]:       #prüfen ob Betrag in existierenden Listeneintrag summiert werden muss oder neuer erstellt werden soll
        #    if filtered_zwischenspeicher[i][0]==COO_row[len(COO_row)-1] and filtered_zwischenspeicher[i][1]==column[len(column)-1]:
        #        value[len(value)-1] += filtered_zwischenspeicher[i][2]
                #print("Betrag muss aufsummiert werden")
                #print("Eintrag COO row",len(COO_row))
        #    else:
        #        COO_row.append(filtered_zwischenspeicher[i][0])
        #        column.append(filtered_zwischenspeicher[i][1])
        #        value.append(filtered_zwischenspeicher[i][2])
                 
        #row.append(0)

        #in CSR Format bringen
        #for i in range(0,len(list(set(np.array(filtered_zwischenspeicher)[:,0])))):    #anzahl an reihen  (testbeispiel = 297) 
        #    plus_row = list(np.array(column)).count(i)
        #    row.append(row[len(row)-1]+plus_row)



        #n_rows = int(max(COO_row)) + 1
        #for i in range(0, n_rows):
        #    count = COO_row.count(i)
        #    row.append(row[-1] + count)

        #COO_row_np = np.array(COO_row)
        #row = np.zeros(int(max(COO_row_np))+2, dtype=int)
        #for r in range(int(max(COO_row_np))+1):
        #    row[r+1] = row[r] + np.sum(COO_row_np == r)

        #n_dof = self.node2dof.size
        #row = np.zeros(n_dof+1, dtype=int)
        #for i in range(n_dof):
        #    row[i+1] = row[i] + np.sum(COO_row_np == i)

        return row, column, value 










    def apply_bc(self,method):

        #noch if bedinung einfüguen (falls druck vorgegeben wird kann der erste teil entfallen)
        face_Jacobian = []
        globalForceVector_COO = []
        globalForceVector = np.zeros(np.shape(self.node2dof)[0]*np.shape(self.node2dof)[1])
   
        
        #shapefunction = self.shape_hex8
        #gaus_point_coordinate = self.gaus_pos_hex8_2
        #gaus_weight = self.gaus_weight_hex8_2


        # if bedigung um anhand der Knoten zu schauen um welche Fläche des Referenzelement es sich handelt und welche Integrationsvariablen für die Flächenberechnung verwendet werden müssen
        
        # Schleife über alle Randbedingungen 
        for i in range(0,len(self.gui.boundary_conditions)):
            item = self.gui.boundary_conditions[i]["type"]

            

            
            #############################################################################################################################################
            # hier noch Fallunterschied dass dieser Teil nur ausgeführt wird wenn es sich um eine Kraftrandbedingung handelt (bei aufgebrachtem Druck ist keine Flächenberechnung nötig/ bei verschiebungs RB ändert das nicht den Kraftvektor  (Eigengewicht muss nur berücksichtigt werden können))
            #############################################################################################################################################
            if item == "Kraftrandbedingung":
                const_pressure_bc, integration_variables, free_coords, face_Jacobian, faces = self.get_constant_pressure(i)
                globalForceVector_COO = self.def_force_vector(faces, integration_variables,free_coords,const_pressure_bc,face_Jacobian,i)

                print("faces an denen die Kraftrandbedingung angewandt wird:", faces)

                for k in range(0,len(globalForceVector_COO)):
                    globalForceVector[globalForceVector_COO[k][1]]+=globalForceVector_COO[k][0]
                    #print("freiheitsgrad des globale Force Vektor eintrags",globalForceVector_COO[k][1])

                print("globaler Force Vector",globalForceVector)


               

                

           # elif item =="Druckrandbedingung":    
           #     integration_variables =
           #     free_coords =
           #    globalForceVector_COO = self.def_force_vector(self, faces, integration_variables,globalForceVector_COO)

           #    for k in range(0,len(globalForceVector_COO)):
           #        globalForceVector[globalForceVector_COO[k][1]]=globalForceVector_COO[k][0]

            elif item =="Lager" or item=="Verschiebungsrandbedingung":
                self.apply_dirichlet_bc(i,item,"penalty")

            self.globalForceVector = globalForceVector

        FEMdata = (self.row, self.column, self.value, np.asarray(globalForceVector))


        print("gesamtkraft: ",sum(globalForceVector))

        return FEMdata 





    def apply_dirichlet_bc(self,i,type,method):

        faces = self.gui.boundary_conditions[i]["faces"]

        if method=="penalty":
            if type=="Lager":
                blocked_DOF = []
                for j in range(0,len(faces)):
                    bearing_type = self.gui.boundary_conditions[i]["tree_items"]["bearing_type"]
                    facenode = list(self.gui.boundary_conditions[i]["facenodes"])[j]        
                    facenodes = [int(x) for x in str(facenode).replace('[','').replace(']','').split()]


                    if bearing_type =="Feste Einspannung":
                        # definieren der zu erzielenden RB  also hier alle FHG von allen facenodes =0 setzen (in anderen Fällen bleiben manche freiheitsgrade von machen Knonten flexibel)
                        for k in range(0,len(facenodes)):
                            blocked_DOF_face =[int(x) for x in str(self.node2dof[facenodes[k]]).replace('[','').replace(']','').split()]
                            blocked_DOF.extend(blocked_DOF_face)
                    #elif bearing_type =="Gelenklager": # hier werden zb direkt die Knoten vorgeben und müssen nicht erst aus den Flächen ermittelt werden
                
                blocked_DOF = list(set(blocked_DOF))  #hier stehen alle gesperrten DOF's drinnen aus allen FHG (achtung nur bei lager sinnvoll bei verschiebung muss wert vorgegeben werden(seperat))


                _lambda = float(1e12)  # oder parametrisierbar

                for dof in blocked_DOF:
                    row_start = self.row[dof]
                    row_end   = self.row[dof + 1]
                    for idx in range(row_start, row_end):
                        if self.column[idx] == dof:
                            self.value[idx] = self.value[idx]+_lambda
                            break


            #elif bearing_type =="Verschiebungsrandbedingung":
                #heraussuchen welche Knoten verschoben werden 
                #aufbringen der Sperrung der DOF Randbedingungen  
 

    def def_force_vector(self,faces,integration_variables,free_coords,const_pressure_bc,face_Jacobian,i):

        shapefunction = self.shape_hex8
        gaus_point_coordinate   = self.gaus_pos_hex8_2
        gaus_weight = self.gaus_weight_hex8_2
        localForceVector = []

        for j in range(0,len(faces)):

            shapefun_use   = []
            forceVectorDOF = []
            
            # prüfen welche Formfunktionen null werden 
            for k in range(0,self.n_SF[0]):                     # nur im Basisbeispiel gültig weil alle elemente lineare hexaeder sind 
                value = shapefunction[k](integration_variables[0],integration_variables[1],integration_variables[2])
                if value !=0:
                    shapefun_use.append(k)

            elementForceVector_COO_face = np.zeros((len(shapefun_use)*3))
            N_T = np.zeros((len(shapefun_use)*3,3)) 
            
            # noch schauen wie flexibilisieren (unterschiedliche anzahl an shapefunktions nötig) um auch flexibel auf andere integrationsschema anwendbar zu sein
            
             
            k=0
            for GP1 in range(0,2):
                for GP2 in range(0,2):
                    integration_variables[free_coords[0]] = gaus_point_coordinate[k][free_coords[0]]
                    integration_variables[free_coords[1]] = gaus_point_coordinate[k][free_coords[1]]
                    k+=1

                    #print(integration_variables)

                    

                    N_T[0][0] = N_T[1][1]  = N_T[2][2]  = shapefunction[shapefun_use[0]](integration_variables[0],integration_variables[1],integration_variables[2])
                    N_T[3][0] = N_T[4][1]  = N_T[5][2]  = shapefunction[shapefun_use[1]](integration_variables[0],integration_variables[1],integration_variables[2])
                    N_T[6][0] = N_T[7][1]  = N_T[8][2]  = shapefunction[shapefun_use[2]](integration_variables[0],integration_variables[1],integration_variables[2])
                    N_T[9][0] = N_T[10][1] = N_T[11][2] = shapefunction[shapefun_use[3]](integration_variables[0],integration_variables[1],integration_variables[2])

                    elementForceVector_COO_face += (N_T@const_pressure_bc*face_Jacobian[j][GP1+2*GP2]*gaus_weight[GP1]*gaus_weight[GP2]).squeeze()

            
            global_facenodes = list(self.gui.boundary_conditions[i]["facenodes"])[j]
            global_facenodes = global_facenodes.strip("[]").split()
            global_facenodes = [int(n) for n in global_facenodes]


            

            forceVectorDOF.extend(self.node2dof[global_facenodes].flatten().tolist())
            localForceVector.extend(zip(elementForceVector_COO_face.tolist(),forceVectorDOF))

            #print("forceVectorDOF",forceVectorDOF)
            #print("localFroceVector",localForceVector)
        

        return localForceVector


    
    




    def get_face_area(self, integration, elementtype,local_element, element_nodes,local_face_id):

        # elementnodes: alle Knotennummer des Elements auf dem die Fläche liegt deren Flächeninhalt berechnet werden soll

        x1 = []
        x2 = []
        x3 = []
        x = []

        dN_dxi   = []
        dN_deta  = []
        dN_dzeta = []
 
        x_d = np.empty((2,3))

        # Knoten der 8 punkte die das element definieren auf dem die fläche liegt 
        for i in range(0,len(element_nodes)):
            x1.append(list(self.nodecoordinate[element_nodes[i]])[0])
            x2.append(list(self.nodecoordinate[element_nodes[i]])[1])
            x3.append(list(self.nodecoordinate[element_nodes[i]])[2])
            x.append((list(float(i) for i in self.nodecoordinate[element_nodes[i]])))


        if elementtype == 'hex8':
            d_shapefunction          = self.d_shape_hex_8
            gaus_point_coordinate    = self.gaus_pos_hex8_2
            gaus_weight              = self.gaus_weight_hex8_2

            if local_face_id in [0,2,5] :
                fx_coordinate = -1
            elif local_face_id in [1,3,4] :
                fx_coordinate = 1
           
            if local_face_id in [0, 1]:
                fixed_coord = 2
                free_coords = (0, 1)

            elif local_face_id in [2, 4]:
                fixed_coord = 1
                free_coords = (0, 2)

            elif local_face_id in [3, 5]:
                fixed_coord = 0
                free_coords = (1, 2)

            print("local face id",local_face_id)
            print("free coord",free_coords)
            print("fixed coords",fixed_coord)
            print("fx coordinate",fx_coordinate)

        #in den elif kommen dann die anderen elementtypen rein und werden d_shapefunction und shapefunction zugeordnet

        # dieser teil nur gültig für lineare hexaeder <a     / evtl mit if bedingung variable blegen --> d_shapefunction = self.d_shape_hex_8 --> in funktion d_shapefunction[shapefunction[0](.....)] 
        # Auswertung aller Formfunktionen an den Gauspunkten abhängig davon welche formfunktionen (immer nur 2 von 3 werden gebraucht) gebraucht werden
        face_area = 0
        face_Jacobian = []

        k=0
        

        for GP_1 in range(0,2):
            for GP_2 in range(0,2):       #2x2=4gauspunkte (wegen 2d)   self.n_GP[0]=8 wegen 3d -->evtl self.n_GP_3d self.n_GP_2d=0 einführen 
                integration_variable = [0.0, 0.0, 0.0]
                integration_variable[fixed_coord]    = fx_coordinate
                integration_variable[free_coords[0]] = gaus_point_coordinate[k][free_coords[0]]
                integration_variable[free_coords[1]] = gaus_point_coordinate[k][free_coords[1]]
                
                k+=1

                integration_variables = np.asarray(integration_variable)
                
                for i in range(0,2):
                    if integration[i] == 1:
                        for shape_function in range(0,self.n_SF[0]):                           
                            dN_dxi.append(d_shapefunction[shape_function][0](integration_variables[0],integration_variables[1],integration_variables[2]))      
                        #print("dN_dxi",dN_dxi)
                        x_d[i][0] = sum(np.asarray(dN_dxi)*np.asarray(x1))       # die unterschiedlichen GP werden 
                        x_d[i][1] = sum(np.asarray(dN_dxi)*np.asarray(x2))
                        x_d[i][2] = sum(np.asarray(dN_dxi)*np.asarray(x3))                           
                        dN_dxi.clear()

                    elif integration[i] == 2:
                        for shape_function in range(0,self.n_SF[0]):
                            dN_deta.append(d_shapefunction[shape_function][1](integration_variables[0],integration_variables[1],integration_variables[2]))
                        #print("dN_deta",dN_deta)
                        x_d[i][0] = sum(np.asarray(dN_deta)*np.asarray(x1))
                        x_d[i][1] = sum(np.asarray(dN_deta)*np.asarray(x2))
                        x_d[i][2] = sum(np.asarray(dN_deta)*np.asarray(x3))
                        dN_deta.clear()

                    elif integration[i] == 3:
                        for shape_function in range(0,self.n_SF[0]):
                            dN_dzeta.append(d_shapefunction[shape_function][2](integration_variables[0],integration_variables[1],integration_variables[2]))
                        #print("dN_dzeta",dN_dzeta)
                        x_d[i][0] = sum(np.asarray(dN_dzeta)*np.asarray(x1))
                        x_d[i][1] = sum(np.asarray(dN_dzeta)*np.asarray(x2))
                        x_d[i][2] = sum(np.asarray(dN_dzeta)*np.asarray(x3))
                        dN_dzeta.clear()


                face_Jacobian.append(np.linalg.norm(np.cross(x_d[0][:],x_d[1][:])))
                face_area += face_Jacobian[-1]*gaus_weight[GP_1]*gaus_weight[GP_2]

                         
        
        return face_area, integration_variables ,free_coords, face_Jacobian
    





    def get_constant_pressure(self,i):

        element_type  = []
        face_area     = []
        face_Jacobian = []


        faces = next(d.get('faces') for d in self.gui.boundary_conditions if 'faces' in d)

        print("faces",faces)

        
        #print("ELEMENTNODES", list(list(self.gui.boundary_conditions[i]["elementnodes"])))

        


        for j in range(0,len(faces)):
            element_type.append(next(d.get('elementtype') for d in self.gui.boundary_conditions if 'elementtype' in d))
            local_element, local_face_id = self.gui.surface_face_mapping[list(faces)[j]]

            print("local element",local_element)
            print("local_face_id",local_face_id)

             
            # Knoten der faces bekommen --> vergleich mit Knoten aus element --> welche freiheitsgrade wurden verwendet un müssen für die Berechnung der Fläceh herangezogen werden 
            # if bedingung später um die anderen Elementtypen erweitern

            
            if str(list(element_type)[j]) == '{<CellType.HEXAHEDRON: 12>}':
                if local_face_id == 0 or  local_face_id == 1:                               # xi und eta wird integriert
                    integration = [1,2]
                elif local_face_id == 2 or local_face_id == 4:                              # xi und zeta wird integriert
                    integration = [1,3]
                elif local_face_id == 3 or local_face_id == 5:                              # eta und zeta wird integriert
                    integration = [2,3]
                elementtype = "hex8"
            else: 
                print("elementtyp nicht gefunden")

            element_nodes = list(list(self.gui.boundary_conditions[i]["elementnodes"])[j][:])
            area, integration_variables, free_coords, face_Jacob  = self.get_face_area(integration,elementtype,local_element, element_nodes, local_face_id)
            face_area.append(area)
            face_Jacobian.append(face_Jacob)

        Full_face_area = sum(np.asarray(face_area))
        const_pressure_bc = (np.asarray(self.gui.boundary_conditions[i]["values"])/Full_face_area).reshape(3,1)

        print("full face area",Full_face_area)
        print("face area",face_area)
        print("constant pressure bc",const_pressure_bc)

        return const_pressure_bc, integration_variables ,free_coords, face_Jacobian, faces
    


