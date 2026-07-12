import sys
import os
import numpy as np
import meshio
import pyvista as pv
from PyQt6.QtWidgets import (
    QApplication, QMainWindow, QWidget, QVBoxLayout, QHBoxLayout,
    QPushButton, QFileDialog, QLabel, QComboBox, QListWidget, QMessageBox,
    QTreeWidget, QTreeWidgetItem, QLineEdit, QInputDialog, QCheckBox, QSizePolicy
)
from PyQt6.QtGui import QFont, QColor, QBrush
from PyQt6.QtCore import Qt
from pyvistaqt import QtInteractor
#from GUI import FEMGUI
import datapreparation
import Solver




class Controller:
    def __init__(self,gui):
        self.gui            = gui
        self.mesh           = None
        self.mesh_window    = None
        self.displacement   = None
        #self.plot_widget = PlotWidget()

        #self.plot_widget.show()


        #self.nodedata   = nodedata      # node2dof, nodecoordinate, node2cell
        #self.FEMdata    = FEMdata       # steifigkeitmatrix, massenmatrix, dämpfungsmatrix, Kraftvektor
        #self.solverdata = solver_data   # welcher algortihmus zb LU Zerlegung oder ähnliches soll genutzt werden?? --> aus GUI entnehmen 




        ##############################################################################################################
        #######  STEIFIGKEITSMATRIX IST VOR UMWANDLUNG IN COO FORMAT UND DANACH CSR FORMAT FÜR JEDES ELEMENT SYYMETRISCH
        #######  PRÜFEN OB COO FALSCH BEFÜLLT WIRD; FALSCH IN CSR UMGEWANDELT WIRD UND OB DIE PENALTY METHODE RICHTIG AUF DIE STIEIFGKEITSMATRIX ANGEWANDT WIRD
        ################################################################################################################
        
    def static_structural(self):

        # später mehrere verschiedene Objekte für anfangsverwaltung der daten, aufstellen des Gleichungssystems, Lösen des Gleichungssystems
        lineare_FEM = datapreparation.def_data(self.gui)   #definiert das objekt der klasse def_data
        nodedata = lineare_FEM.nodematrices() 
        lineare_FEM.mat_model()
        lineare_FEM.def_stiffness()
        FEMdata = lineare_FEM.apply_bc("penalty")   # oder None wenn keine dirichlet RB/ eigentlich aus GUI herausziehen
        
        #solverdata = ("static","linear","Cholesky")  # nur vorübergehen später aus GUI herausziehen
        #solver = Solver.solveFEM(nodedata,FEMdata,solverdata)
        solver = Solver.solveFEM(lineare_FEM)
        self.displacements = solver.solve_static_structural()
        print("die displacements wurden berechnet")

        
        
        #solveFEM.static_structural()



    # def modal_analysis(self):
    # def harmonic analysis(self):