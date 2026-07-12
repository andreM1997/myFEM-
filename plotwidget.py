import pyvista as pv
from PyQt6.QtWidgets import (
    QApplication, QMainWindow, QWidget, QVBoxLayout, QHBoxLayout,
    QPushButton, QFileDialog, QLabel, QComboBox, QListWidget, QMessageBox,
    QTreeWidget, QTreeWidgetItem, QLineEdit, QInputDialog, QCheckBox, QSizePolicy,QSlider
)
from PyQt6.QtGui import QFont, QColor, QBrush
from PyQt6.QtCore import Qt
from pyvistaqt import QtInteractor
import numpy as np




class PlotWidget(QWidget):

    def __init__(self,parent):

        super().__init__(parent)
        layout = QVBoxLayout(self)
        self.plotter = QtInteractor(self)
        layout.addWidget(self.plotter.interactor)

        self.plotter.set_background("white")
        self.plotter.add_axes()

        # Slider für Verschiebungsskalierung
        self.disp_slider = QSlider(Qt.Orientation.Horizontal)
        layout.addWidget(self.disp_slider)
        self.disp_slider.setRange(1, 100)  # Faktor 0.01 bis 1.0 z.B.
        self.disp_slider.valueChanged.connect(self.update_deformed_mesh)
        self.disp_slider.hide()  # nur zeigen, wenn der Deformationsplot aktiv ist

        self.current_mesh = None
        self.current_disp = None

    def show_mesh(self, mesh): 
        self.plotter.clear()
        self.plotter.add_mesh(mesh)
        self.current_mesh = mesh
        self.current_disp = None
        self.disp_slider.hide()  # Slider nur für Verschiebung relevant

    def show_deformed_mesh(self, mesh, displacements,scalar_field,item):
        self.plotter.clear()
        self.current_mesh = mesh
        self.current_item = item
        self.current_disp = displacements.reshape(int(len(displacements)/3),3)
        self.scalar_field = scalar_field.reshape(int(len(displacements)/3),3)   #Spannungen, x,y,z Verschiebungen, hauptspannung..... 
        self.update_deformed_mesh()
        self.disp_slider.show()  # Slider einblenden


    def update_deformed_mesh(self):

        if self.current_mesh is None or self.current_disp is None:
            return
        factor = 10*self.disp_slider.value()**2
        deformed_mesh = self.current_mesh.copy()
        deformed_mesh.points += self.current_disp * factor

        print("Faktor des disk sliders",factor)


        
        # hier muss dann irgendwie noch verlinkung zur berechnung von Spannungen und ähnlichem Erfolgen 
        if self.current_item =="x Verschiebung":
            coord = 0
        elif self.current_item =="y Verschiebung":
            coord = 1
        elif self.current_item =="z Verschiebung":
            coord = 2
        else:
            print("keine der Verfügbaren Verschiebungsoptionen wurde ausgewählt hier sollte dann weiter zu funktionen geleitete werden die Spannungen oder ähnliches zurück geben")

      
        deformed_mesh['Data'] = self.scalar_field[:,coord]


        print("min ux:", np.min(self.current_disp[:,0]))
        print("max ux:", np.max(self.current_disp[:,0]))
        print("min uy:", np.min(self.current_disp[:,1]))
        print("max uy:", np.max(self.current_disp[:,1]))
        print("min uz:", np.min(self.current_disp[:,2]))
        print("max uz:", np.max(self.current_disp[:,2]))


        
        self.plotter.clear()
        self.plotter.add_mesh(deformed_mesh,scalars='Data',cmap="jet",show_edges=True,smooth_shading=True,scalar_bar_args={'title': self.current_item})