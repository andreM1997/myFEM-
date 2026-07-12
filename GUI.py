import sys
import numpy as np
import meshio
import pyvista as pv
from PyQt6.QtWidgets import (
    QApplication, QMainWindow, QWidget, QVBoxLayout, QHBoxLayout,
    QPushButton, QFileDialog, QLabel, QComboBox, QListWidget, QMessageBox,
    QTreeWidget, QTreeWidgetItem, QLineEdit, QInputDialog, QCheckBox, QSizePolicy,QMenu
)
from PyQt6.QtGui import QFont, QColor, QBrush
from PyQt6.QtCore import Qt
from pyvistaqt import QtInteractor
from plotwidget import PlotWidget




VTK_CELL_TYPES = {
    "tetra": pv.CellType.TETRA,
    "hexahedron": pv.CellType.HEXAHEDRON,
    "wedge": pv.CellType.WEDGE,
    "pyramid": pv.CellType.PYRAMID,
}

ELEMENT_FACES = {
    pv.CellType.TETRA: [[0,1,2],[0,1,3],[0,2,3],[1,2,3]],
    pv.CellType.HEXAHEDRON: [[0,1,2,3],[4,5,6,7],[0,1,5,4],[1,2,6,5],[2,3,7,6],[3,0,4,7]],
    pv.CellType.WEDGE: [[0,1,2],[3,4,5],[0,1,4,3],[1,2,5,4],[2,0,3,5]],
    pv.CellType.PYRAMID: [[0,1,2,3],[0,1,4],[1,2,4],[2,3,4],[3,0,4]],
}



class FEMGUI(QMainWindow):
    def __init__(self,controller):
        super().__init__()
        self.controller=controller
        self.setWindowTitle("myFEM")
        self.setGeometry(60, 60, 1800, 900)

        # Meshdaten
        self.mesh           = None
        self.mesh_points    = None
        self.mesh_cells     = None
        self.nodes_polydata = None
        self.surface_faces  = None
        self.meshio_mesh    = None
        self.surface_face_mapping = {}
        self.elements       = []
        self.el_type        = []
        self.el_type_num    = []
        self.n_SF           = []
        self.n_elements     = []
        self.n_el_node      = []

        # Selektion
        self.selected_nodes = set()
        self.selected_faces = set()
        self.node_selection_actors = {}
        self.face_selection_actors = {}

        # Randbedingungen
        self.boundary_conditions = []
        self.current_bc_index = None
        self.selection_mode = "Knoten"
        self.arrow_actor         = []
        self.arrow_nodes         = []

        # Drag detection
        self.mouse_pressed_pos = None
        self.drag_threshold = 5  # Pixel

        # Materialwerte
        self.mat_values = {"E": [], "v": []}     #anpassen um auch andere materialmodelle abbilden zu können
        self.mat_combo = None
        self.mat_combo2 = None

        #Ergebnisse
        self.evaluated_result = []

        self.init_ui()

    def init_ui(self):
        central = QWidget()
        self.setCentralWidget(central)
        main_layout = QHBoxLayout(central)

        # Steuerung links
        ctrl_layout = QVBoxLayout()
        main_layout.addLayout(ctrl_layout, 0)

        self.load_btn = QPushButton("MSH Datei laden")
        self.load_btn.clicked.connect(self.on_load_mesh)
        ctrl_layout.addWidget(self.load_btn)

        self.plot_widget = PlotWidget(parent=self)
        main_layout.addWidget(self.plot_widget, 1)
        self.plot_widget.plotter.interactor.mousePressEvent = self.mousePressEvent
        self.plot_widget.plotter.interactor.mouseReleaseEvent = self.mouseReleaseEvent


        self.header_tree = QTreeWidget()
        self.header_tree.setHeaderLabels([" "])
        self.header_tree.setRootIsDecorated(True)
        self.header_tree.setIndentation(10)
        ctrl_layout.addWidget(self.header_tree)

        self.bc_root_item = QTreeWidgetItem(self.header_tree, ["Randbedingungen"])
        self.bc_root_item.setFont(0, QFont("", weight=QFont.Weight.Bold))
        self.bc_root_item.setBackground(0, QBrush(QColor(200, 200, 250)))
        self.bc_root_item.setExpanded(True)
        self.header_tree.addTopLevelItem(self.bc_root_item)

        self.add_bc_btn = QPushButton("+ Neue Randbedingung")
        self.add_bc_btn.clicked.connect(self.add_boundary_condition)
        self.header_tree.setItemWidget(QTreeWidgetItem(self.bc_root_item, [""]), 0, self.add_bc_btn)

        self.header_tree.setContextMenuPolicy(Qt.ContextMenuPolicy.CustomContextMenu)
        self.header_tree.customContextMenuRequested.connect(self.show_context_menu)
        ctrl_layout.addWidget(self.header_tree)


        mat_item = QTreeWidgetItem(self.header_tree, ["Materialdefinition"])
        mat_item.setFont(0, QFont("", weight=QFont.Weight.Bold))
        mat_item.setBackground(0, QBrush(QColor(200, 200, 250)))
        mat_item.setExpanded(True)

        # --- Dropdowns für Materialauswahl ---
        sel_mat_item = QTreeWidgetItem(mat_item)
        self.mat_dir_model = ["isotrop", "orthotrop", "transversal isotrop"]
        self.mat_combo = QComboBox()
        self.mat_combo.addItems(self.mat_dir_model)
        self.mat_combo.setCurrentIndex(0)
        self.header_tree.setItemWidget(sel_mat_item, 0, self.mat_combo)

        sel_mat_item2 = QTreeWidgetItem(mat_item)
        self.mat_model = ["linear elastisch", "hyperelastisch", "linear elastisch ideal plastisch"]
        self.mat_combo2 = QComboBox()
        self.mat_combo2.addItems(self.mat_model)
        self.mat_combo2.setCurrentIndex(0)
        self.header_tree.setItemWidget(sel_mat_item2, 0, self.mat_combo2)

        # --- Parameterfeld-Container ---
        self.mat_params_widget = QWidget()
        self.mat_params_layout = QVBoxLayout(self.mat_params_widget)
        self.mat_params_layout.setContentsMargins(2, 2, 2, 2)
        self.mat_params_widget.setSizePolicy(QSizePolicy.Policy.Expanding, QSizePolicy.Policy.Fixed)
        self.mat_params_widget.setMinimumHeight(120)

        self.mat_params_tree_item = QTreeWidgetItem(mat_item, [""])  # leer, kein "Parameter"-Label
        self.header_tree.setItemWidget(self.mat_params_tree_item, 0, self.mat_params_widget)
        ctrl_layout.addWidget(self.header_tree)
        
        self.result_root_item = QTreeWidgetItem(self.header_tree, ["Ergebnisse"])
        self.result_root_item.setFont(0, QFont("", weight=QFont.Weight.Bold))
        self.result_root_item.setBackground(0, QBrush(QColor(200, 200, 250)))
        self.result_root_item.setExpanded(True)
        self.header_tree.addTopLevelItem(self.result_root_item)

        self.add_result_btn = QPushButton("+ Neuer Analyseblock")
        self.add_result_btn.clicked.connect(self.add_result_evaluation)
        self.header_tree.setItemWidget(QTreeWidgetItem(self.result_root_item, [""]), 0, self.add_result_btn)
        self.header_tree.itemClicked.connect(self.show_results)




        # Play Button Berechnung
        self.playbutton = QPushButton("Berechnung durchführen")
        self.playbutton.clicked.connect(lambda: aufrufauswahl())
        ctrl_layout.addWidget(self.playbutton)



        # --- Update-Funktion für Material ---
        def update_material():
            self.define_material(self.mat_params_layout)
            # WICHTIG: Widget erneut zuweisen, damit neue Felder sofort angezeigt werden
            self.header_tree.setItemWidget(self.mat_params_tree_item, 0, self.mat_params_widget)
            #self.mat_root_item.setExpanded(True)
            mat_item.setExpanded(True)

        self.mat_combo.currentIndexChanged.connect(lambda _: update_material())
        self.mat_combo2.currentIndexChanged.connect(lambda _: update_material())

        # Initiales Aufbauen
        update_material()
        

        # definition welche Berechnungen ausgeführt werden - statische FEM, Modalanalyse, transiente Analyse... 
        # (momentan auch unterteilung anhand materialmodell das aber eher im Materialmodell selber differenzieren)
        def aufrufauswahl():
            #hier wird über if schleifen festgestellt welche methode des controllers aufgerufen werden soll
            if self.mat_combo2.currentText()=="linear elastisch":
                self.controller.static_structural()
                print("lineare FEM wird ausgeführt")
            elif self.mat_combo2.currentText()=="hyperelastisch":
                print("hyperelastische Berechnung wird ausgeführt")



       


    def show_context_menu(self,pos):
        item = self.header_tree.itemAt(pos)
        bc_index = item.data(0, Qt.ItemDataRole.UserRole)
        # Wurde nicht auf ein Item geklickt?
        if item is None:
            return
        if item.parent() != self.bc_root_item:
            return
        menu = QMenu(self)
        delete_action = menu.addAction("Randbedingung löschen")
        action = menu.exec(self.header_tree.viewport().mapToGlobal(pos))
        if action == delete_action:
            self.delete_bc(bc_index)


    def get_face_nodes(self, surf_cell_idx):
        element_id, local_face_id = self.surface_face_mapping[surf_cell_idx]

        cell_type, pts_ids = self.elements[element_id]

        face_nodes_global = pts_ids[ELEMENT_FACES[cell_type][local_face_id]]
        return face_nodes_global

    # ------------------------ MOUSE EVENTS ------------------------
    def mousePressEvent(self, event):
        if event.button() == Qt.MouseButton.LeftButton:
            self.mouse_pressed_pos = event.pos()
        QtInteractor.mousePressEvent(self.plot_widget.plotter.interactor, event)

    def mouseReleaseEvent(self, event):
        if event.button() == Qt.MouseButton.LeftButton:
            if self.mouse_pressed_pos is not None:
                distance = (event.pos() - self.mouse_pressed_pos).manhattanLength()
                if distance < self.drag_threshold:
                    click_point = self.plot_widget.plotter.pick_mouse_position()
                    if click_point is not None:
                        if self.selection_mode == "Knoten":
                            self._on_point_picked(click_point, None, self.mesh)
                        else:
                            self._on_face_clicked(click_point)
        self.mouse_pressed_pos = None
        QtInteractor.mouseReleaseEvent(self.plot_widget.plotter.interactor, event)

    # ------------------------ MESH LADEN und DATEN ZUWEISEN------------------------
    def on_load_mesh(self):
        path, _ = QFileDialog.getOpenFileName(self, "MSH Datei auswählen", "", "MSH Files (*.msh)")
        if not path:
            return
        try:
            self.meshio_mesh = meshio.read(path)                                          ###### self
            self.mesh, self.elements = self._meshio_to_pyvista(self.meshio_mesh)

        except Exception as e:
            QMessageBox.critical(self, "Fehler", f"Mesh konnte nicht gelesen werden:\n{e}")
            return

        self.mesh_points    = np.asarray(self.mesh.points)
        self.n_elements     = 0 #müsste oben im initiator stehen
        

        #ÜBERLEGEN OB DIESER PART NICHT EHER IN NODEMATRICES IN DATAPREPARATION GEHÖRT HAT JA NICHTS MEHR MIT DEM EINLESEN ZU TUN 
        for i in range(0, len(self.meshio_mesh.cells)):
            if self.meshio_mesh.cells[i].type in ('hexahedron','tetrahedron','wedgde'):   #funktioniert nur mit 3d für 2d wäre separate schleife nötig
                if self.mesh_cells is None:
                    self.mesh_cells = self.meshio_mesh.cells[i].data
                    self.el_type_num.append(len(self.meshio_mesh.cells[i].data))

                    if self.meshio_mesh.cells[i].type == 'hexahedron':
                        if np.shape(self.meshio_mesh.cells[i].data)[1]==8:
                           self.el_type.append('hex8')
                           self.n_SF.append(8)
                           self.n_el_node.append(8)
                        elif np.shape(self.meshio_mesh.cells[i].data)[1]==27: 
                           self.el_type.append('hex27')
                           self.n_SF.append(27)
                           self.n_el_node.append(27)

                    elif self.meshio_mesh.cells[i].type == 'tetrahedron':
                        if  np.shape(self.meshio_mesh.cells[i].data)[1]==4:
                            self.el_type.append('tet4')
                            self.n_SF.append(4)
                            self.n_el_node.append(4)
                        elif  np.shape(self.meshio_mesh.cells[i].data)[1]==10:
                            self.el_type.append('tet10??')
                            self.n_SF.append(10)
                            self.n_el_node.append(10)


                    elif self.meshio_mesh.cells[i].type == 'wedge':
                        print(np.shape(self.meshio_mesh.cells[i].data)[1])
                        print(self.meshio_mesh.cells[i].type)
                        print("Es wurde im Meshing ein Element verwendet das in dieser FEM Software nicht unterstützt wird")
                        print(i)
                else: 
                    if np.shape(self.meshio_mesh.cells[i].data)[1] > np.shape(self.mesh_cells)[1]:
                        self.mesh_cells = np.hstack(self.mesh_cells, np.zeros((np.shape(self.mesh_cells)[0],np.shape(self.meshio_mesh.cells[i].data)[1]-np.shape(self.mesh_cells)[1])))
                    elif np.shape(self.mesh_cells)[1] > np.shape(self.meshio_mesh.cells[i].data)[1]:
                        self.mesh_cells = np.hstack(self.mesh_cells, np.zeros((np.shape(self.mesh_cells)[0],np.shape(self.mesh_cells)[1]-np.shape(self.meshio_mesh.cells[i].data)[1])))
                    self.mesh_cells = np.append(self.mesh_cells,self.meshio_mesh.cells[i].data)
                
                self.n_elements += len(self.meshio_mesh.cells[i].data)  
              





        
        self.nodes_polydata = pv.PolyData(self.mesh.points)

       
        self._prepare_surface_faces()

        # Reset selections
        self._clear_plotter_keep_plotter()
        self.selected_nodes.clear()
        self.selected_faces.clear()
        self.node_selection_actors.clear()
        self.face_selection_actors.clear()

        self.plot_widget.plotter.add_mesh(self.mesh, color="lightgrey", show_edges=True, opacity=1.0, name="volume")
        self.plot_widget.plotter.add_mesh(self.nodes_polydata, color="darkgrey", point_size=12, render_points_as_spheres=True, name="nodes_points")
        self.plot_widget.plotter.reset_camera()
        self.plot_widget.plotter.render()
        self._setup_picking()

    # Meshio to PyVista
    def _meshio_to_pyvista(self, meshio_mesh):
        points = meshio_mesh.points
        all_cells = []
        all_cell_types = []
        elements = []

        
        

        for cell_block in meshio_mesh.cells:
            #print("cell_block.data",cell_block.data)
            ctype = cell_block.type
            if ctype not in VTK_CELL_TYPES:
                continue
            cells = np.asarray(cell_block.data, dtype=int)
            if cells.size == 0:
                continue
            n_points = cells.shape[1]

            #print("n_points",n_points)
            pv_cells = np.hstack([np.full((cells.shape[0], 1), n_points, dtype=int), cells])
            all_cells.append(pv_cells)
            all_cell_types.append(np.full(cells.shape[0], VTK_CELL_TYPES[ctype], dtype=np.uint8))
            for cell in cells:
                elements.append((VTK_CELL_TYPES[ctype], cell))

        if not all_cells:
            raise ValueError("Keine unterstützten 3D-Elementtypen gefunden.")

        cells_concat = np.vstack(all_cells).ravel()
        cell_types_concat = np.concatenate(all_cell_types)
        return pv.UnstructuredGrid(cells_concat, cell_types_concat, points), elements








    # ------------------------ OBERFLÄCHEN ------------------------
    def _prepare_surface_faces(self):
        self.surface_faces = pv.PolyData()
        self.surface_face_mapping.clear()
        surf_cells = []
        cell_counter = 0

        for elem_id, (cell_type, pts_ids) in enumerate(self.elements):
            if cell_type not in ELEMENT_FACES:
                continue
            faces = ELEMENT_FACES[cell_type]

            for local_face_id, face_nodes in enumerate(faces):
                face_pts = pts_ids[face_nodes]
                surf_cells.append([len(face_pts)] + list(face_pts))
                self.surface_face_mapping[cell_counter] = (elem_id, local_face_id)
                cell_counter += 1
        if surf_cells:
            cells_np = np.array(surf_cells).ravel()
            self.surface_faces = pv.PolyData(self.mesh.points, cells_np)



    







    # ------------------------ PICKING ------------------------
    def _setup_picking(self):
        self._disable_picking()

    def _disable_picking(self):
        try:
            self.plot_widget.plotter.disable_picking()
        except Exception:
            pass

    # ------------------------ RANDBEDINGUNGEN ------------------------
    def add_boundary_condition(self):
        items = ["Kraftrandbedingung", "Verschiebungsrandbedingung", "Lager"]
        item, ok = QInputDialog.getItem(self, "Randbedingung auswählen", "Typ der Randbedingung:", items, 0, False)
        if not (ok and item):
            return


        #facenodes evtl nicht notwendig da eh alle nodes des elements a dem die fläche berechnet wird gebraucht werden 
        bc = {"type": item, "nodes": set(), "faces": set(), "elementtype":set(),"facenodes":set(),"elementnodes":set(),
              "values": [0,0,0] if item != "Lager" else None,
              "tree_items": {},"status":False}
        self.boundary_conditions.append(bc)
        index = len(self.boundary_conditions) - 1
        self.current_bc_index = index


        self.bc_item = QTreeWidgetItem(self.bc_root_item, [item])
        self.bc_item.setData(0, Qt.ItemDataRole.UserRole, index)
        self.bc_item.setFont(0, QFont("", weight=QFont.Weight.Bold))
        self.bc_item.setBackground(0, QBrush(QColor(255, 0, 0)))
        self.bc_item.setExpanded(True)



        

        # Dropdown Auswahlmodus
        #sel_item = QTreeWidgetItem(bc_item, ["Auswahlmodus"])
        #combo = QComboBox()
        #combo.addItems(["Knoten", "Fläche"])
        #combo.currentIndexChanged.connect(lambda idx, bc_idx=index: self.set_bc_selection_mode(bc_idx, idx))
        #self.bc_header_tree.setItemWidget(sel_item, 0, combo)

        # Dropdown Auswahlmodus2
        #sel_item2 = QTreeWidgetItem(bc_item, ["Auswahlmodus"])
        #combo2 = QComboBox()
        if item == "Kraftrandbedingung":
            combo2 = QComboBox()
            combo2.addItems(["Kraft", "Drehmoment"])
            sel_item2 = QTreeWidgetItem(self.bc_item)
            self.header_tree.setItemWidget(sel_item2, 0, combo2)

            # Dropdown Auswahlmodus
            combo = QComboBox()
            combo.addItems(["Knoten", "Fläche"])
            sel_item = QTreeWidgetItem(self.bc_item)
            combo.currentIndexChanged.connect(lambda idx, bc_idx=index: self.set_bc_selection_mode(bc_idx, idx))
            self.header_tree.setItemWidget(sel_item, 0, combo)
        #else:
        #    combo2.addItems(["Translation", "Rotation"])
        #self.bc_header_tree.setItemWidget(sel_item2, 0, combo2)
        elif item=="Verschiebungsrandbedingung":
            combo = QComboBox()
            combo.addItems(["Translation", "Rotation"])
            sel_item2 = QTreeWidgetItem(self.bc_item)
            self.header_tree.setItemWidget(sel_item2, 0, combo)

        # Zuordnung
        list_item = QTreeWidgetItem(self.bc_item)
        list_widget = QListWidget()
        list_widget.setFixedHeight(100)
        self.header_tree.setItemWidget(list_item, 0, list_widget)

        # Werte / Locks
        val_item = QTreeWidgetItem(self.bc_item, [""])
        val_widget = QWidget()
        layout = QVBoxLayout(val_widget)
        layout.setSpacing(2)
        layout.setContentsMargins(2, 2, 2, 2)

        if item == "Lager":

            # Dropdown Auswahlmodus
            combo = QComboBox()
            combo.addItems(["Feste Einspannung","Festlager","Loslager","Gelenklager"])
            sel_item = QTreeWidgetItem(self.bc_item)
            self.header_tree.setItemWidget(sel_item, 0, combo)
            self.lock_checkboxes={}

            if combo.currentText()=="Feste Einspannung":
                self.set_bc_selection_mode(index, idx=1)    #kann nur Flächen auswählen (evtl noch machen das diese immer grün sind)/für gelenklager zb nur Knoten auswählbar
                self.lock_checkboxes[len(self.lock_checkboxes)-1] = {('Tx:',True),('Ty:',True),('Tz:',True),('Rx:',True),('Ry:',True),('Rz:',True)}   # indexierung noch anders da rückwirkend keine Änderungen mehr gemacht werden können only new bearings can be added correctly
                bc["tree_items"]["bearing_type"] = combo.currentText()
                print("welches lager haben wir",combo.currentText())

                
               









            #self.locks = {}
            #self.lock_checkboxes = {}
            #for dof in ["Tx","Ty","Tz","Rx","Ry","Rz"]:
            #    cb = QCheckBox(dof)
            #    layout.addWidget(cb)

            #    self.lock_checkboxes[dof] = cb
            #    self.locks[dof] = cb.isChecked()
            #    cb.stateChanged.connect(lambda _ : self.update_bc_value(bc,item))
            #bc["tree_items"]["locks"] = self.locks
            #self.x_input = self.y_input = self.z_input = None
            #cb.stateChanged.connect(lambda _ : self.update_bc_value(bc,item))
            
        else:


            if item == "Kraftrandbedingung":
                Qlabel =["Fx:","Fy:","Fz:"]
            elif item == "Verschiebungsrandbedingung":
                Qlabel =["ux:","uy:","uz:"]


            self.x_input = QLineEdit("0"); layout.addWidget(QLabel(Qlabel[0])); layout.addWidget(self.x_input)
            self.y_input = QLineEdit("0"); layout.addWidget(QLabel(Qlabel[1])); layout.addWidget(self.y_input)
            self.z_input = QLineEdit("0"); layout.addWidget(QLabel(Qlabel[2])); layout.addWidget(self.z_input)

            self.x_input.textChanged.connect(lambda _ : self.update_arrow())
            self.y_input.textChanged.connect(lambda _ : self.update_arrow())
            self.z_input.textChanged.connect(lambda _ : self.update_arrow())


            self.x_input.editingFinished.connect(lambda : self.update_bc_value(bc,item))
            self.y_input.editingFinished.connect(lambda : self.update_bc_value(bc,item))
            self.z_input.editingFinished.connect(lambda : self.update_bc_value(bc,item))

        
        self.header_tree.setItemWidget(val_item, 0, val_widget)

        # store references
        #bc["tree_items"]["sel"]             = combo
        bc["tree_items"]["list"]            = list_widget
        bc["tree_items"]["tree_item_ref"]   = self.bc_item



    def set_bc_selection_mode(self, bc_index, idx):
        self.current_bc_index = bc_index
        self.selection_mode = "Knoten" if idx == 0 else "Fläche"




    def delete_bc(self, bc_index):

        bc = self.boundary_conditions[bc_index]

        # Ausgewählte Flächen entfernen
        for surf_cell_idx in list(bc["faces"]):

            # face Actor entfernen
            if surf_cell_idx in self.face_selection_actors:
                self.plot_widget.plotter.remove_actor(
                    self.face_selection_actors[surf_cell_idx]
                )
                del self.face_selection_actors[surf_cell_idx]

            # Fläche aus den global ausgewählten Flächen entfernen
            self.selected_faces.discard(surf_cell_idx)
        bc["faces"].clear()

        # arrow Actor entfernen
        if bc_index < len(self.arrow_actor):
            try:
                self.plot_widget.plotter.remove_actor(self.arrow_actor[bc_index])
                self.plot_widget.plotter.render()
            except Exception:
                pass

        # tree item entfernen
        bc_item = bc["tree_items"]["tree_item_ref"]
        parent_item = bc_item.parent()
        parent_item.removeChild(bc_item)

        # entfernen aus boundary_conditions, arrow_actor und arrow_nodes liste
        self.boundary_conditions.pop(bc_index)
        try:
            self.arrow_actor.pop(bc_index)
            self.arrow_nodes.pop(bc_index)
        except Exception:
            pass 


        # Update indices in tree items
        for i, bc in enumerate(self.boundary_conditions):
            bc["tree_items"]["tree_item_ref"].setData(0, Qt.ItemDataRole.UserRole, i)

        #  current_bc_index zurücksetzen 
        if self.current_bc_index == bc_index:
            self.current_bc_index = None
        

    # Wert wird automatisch durch connect geupdatet
    def update_bc_value(self,bc,item): 

        if item != 'Lager': 
            bc["values"] = [float(self.x_input.text()),float(self.y_input.text()),float(self.z_input.text())]
        else:
            for dof, cb in self.lock_checkboxes.items():
                self.locks[dof] = cb.isChecked()
                bc["tree_items"]["locks"] = self.locks

        self.check_widget_status(bc)

        


    def check_widget_status(self, bc):
        if bc["faces"] and bc["type"]!="Lager":
            if bc["values"]!=[0,0,0]:
                self.bc_item.setBackground(0, QBrush(QColor(0, 255, 0)))
                bc["status"]=True
            else:
                self.bc_item.setBackground(0, QBrush(QColor(255, 0, 0)))
                bc["status"]=False
        elif bc["faces"]:
            self.bc_item.setBackground(0, QBrush(QColor(0, 255, 0)))
            bc["status"]=True
        else:
            self.bc_item.setBackground(0, QBrush(QColor(255, 0, 0)))
            bc["status"]=False






    def update_arrow(self):

        # Leeren String abfangen
        try:
            x_val = float(self.x_input.text())
        except ValueError:
            x_val = 0.0
        try:
            y_val = float(self.y_input.text())
        except ValueError:
            y_val = 0.0
        try:
            z_val = float(self.z_input.text())
        except ValueError:
            z_val = 0.0

        # Nullpunkt des Pfeils ermitteln
        try:
            coordinate_data = np.zeros((len(self.arrow_nodes[self.current_bc_index]["nodes"]),3))
            maxmin = np.zeros((3,2))
            for i in range(0,len(self.arrow_nodes[self.current_bc_index]["nodes"])):
                coordinate_data[i,:] = self.mesh_points[self.arrow_nodes[self.current_bc_index]["nodes"][i]]  #achtung noch korrektur vornehmen da sich evtl lagerbedingungen dazwischen mischen können

            maxmin[0,0] = np.max(coordinate_data[:,0])
            maxmin[1,0] = np.max(coordinate_data[:,1])
            maxmin[2,0] = np.max(coordinate_data[:,2])
            maxmin[0,1] = np.min(coordinate_data[:,0])
            maxmin[1,1] = np.min(coordinate_data[:,1])
            maxmin[2,1] = np.min(coordinate_data[:,2])

            x_start = (maxmin[0,0] + maxmin[0,1])/2
            y_start = (maxmin[1,0] + maxmin[1,1])/2
            z_start = (maxmin[2,0] + maxmin[2,1])/2
        except:
            return
       
        # Pfeil löschen damit sich der Pfeil bei Änderung der Werte neu ausrichtet
        try:
            self.arrow_actor[self.current_bc_index] 
            self.plot_widget.plotter.remove_actor(self.arrow_actor[self.current_bc_index])
            self.plot_widget.plotter.render()
        except:
            pass 
        
        # Pfeil erstellen
        try:
            self.arrow_actor[self.current_bc_index] = self.plot_widget.plotter.add_arrows(
                np.array([[x_start, y_start, z_start]]),
                np.array([[x_val, y_val, z_val]]),
                mag=10.0,
                color="black"
            )
        except:
            self.arrow_actor.append(self.plot_widget.plotter.add_arrows(
                np.array([[x_start, y_start, z_start]]),
                np.array([[x_val, y_val, z_val]]),
                mag=10.0,
                color="black"
            ))

        self.plot_widget.plotter.render()

        

    # ------------------------ MATERIALDEFINITION ------------------------
    def define_material(self,layout):
        # Alte Widgets löschen
        self.clear_layout(layout)

        E_Modul = []
        v = []

        if self.mat_combo.currentText() == "isotrop" and self.mat_combo2.currentText() == "linear elastisch":
            E_Modul.append(QLineEdit("0"))
            v.append(QLineEdit("0"))
            layout.addWidget(QLabel("E Modul:"))
            layout.addWidget(E_Modul[0])
            layout.addWidget(QLabel("Querkontraktionszahl:"))
            layout.addWidget(v[0])
        elif self.mat_combo.currentText() == "orthotrop" and self.mat_combo2.currentText() == "linear elastisch":
            left_layout = QVBoxLayout()
            right_layout = QVBoxLayout()
            row_layout = QHBoxLayout()
            for i, label in enumerate(["x","y","z"]):
                e = QLineEdit("0"); E_Modul.append(e); left_layout.addWidget(QLabel(f"E Modul {label}:")); left_layout.addWidget(e)
                vi = QLineEdit("0"); v.append(vi); right_layout.addWidget(QLabel(f"Querkontraktionszahl {label}:")); right_layout.addWidget(vi)
            row_layout.addLayout(left_layout)
            row_layout.addLayout(right_layout)
            layout.addLayout(row_layout)
        else:
            E_Modul.append(QLineEdit("0"))
            layout.addWidget(QLabel("E Modul (allgemein):"))
            layout.addWidget(E_Modul[0])

        # Werte speichern
        self.mat_values = {"E": E_Modul, "v": v}

        # **Neu anhängen + sofort sichtbar machen**
        self.mat_params_widget.adjustSize()
        self.mat_params_widget.updateGeometry()
        self.mat_params_widget.update()
        self.header_tree.setItemWidget(self.mat_params_tree_item, 0, self.mat_params_widget)  # <- wichtig
        self.header_tree.viewport().update()




    # ------------------------ AUSWERTUNG------------------------
    def add_result_evaluation(self):
        items = ["x Verschiebung","y Verschiebung","z Verschiebung", "van Mises Vergleichsspannung", "Hauptnormalspnnungen"]
        item, ok = QInputDialog.getItem(self, "Ergebnis auswählen", "Typ des Ergebnisses:", items, 0, False)
        if not (ok and item):
            return
        

        result = {"type": item}
        self.evaluated_result.append(result)
        
        result_item = QTreeWidgetItem(self.result_root_item, [item])
        result_item.setFont(0, QFont("", weight=QFont.Weight.Bold))
        result_item.setBackground(0, QBrush(QColor(220, 220, 220)))
        result_item.setExpanded(True)
        
       


    def clear_layout(self, layout):
        if layout is None:
            return
        while layout.count():
            item = layout.takeAt(0)
            if item is None:
                continue
            widget = item.widget()
            child_layout = item.layout()
            if widget:
                widget.deleteLater()
            elif child_layout:
                self.clear_layout(child_layout)

    # -------------------- BERECHNUNG STARTEN -------------------
    def starte_berechnung(self):
        if self.mat_model =="linear elastisch":
            self.controller.lineare_FEM()
        

    # ------------------------ SELEKTION ------------------------
    def _on_point_picked(self, point, actor, mesh):
        if self.mesh is None or point is None or self.current_bc_index is None:
            return
        pid = int(self.mesh.find_closest_point(np.array(point)))
        if pid in self.selected_nodes:
            self._deselect_node(pid)
        else:
            self._select_node(pid)

    def _on_face_clicked(self, click_point):
        if self.surface_faces is None or click_point is None or self.current_bc_index is None:
            return
        pt = np.array(click_point)
        centers = self.surface_faces.cell_centers().points
        if centers.size == 0:
            return
        dists = np.linalg.norm(centers - pt.reshape(1,3), axis=1)
        surf_cell_idx = int(np.argmin(dists))
        if surf_cell_idx in self.selected_faces:
            self._deselect_face(surf_cell_idx)
        else:
            self._select_face(surf_cell_idx)


    # Node select/deselect
    def _select_node(self, pid):
        if pid < 0 or pid >= len(self.mesh.points):
            return
        self.selected_nodes.add(pid)
        coord = self.mesh.points[pid]
        radius = max(self.mesh.length/50.0,1e-6)
        sphere = pv.Sphere(radius=radius, center=coord, theta_resolution=20, phi_resolution=20)
        actor = self.plot_widget.plotter.add_mesh(sphere, color="red", opacity=1.0)
        self.node_selection_actors[pid] = actor
        if self.current_bc_index is not None:
            bc = self.boundary_conditions[self.current_bc_index]
            bc["nodes"].add(pid)
            list_widget = bc["tree_items"]["list"]
            list_widget.addItem(f"Knoten {pid}: {coord[0]:.3f},{coord[1]:.3f},{coord[2]:.3f}")
        self.plot_widget.plotter.render()

    def _deselect_node(self, pid):
        if pid not in self.selected_nodes:
            return
        self.selected_nodes.remove(pid)
        actor = self.node_selection_actors.pop(pid, None)
        if actor:
            try: self.plot_widget.plotter.remove_actor(actor)
            except: pass
        if self.current_bc_index is not None:
            bc = self.boundary_conditions[self.current_bc_index]
            bc["nodes"].discard(pid)
            list_widget = bc["tree_items"]["list"]
            for i in range(list_widget.count()):
                if list_widget.item(i).text().startswith(f"Knoten {pid}:"):
                    list_widget.takeItem(i)
                    break
        self.plot_widget.plotter.render()


    # Face select/deselect
    def _select_face(self, surf_cell_idx):
        if surf_cell_idx in self.selected_faces:
            self._deselect_face(surf_cell_idx)
            return
        self.selected_faces.add(surf_cell_idx)
        element_id, local_face_id = self.surface_face_mapping.get(surf_cell_idx, (None, None))
        bc = self.boundary_conditions[self.current_bc_index]
        bc_type = bc["type"]

        try:
            face_poly = self.surface_faces.extract_cells(surf_cell_idx).copy()
            offset_small = 1e-6 * max(1.0, getattr(self.mesh,"length",1.0))
            face_poly.points[:] += offset_small
            color = "green" if bc_type=="Lager" else "red"
            actor = self.plot_widget.plotter.add_mesh(face_poly, color=color, opacity=1.0, reset_camera=False, show_edges=False)
            try:
                prop = actor.GetProperty()
                prop.SetPolygonOffset(1.0,1.0)
                try: prop.PolygonOffsetOn()
                except: pass
            except: pass
            self.face_selection_actors[surf_cell_idx] = actor
            bc["faces"].add(surf_cell_idx)
 
            list_widget = bc["tree_items"]["list"]
            list_widget.addItem(f"Fläche {element_id},{local_face_id}")
            self.plot_widget.plotter.render()

        except Exception as e:
            print("Fehler beim Hervorheben der Fläche:", e)

        self.check_widget_status(bc)







        for surf_idx in bc["faces"]:
            faces = surf_idx
            element_id = self.surface_face_mapping[faces][0]
            element_type = self.elements[element_id][0]
            self.boundary_conditions[self.current_bc_index]["elementtype"].add(element_type)
            self.boundary_conditions[self.current_bc_index]["elementnodes"].add(tuple(int(i) for i in self.mesh_cells[element_id]))

        #for i in range(0,len(self.boundary_conditions)):

        if self.boundary_conditions[self.current_bc_index]["type"]!="Lager":
            all_nodes2 = np.array([],dtype=int)
            for surf_idx in bc["faces"]:
                self.boundary_conditions[self.current_bc_index]["facenodes"].add(str(self.get_face_nodes(surf_idx)))
                nodes = np.array(self.get_face_nodes(surf_idx))
                all_nodes2 = np.append(all_nodes2,nodes)  # besser als list formulieren da np jedes mal neu aufsetzen musss mit .append 
            try:
                self.arrow_nodes[self.current_bc_index]["nodes"] = all_nodes2
            except:
                self.arrow_nodes.append({})
                self.arrow_nodes[self.current_bc_index]["nodes"] = all_nodes2

            self.update_arrow()
                
        else:
            for surf_idx in bc["faces"]:
                self.boundary_conditions[self.current_bc_index]["facenodes"].add(str(self.get_face_nodes(surf_idx)))
                nodes = np.array(self.get_face_nodes(surf_idx)[0])
            try:
                self.arrow_nodes[self.current_bc_index]["nodes"] = None
            except:
                self.arrow_nodes.append({})
                self.arrow_nodes[self.current_bc_index]["nodes"] = None


            



    def _deselect_face(self, surf_cell_idx):
        if surf_cell_idx not in self.selected_faces:
            return
        self.selected_faces.remove(surf_cell_idx)
        actor = self.face_selection_actors.pop(surf_cell_idx,None)
        if actor:
            try: self.plot_widget.plotter.remove_actor(actor)
            except: pass
        if self.current_bc_index is not None:
            bc = self.boundary_conditions[self.current_bc_index]
            bc["faces"].discard(surf_cell_idx)
            list_widget = bc["tree_items"]["list"]
            elem_local = f"{self.surface_face_mapping[surf_cell_idx][0]},{self.surface_face_mapping[surf_cell_idx][1]}"
            for i in range(list_widget.count()):
                if elem_local in list_widget.item(i).text():
                    list_widget.takeItem(i)
                    break
        self.update_arrow()
        self.check_widget_status(bc)




        for surf_idx in bc["faces"]:
            faces = surf_idx
            element_id = self.surface_face_mapping[faces][0]
            element_type = self.elements[element_id][0]
            self.boundary_conditions[self.current_bc_index]["elementtype"].add(element_type)
            self.boundary_conditions[self.current_bc_index]["elementnodes"].add(tuple(int(i) for i in self.mesh_cells[element_id]))
           
        

        #for i in range(0,len(self.boundary_conditions)):

        if self.boundary_conditions[self.current_bc_index]["type"]!="Lager":
            all_nodes2 = np.array([],dtype=int)
            for surf_idx in bc["faces"]:
                self.boundary_conditions[self.current_bc_index]["facenodes"].add(str(self.get_face_nodes(surf_idx)))
                nodes = np.array(self.get_face_nodes(surf_idx))
                all_nodes2 = np.append(all_nodes2,nodes)  # besser als list formulieren da np jedes mal neu aufsetzen musss mit .append 
            try:
                self.arrow_nodes[self.current_bc_index]["nodes"] = all_nodes2
            except:
                self.arrow_nodes.append({})
                self.arrow_nodes[self.current_bc_index]["nodes"] = all_nodes2

            self.update_arrow()
                
        else:
            for surf_idx in bc["faces"]:
                self.boundary_conditions[self.current_bc_index]["facenodes"].add(str(self.get_face_nodes(surf_idx)))
                nodes = np.array(self.get_face_nodes(surf_idx)[0])
            try:
                self.arrow_nodes[self.current_bc_index]["nodes"] = None
            except:
                self.arrow_nodes.append({})
                self.arrow_nodes[self.current_bc_index]["nodes"] = None
        






        self.update_arrow()
        self.plot_widget.plotter.render()

        

    # ------------------------ UTILITY ------------------------
    def _clear_plotter_keep_plotter(self):
        try:
            self.plot_widget.plotter.clear()
        except Exception:
            pass


    # -------------------- PLOT RESULTS ------------------------
    def show_results(self,item,column):

        print("ich habe das item ",item.text(column)," ausgewählt")

        if item.text(column)=="x Verschiebung" or item.text(column)=="y Verschiebung" or item.text(column)=="z Verschiebung" :
            print(item.text(column))
            self.plot_widget.show_deformed_mesh(self.mesh,self.controller.displacements,self.controller.displacements,item.text(column))
        elif item.text(column)=="Spannungen":
            #try:
            #    self.plot_widget.show_deformed_mesh(self.controller.displacements,self.controller.van_Mises)
            #except:
                #controller muss dann postprocessor auffordern zb. die spannungen noch nachträglich aus den verschiebungen zu berechnen    
            pass

        # je nach ausgewähltem item muss Fallunterscheidung durch if bedingung sein 



        



