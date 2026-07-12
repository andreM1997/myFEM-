# gui_fem_bc_final_material_fixed.py
import sys
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
    def __init__(self):
        super().__init__()
        self.setWindowTitle("myFEM")
        self.setGeometry(60, 60, 1800, 900)

        # Meshdaten
        self.mesh = None
        self.mesh_points = None
        self.nodes_polydata = None
        self.surface_faces = None
        self.surface_face_mapping = {}
        self.elements = []

        # Selektion
        self.selected_nodes = set()
        self.selected_faces = set()
        self.node_selection_actors = {}
        self.face_selection_actors = {}

        # Randbedingungen
        self.boundary_conditions = []
        self.current_bc_index = None
        self.selection_mode = "Knoten"

        # Drag detection
        self.mouse_pressed_pos = None
        self.drag_threshold = 5  # Pixel

        # Materialwerte
        self.mat_values = {"E": [], "v": []}

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

        # -----------------------
        # Randbedingungen Tree
        # -----------------------
        self.bc_header_tree = QTreeWidget()
        self.bc_header_tree.setHeaderLabels(["Randbedingungen"])
        self.bc_header_tree.setColumnCount(1)
        self.bc_header_tree.setRootIsDecorated(True)
        self.bc_header_tree.setIndentation(10)
        ctrl_layout.addWidget(self.bc_header_tree)

        self.bc_root_item = QTreeWidgetItem(self.bc_header_tree, ["Randbedingungen"])
        self.bc_root_item.setFont(0, QFont("", weight=QFont.Weight.Bold))
        self.bc_root_item.setBackground(0, QBrush(QColor(200, 200, 250)))
        self.bc_root_item.setExpanded(True)
        self.bc_header_tree.addTopLevelItem(self.bc_root_item)

        self.add_bc_btn = QPushButton("+ Neue Randbedingung")
        self.add_bc_btn.clicked.connect(self.add_boundary_condition)
        self.bc_header_tree.setItemWidget(QTreeWidgetItem(self.bc_root_item, [""]), 0, self.add_bc_btn)

    # -----------------------
    # Material Tree direkt unter Randbedingungen
    # -----------------------
        self.mat_header_tree = QTreeWidget()
        self.mat_header_tree.setHeaderLabels(["Material"])
        self.mat_header_tree.setColumnCount(1)
        self.mat_header_tree.setRootIsDecorated(True)
        self.mat_header_tree.setIndentation(10)
        ctrl_layout.addWidget(self.mat_header_tree)

        self.mat_root_item = QTreeWidgetItem(self.mat_header_tree, ["Material"])
        self.mat_root_item.setFont(0, QFont("", weight=QFont.Weight.Bold))
        self.mat_root_item.setBackground(0, QBrush(QColor(200, 200, 250)))
        self.mat_root_item.setExpanded(True)
        self.mat_header_tree.addTopLevelItem(self.mat_root_item)

        mat_item = QTreeWidgetItem(self.mat_root_item, ["Materialdefinition"])
        mat_item.setFont(0, QFont("", weight=QFont.Weight.Bold))
        mat_item.setBackground(0, QBrush(QColor(220, 220, 220)))
        mat_item.setExpanded(True)

    # --- Dropdowns für Materialauswahl ---
        sel_mat_item = QTreeWidgetItem(mat_item, ["Richtung"])
        mat_dir_model = ["isotrop", "orthotrop", "transversal isotrop"]
        self.mat_combo = QComboBox()
        self.mat_combo.addItems(mat_dir_model)
        self.mat_combo.setCurrentIndex(0)
        self.mat_header_tree.setItemWidget(sel_mat_item, 0, self.mat_combo)

        sel_mat_item2 = QTreeWidgetItem(mat_item, ["Modell"])
        mat_model = ["linear elastisch", "hyperelastisch", "linear elastisch ideal plastisch"]
        self.mat_combo2 = QComboBox()
        self.mat_combo2.addItems(mat_model)
        self.mat_combo2.setCurrentIndex(0)
        self.mat_header_tree.setItemWidget(sel_mat_item2, 0, self.mat_combo2)

    # --- Parameterfeld-Container ---
        self.mat_params_widget = QWidget()
        self.mat_params_layout = QVBoxLayout(self.mat_params_widget)
        self.mat_params_layout.setContentsMargins(2, 2, 2, 2)
        self.mat_params_widget.setSizePolicy(QSizePolicy.Policy.Expanding, QSizePolicy.Policy.Fixed)
        self.mat_params_widget.setMinimumHeight(120)

        self.mat_params_tree_item = QTreeWidgetItem(mat_item, [""])  # leer, kein "Parameter"-Label
        self.mat_header_tree.setItemWidget(self.mat_params_tree_item, 0, self.mat_params_widget)

    # --- Update-Funktion für Material ---
        def update_material():
            direction = self.mat_combo.currentText()
            model = self.mat_combo2.currentText()
            self.define_material(direction, model, self.mat_params_layout)
            # WICHTIG: Widget erneut zuweisen, damit neue Felder sofort angezeigt werden
            self.mat_header_tree.setItemWidget(self.mat_params_tree_item, 0, self.mat_params_widget)
            self.mat_root_item.setExpanded(True)
            mat_item.setExpanded(True)

        self.mat_combo.currentIndexChanged.connect(lambda _: update_material())
        self.mat_combo2.currentIndexChanged.connect(lambda _: update_material())

        # Initiales Aufbauen
        update_material()

        # ------------------------ PLOTTER ------------------------
        plot_layout = QVBoxLayout()
        main_layout.addLayout(plot_layout, 1)
        self.plotter = QtInteractor(self)
        plot_layout.addWidget(self.plotter.interactor)
        self.plotter.set_background("white")
        self.plotter.add_axes()
        self.plotter.show_bounds(grid='front', location='outer')
        self.plotter.interactor.mousePressEvent = self.mousePressEvent
        self.plotter.interactor.mouseReleaseEvent = self.mouseReleaseEvent

    # ------------------------ MOUSE EVENTS ------------------------
    def mousePressEvent(self, event):
        if event.button() == Qt.MouseButton.LeftButton:
            self.mouse_pressed_pos = event.pos()
        QtInteractor.mousePressEvent(self.plotter.interactor, event)

    def mouseReleaseEvent(self, event):
        if event.button() == Qt.MouseButton.LeftButton:
            if self.mouse_pressed_pos is not None:
                distance = (event.pos() - self.mouse_pressed_pos).manhattanLength()
                if distance < self.drag_threshold:
                    click_point = self.plotter.pick_mouse_position()
                    if click_point is not None:
                        if self.selection_mode == "Knoten":
                            self._on_point_picked(click_point, None, self.mesh)
                        else:
                            self._on_face_clicked(click_point)
        self.mouse_pressed_pos = None
        QtInteractor.mouseReleaseEvent(self.plotter.interactor, event)

    # ------------------------ MESH LADEN ------------------------
    def on_load_mesh(self):
        path, _ = QFileDialog.getOpenFileName(self, "MSH Datei auswählen", "", "MSH Files (*.msh)")
        if not path:
            return
        try:
            meshio_mesh = meshio.read(path)
            self.mesh, self.elements = self._meshio_to_pyvista(meshio_mesh)
        except Exception as e:
            QMessageBox.critical(self, "Fehler", f"Mesh konnte nicht gelesen werden:\n{e}")
            return

        self.mesh_points = np.asarray(self.mesh.points)
        self.nodes_polydata = pv.PolyData(self.mesh.points)
        self._prepare_surface_faces()

        # Reset selections
        self._clear_plotter_keep_plotter()
        self.selected_nodes.clear()
        self.selected_faces.clear()
        self.node_selection_actors.clear()
        self.face_selection_actors.clear()

        self.plotter.add_mesh(self.mesh, color="lightgrey", show_edges=True, opacity=1.0, name="volume")
        self.plotter.add_mesh(self.nodes_polydata, color="darkgrey", point_size=12, render_points_as_spheres=True, name="nodes_points")
        self.plotter.reset_camera()
        self.plotter.render()
        self._setup_picking()

    # Meshio to PyVista
    def _meshio_to_pyvista(self, meshio_mesh):
        points = meshio_mesh.points
        all_cells = []
        all_cell_types = []
        elements = []

        for cell_block in meshio_mesh.cells:
            ctype = cell_block.type
            if ctype not in VTK_CELL_TYPES:
                continue
            cells = np.asarray(cell_block.data, dtype=int)
            if cells.size == 0:
                continue
            n_points = cells.shape[1]
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
            self.plotter.disable_picking()
        except Exception:
            pass

    # ------------------------ RANDBEDINGUNGEN ------------------------
    def add_boundary_condition(self):
        items = ["Kraftrandbedingung", "Verschiebungsrandbedingung", "Lager"]
        item, ok = QInputDialog.getItem(self, "Randbedingung auswählen", "Typ der Randbedingung:", items, 0, False)
        if not (ok and item):
            return

        bc = {"type": item, "nodes": set(), "faces": set(),
              "values": [0, 0, 0] if item != "Lager" else None,
              "tree_items": {}}
        self.boundary_conditions.append(bc)
        index = len(self.boundary_conditions) - 1
        self.current_bc_index = index

        bc_item = QTreeWidgetItem(self.bc_root_item, [item])
        bc_item.setFont(0, QFont("", weight=QFont.Weight.Bold))
        bc_item.setBackground(0, QBrush(QColor(220, 220, 220)))
        bc_item.setExpanded(True)

        # Dropdown Auswahlmodus
        sel_item = QTreeWidgetItem(bc_item, ["Auswahlmodus"])
        combo = QComboBox()
        combo.addItems(["Knoten", "Fläche"])
        combo.currentIndexChanged.connect(lambda idx, bc_idx=index: self.set_bc_selection_mode(bc_idx, idx))
        self.bc_header_tree.setItemWidget(sel_item, 0, combo)

        # Dropdown Auswahlmodus2
        sel_item2 = QTreeWidgetItem(bc_item, ["Auswahlmodus"])
        combo2 = QComboBox()
        if item == "Kraftrandbedingung":
            combo2.addItems(["Kraft", "Drehmoment"])
        else:
            combo2.addItems(["Translation", "Rotation"])
        self.bc_header_tree.setItemWidget(sel_item2, 0, combo2)

        # Zuordnung
        list_item = QTreeWidgetItem(bc_item, ["Zuordnung"])
        list_widget = QListWidget()
        list_widget.setFixedHeight(100)
        self.bc_header_tree.setItemWidget(list_item, 0, list_widget)

        # Werte / Locks
        val_item = QTreeWidgetItem(bc_item, [""])
        val_widget = QWidget()
        layout = QVBoxLayout(val_widget)
        layout.setSpacing(2)
        layout.setContentsMargins(2, 2, 2, 2)

        if item == "Lager":
            locks = {}
            for dof in ["Tx","Ty","Tz","Rx","Ry","Rz"]:
                cb = QCheckBox(dof)
                layout.addWidget(cb)
                locks[dof] = cb
            bc["tree_items"]["locks"] = locks
            x_input = y_input = z_input = None
        else:
            x_input = QLineEdit("0"); layout.addWidget(QLabel("x:")); layout.addWidget(x_input)
            y_input = QLineEdit("0"); layout.addWidget(QLabel("y:")); layout.addWidget(y_input)
            z_input = QLineEdit("0"); layout.addWidget(QLabel("z:")); layout.addWidget(z_input)

        self.bc_header_tree.setItemWidget(val_item, 0, val_widget)

        # store references
        bc["tree_items"]["sel"] = combo
        bc["tree_items"]["list"] = list_widget
        bc["tree_items"]["x"] = x_input
        bc["tree_items"]["y"] = y_input
        bc["tree_items"]["z"] = z_input
        bc["tree_items"]["tree_item_ref"] = bc_item

    def set_bc_selection_mode(self, bc_index, idx):
        self.current_bc_index = bc_index
        self.selection_mode = "Knoten" if idx == 0 else "Fläche"

    # ------------------------ MATERIAL ------------------------
    def define_material(self, mat_dir_model, mat_model, layout):
        # Alte Widgets löschen
        self.clear_layout(layout)

        E_Modul = []
        v = []

        if mat_dir_model == "isotrop" and mat_model == "linear elastisch":
            E_Modul.append(QLineEdit("0"))
            v.append(QLineEdit("0"))
            layout.addWidget(QLabel("E Modul:"))
            layout.addWidget(E_Modul[0])
            layout.addWidget(QLabel("Querkontraktionszahl:"))
            layout.addWidget(v[0])
        elif mat_dir_model == "orthotrop" and mat_model == "linear elastisch":
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
        self.mat_header_tree.setItemWidget(self.mat_params_tree_item, 0, self.mat_params_widget)  # <- wichtig
        self.mat_header_tree.viewport().update()
        self.mat_root_item.setExpanded(True)


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
        actor = self.plotter.add_mesh(sphere, color="red", opacity=1.0)
        self.node_selection_actors[pid] = actor
        if self.current_bc_index is not None:
            bc = self.boundary_conditions[self.current_bc_index]
            bc["nodes"].add(pid)
            list_widget = bc["tree_items"]["list"]
            list_widget.addItem(f"Knoten {pid}: {coord[0]:.3f},{coord[1]:.3f},{coord[2]:.3f}")
        self.plotter.render()

    def _deselect_node(self, pid):
        if pid not in self.selected_nodes:
            return
        self.selected_nodes.remove(pid)
        actor = self.node_selection_actors.pop(pid, None)
        if actor:
            try: self.plotter.remove_actor(actor)
            except: pass
        if self.current_bc_index is not None:
            bc = self.boundary_conditions[self.current_bc_index]
            bc["nodes"].discard(pid)
            list_widget = bc["tree_items"]["list"]
            for i in range(list_widget.count()):
                if list_widget.item(i).text().startswith(f"Knoten {pid}:"):
                    list_widget.takeItem(i)
                    break
        self.plotter.render()

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
            actor = self.plotter.add_mesh(face_poly, color=color, opacity=1.0, reset_camera=False, show_edges=False)
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
            self.plotter.render()
        except Exception as e:
            print("Fehler beim Hervorheben der Fläche:", e)

    def _deselect_face(self, surf_cell_idx):
        if surf_cell_idx not in self.selected_faces:
            return
        self.selected_faces.remove(surf_cell_idx)
        actor = self.face_selection_actors.pop(surf_cell_idx,None)
        if actor:
            try: self.plotter.remove_actor(actor)
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
        self.plotter.render()

    # ------------------------ UTILITY ------------------------
    def _clear_plotter_keep_plotter(self):
        try:
            self.plotter.clear()
        except Exception:
            pass

# ------------------------ MAIN ------------------------
def main():
    app = QApplication(sys.argv)
    window = FEMGUI()
    window.show()
    sys.exit(app.exec())

if __name__ == "__main__":
    main()
