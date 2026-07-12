import gmsh


# Initialisierung
gmsh.initialize()

# Mesh-Datei öffnen (ersetze durch deinen Dateinamen)
mesh_file = "C:/Users/andre/OneDrive/Desktop/myFEM/wuerfel.msh"
gmsh.open(mesh_file)

# Alle Flächen (dim=2) abrufen
surfaces = gmsh.model.getEntities(dim=2)  # dim=2 → surfaces

# Physical Groups für jede Fläche erstellen
for i, (dim, tag) in enumerate(surfaces):
    # Sicherstellen, dass tag ein Integer ist
    tag_int = int(tag)
    # Physical Group erzeugen
    # Name optional, nur Gmsh >= 4.7 unterstützt den Namen
    try:
        gmsh.model.addPhysicalGroup(dim, [tag_int], f"Surface_{tag_int}")
    except TypeError:
        # Fallback: Name weglassen, falls TypeError wegen Name auftritt
        gmsh.model.addPhysicalGroup(dim, [tag_int])

# Optional: alle Physical Groups ausgeben
print("Alle Physical Groups:")
for dim, tag in gmsh.model.getPhysicalGroups():
    try:
        name = gmsh.model.getPhysicalName(dim, tag)
    except:
        name = "<kein Name>"
    print(f"dim={dim}, tag={tag}, name={name}")

# Mesh mit Physical Groups speichern
gmsh.write("C:/Users/andre/OneDrive/Desktop/myFEM/dein_mesh_with_physical.msh")

# Gmsh beenden
gmsh.finalize()
