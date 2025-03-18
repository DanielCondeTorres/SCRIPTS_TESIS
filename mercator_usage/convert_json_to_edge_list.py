import json
import networkx as nx
import matplotlib.pyplot as plt

# Cargar el JSON usando 'utf-8-sig' para manejar el BOM
with open("file.json", "r", encoding="utf-8-sig") as f:
    data = json.load(f)

# Extraer la sección de red
network = data.get("network", {})
nodes = network.get("items", [])
edges = network.get("edges", network.get("links", []))

# Crear el grafo (usa nx.DiGraph() si los enlaces son dirigidos)
G = nx.Graph()

# Agregar los nodos al grafo, usando el 'id' y guardando la etiqueta
for node in nodes:
    G.add_node(node["id"], label=node.get("label", str(node["id"])))

# Agregar los enlaces utilizando las claves "source_id", "target_id" y "strength" como peso
for edge in edges:
    if "source_id" in edge and "target_id" in edge:
        source_id = edge["source_id"]
        target_id = edge["target_id"]
        weight = edge.get("strength", 1)  # Se asigna 1 si no se encuentra 'strength'
        G.add_edge(source_id, target_id, weight=weight)
    else:
        print("No se encontró 'source_id' ni 'target_id' en el enlace:", edge)

# Guardar la edge list en un archivo de texto llamado "edgelist.txt"
nx.write_edgelist(G, "edgelist.txt", data=True, delimiter=" ")

# Dibujar el grafo utilizando una disposición de fuerza de resorte
pos = nx.spring_layout(G)
labels = nx.get_node_attributes(G, "label")

plt.figure(figsize=(12, 12))
nx.draw_networkx_nodes(G, pos, node_size=500, node_color="skyblue")
nx.draw_networkx_edges(G, pos, width=1.0, alpha=0.7)
nx.draw_networkx_labels(G, pos, labels, font_size=10)
plt.title("Grafo generado a partir del JSON de VOSviewer")
plt.axis("off")
plt.show()

