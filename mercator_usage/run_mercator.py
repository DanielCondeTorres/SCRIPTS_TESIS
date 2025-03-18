import mercator

# Especifica la ruta de tu archivo edgelist (ej: "edgelist.txt")
edgelist_filename = "edgelist_networkx.txt"

# Genera la incrustación hiperbólica
embedding = mercator.embed(edgelist_filename)

# Guarda los resultados (si es necesario)
with open("embedding_result.txt", "w") as f:
    for node, coords in embedding.items():
        f.write(f"{node} {coords[0]} {coords[1]}\n")
