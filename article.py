import pandas as pd
import numpy as np
import networkx as nx
import nearmincut as nmc
import matplotlib.pyplot as plt
from pymatgen.core import Structure
from matplotlib.animation import FuncAnimation, PillowWriter, FFMpegWriter


def euclidean_distance(p1, p2):
    dist = np.sqrt(np.sum((p1-p2)**2, axis=0))
    return dist

def get_cut_edges(G, partition):
    """
    Input: G es el grafo del material. Partition es la 2-particion (X, Y con X intersec Y = vacio) de nodos que genera el algoritmo min_cut
    Output: lista de arcos que se deben eliminar con la 2-particion dada
    """
    H1 = G.subgraph(partition[0])  # subgrafo inducido por el subconjunto X de la particion
    # subgrafo inducido por el subconjunto Y de la particion
    H2 = G.subgraph(partition[1])
    G2 = nx.Graph()
    # se crea un grafo a partir de los subgrafos inducidos
    G2.add_nodes_from(G.nodes)
    G2.add_edges_from(H1.edges)
    G2.add_edges_from(H2.edges)
    # los arcos del grafos G-G2 son los arcos que se deben eliminar con la 2-particion dada
    R = nx.difference(G, G2)
    cut_edges = list(R.edges())
    return cut_edges

def plot_2d_graph(G, df, color_map):
    pos = nx.circular_layout(G)
    nx.draw(G, pos, labels={node: df['Species'].to_list()[
            node] for node in G.nodes()}, node_color=color_map)
    labels = {(x[0], x[1]): round(x[2], 4)
              for x in list(G.edges.data("capacity"))}
    nx.draw_networkx_edge_labels(
        G, pos, edge_labels=labels, font_color='red', font_size=7)
    plt.show()
    return

def plot_3d_graph(G, cartesian_coords, color_set, plot_cut, cut_edges=None):
    coord_x, coord_y, coord_z = map(list, zip(*cartesian_coords))
    label_atoms = list(nx.get_node_attributes(G, "Species").values())
    fig = plt.figure()
    ax = fig.add_subplot(projection='3d')
    # Plot nodes
    for x, y, z, l in zip(coord_x, coord_y, coord_z, label_atoms):
        ax.scatter(x, y, z, label=l, color=color_set[l], s=100)
    # Plot edges
    for e in G.edges():
        u, v = e[0], e[1]
        x_i = [coord_x[u], coord_x[v]]
        y_i = [coord_y[u], coord_y[v]]
        z_i = [coord_z[u], coord_z[v]]
        color = 'black'
        if plot_cut:  # plot edges in nearmincut algorithm
            if (u, v) in cut_edges or (v, u) in cut_edges:
                color = 'red'
        ax.plot3D(x_i, y_i, z_i, color=color, linewidth=0.5)
    # Set the axis labels
    handles, labels = plt.gca().get_legend_handles_labels()
    by_label = dict(zip(labels, handles))
    ax.set_xlabel('x-axis')
    ax.set_ylabel('y-axis')
    ax.set_zlabel('z-axis')
    ax.legend(by_label.values(), by_label.keys(), loc='best')
    plt.show()
    return

# Read the CIF File
filename = "CrN.cif"
structure = Structure.from_file(filename)  # unit cell
super_cell = structure.copy()
super_cell.make_supercell([2, 2, 2])

df = super_cell.as_dataframe()
label_set = dict(zip(df['Species'].unique(), ['Cr', 'N']))
color_set = dict(zip(['Cr', 'N'], ['purple', 'yellow']))

# Create graph given a spherical neighborhood
G = nx.Graph()
for i, site in enumerate(super_cell):
    G.add_node(i)
    G.nodes[i]["Species"] = label_set[site.species]
    G.nodes[i]["position"] = (site.x, site.y, site.z)

for i, site in enumerate(super_cell):
    neighbors = [(n.index, n.nn_distance)
                 for n in super_cell.get_neighbors(site, r=3)]
    for n in neighbors:
        G.add_edge(i, n[0], weight=n[1])

# Plot graph (spherical neighborhood)
d = nx.get_node_attributes(G, "position")
plot_3d_graph(G, d.values(), color_set, False)

# Create graph using euclidean distance
lattice = super_cell.lattice
fractional_coords = super_cell.frac_coords
# Convert the fractional coordinates to Cartesian coordinates using the lattice vectors
cartesian_coords = lattice.get_cartesian_coords(fractional_coords)
distances = []
N = len(cartesian_coords)
for i in range(N):
    p1 = cartesian_coords[i]
    dist_i = {}
    for j in range(N):
        p2 = cartesian_coords[j]
        if j != i:
            dist_i[j] = euclidean_distance(p1, p2)
    distances.append(dist_i)

G2 = nx.Graph()
for i, site in enumerate(super_cell):
    G2.add_node(i)
    G2.nodes[i]["Species"] = label_set[site.species]
for i in range(N):
    for key, value in distances[i].items():
        if value <= 2.5:  # metric for connection of 2 atoms
            eV = 9.462 #energy between nearest atoms
            G2.add_edge(i, key, capacity=value)

plot_3d_graph(G2, cartesian_coords, color_set, False)

# Max-flow min-cut networkx
s, t = 0, 1
cut_value, partition = nx.minimum_cut(G2, s, t)
print(cut_value)
cut_edges = get_cut_edges(G2, partition)

plot_3d_graph(G2, cartesian_coords, color_set, True, cut_edges)
"""
#near Max-flow min-cut
kwargs={'eps':0.125}
s, t = 0, 2
sol = nmc.near_min_cuts(G,s,t,**kwargs)
cut_edges = sol[0][1]
"""