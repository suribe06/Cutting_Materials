import pandas as pd
import numpy as np
import networkx as nx
import nearmincut as nmc
import matplotlib.pyplot as plt
from pymatgen.core import Structure
from interpo import energyFunctionCrN, energyFunctionCrCr, energyFunctionNN

def euclidean_distance(p1, p2):
    dist = np.sqrt(np.sum((p1-p2)**2, axis=0))
    return dist

def getDistanceMatrix(cartesian_coords):
    N = len(cartesian_coords)
    distances = np.zeros((N,N))
    for i in range(N):
        for j in range(i, N):
            p1 = cartesian_coords[i]
            p2 = cartesian_coords[j]
            distances[i][j] = distances[j][i] = euclidean_distance(p1, p2)
    return distances

def createSphericalNeighborhoodGraph(super_cell):
    # Create graph given a spherical neighborhood
    G = nx.Graph()
    for i, site in enumerate(super_cell):
        G.add_node(i)
        G.nodes[i]["Species"] = label_set[site.species]
    radius = 3
    for i, site in enumerate(super_cell):
        neighbors = [(n.index, n.nn_distance) for n in super_cell.get_neighbors(site, r=radius)]
        for n in neighbors:
            G.add_edge(i, n[0], weight=n[1])
    return G

def createDistanceGraph(cartesian_coords, label_set):
    cutoff_distance = 2.5  # Umbral de distancia para conexión
    N = len(cartesian_coords)
    distances = getDistanceMatrix(cartesian_coords)
    G = nx.DiGraph()
    for i, site in enumerate(super_cell):
        G.add_node(i)
        G.nodes[i]["Species"] = label_set[site.species]
    for i in range(N):
        for j in range(i+1, N):
            dist = distances[i][j]
            if dist <= cutoff_distance:  # metric for connection of 2 atoms
                eV = 9.462 #energy between nearest atoms
                G.add_edge(i, j, capacity=eV, distance=dist)
                G.add_edge(j, i, capacity=eV, distance=dist)
    return G

def createEnergyGraph(cartesian_coords, df, label_set):
    N = len(cartesian_coords)
    energy_matrix = np.zeros((N, N), dtype=float)
    distances = getDistanceMatrix(cartesian_coords)
    f = energyFunctionCrCr()
    g = energyFunctionNN()
    h = energyFunctionCrN()
    for i in range(N):
        for j in range(i+1, N):
            species_i = label_set[df.iloc[i]['Species']]
            species_j = label_set[df.iloc[j]['Species']]
            distance = distances[i,j]
            if species_i == 'Cr' and species_j == 'Cr': energy = f(distance)
            elif species_i == 'N' and species_j == 'N': energy = g(distance)
            else: energy = h(distance)
            energy_matrix[i,j] = energy_matrix[j,i] = energy
    G = nx.DiGraph()
    for i, site in enumerate(super_cell):
        G.add_node(i)
        G.nodes[i]["Species"] = label_set[site.species]
    for i in range(N):
        for j in range(i+1, N):
            G.add_edge(i, j, capacity=energy_matrix[i,j])
            G.add_edge(j, i, capacity=energy_matrix[i,j])
    return G

def getCutEdges(G, partition):
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

def plotGraph3D(G, cartesian_coords, color_set, plot_cut, cut_edges=None):
    coord_x, coord_y, coord_z = map(list, zip(*cartesian_coords))
    label_atoms = list(nx.get_node_attributes(G, "Species").values())
    fig = plt.figure()
    ax = fig.add_subplot(projection='3d')
    # Plot nodes
    for x, y, z, l in zip(coord_x, coord_y, coord_z, label_atoms):
        ax.scatter(x, y, z, label=l, color=color_set[l], s=100)
    # Plot edges
    plotted_edges = set()
    edges = G.edges()
    for u, v in edges:
        # Check if the edge has been plotted in the opposite direction
        if (v, u) in plotted_edges:
            continue
        x_i = [coord_x[u], coord_x[v]]
        y_i = [coord_y[u], coord_y[v]]
        z_i = [coord_z[u], coord_z[v]]
        color = 'black'
        if plot_cut:  # plot edges in nearmincut algorithm
            if (u, v) in cut_edges or (v, u) in cut_edges:
                color = 'red'
        ax.plot3D(x_i, y_i, z_i, color=color, linewidth=0.5)
        plotted_edges.add((u, v))
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

# Plot graph (spherical neighborhood)
#G1 = createSphericalNeighborhoodGraph(super_cell)
#d = nx.get_node_attributes(G1, "position")
#plotGraph3D(G1, d.values(), color_set, False)

# Create graph using euclidean distance
lattice = super_cell.lattice
fractional_coords = super_cell.frac_coords
# Convert the fractional coordinates to Cartesian coordinates using the lattice vectors
cartesian_coords = lattice.get_cartesian_coords(fractional_coords)
G2 = createDistanceGraph(cartesian_coords, label_set)
plotGraph3D(G2, cartesian_coords, color_set, False)

# Create graph using energy
#G3 = createEnergyGraph(cartesian_coords, df, label_set)
#plotGraph3D(G3, cartesian_coords, color_set, False)

"""
# Max-flow min-cut networkx
s, t = 0, 1
cut_value, partition = nx.minimum_cut(G2, s, t)
print(cut_value)
cut_edges = getCutEdges(G2, partition)

plotGraph3D(G2, cartesian_coords, color_set, True, cut_edges)

#near Max-flow min-cut
kwargs={'eps':0.125}
s, t = 0, 2
sol = nmc.near_min_cuts(G,s,t,**kwargs)
cut_edges = sol[0][1]
"""