import pandas as pd
import numpy as np
import networkx as nx
import nearmincut as nmc
import matplotlib.pyplot as plt
from pymatgen.core import Structure

def get_space_distro_matrix(filename):
    """
    Input: Nombre del archivo .cif de la estructura
    Output: dataframe de la distribucion en el espacio de los atomos
    """
    material = Structure.from_file(filename).as_dict()
    df = pd.DataFrame()
    atoms, a, b, c = [], [], [], []
    for d in material['sites']:
        atoms.append(d['label'])
        space_distro = d['abc'] #d['xyz']
        a.append(space_distro[0])
        b.append(space_distro[1])
        c.append(space_distro[2])
    df['Atom'], df['a'], df['b'], df['c']  = atoms, a, b, c
    return df

def space_euclidean_distance(p1, p2):
    """
    Input: 2 puntos
    Output: Distancia euclidania entre esos 2 puntos
    """
    dist = np.sqrt(np.sum((p1-p2)**2, axis=0))
    return dist

def get_cut_edges(G, partition):
    """
    Input: G es el grafo del material. Partition es la 2-particion (X, Y con X intersec Y = vacio) de nodos que genera el algoritmo min_cut
    Output: lista de arcos que se deben eliminar con la 2-particion dada
    """
    H1 = G.subgraph(partition[0]) #subgrafo inducido por el subconjunto X de la particion
    H2 = G.subgraph(partition[1]) #subgrafo inducido por el subconjunto Y de la particion
    G2 = nx.Graph()
    G2.add_nodes_from(G.nodes) #se crea un grafo a partir de los subgrafos inducidos
    G2.add_edges_from(H1.edges)
    G2.add_edges_from(H2.edges)
    R = nx.difference(G, G2) #los arcos del grafos G-G2 son los arcos que se deben eliminar con la 2-particion dada
    cut_edges = list(R.edges())
    return cut_edges

def plot_2d_graph(G, df, color_map):
    pos = nx.circular_layout(G)
    nx.draw(G, pos, labels={node: df['Atom'].to_list()[node] for node in G.nodes()}, node_color=color_map)
    labels={(x[0], x[1]) : round(x[2], 4) for x in list(G.edges.data("capacity"))}
    nx.draw_networkx_edge_labels(G, pos, edge_labels=labels, font_color='red', font_size=7)
    plt.show()
    return

def plot_3d_graph(G, df, color_set, plot_cut, cut_edges=None):
    fig = plt.figure(constrained_layout=True)
    ax = fig.add_subplot(projection="3d")

    for index, row in df.iterrows():
        x, y, z = row['a'], row['b'], row['c']
        atom = row['Atom']
        ax.scatter(x, y, z, color=color_set[atom], label=atom)

    for e in G.edges():
        u, v = e[0], e[1]
        x_i = [df['a'].to_list()[u], df['a'].to_list()[v]]
        y_i = [df['b'].to_list()[u], df['b'].to_list()[v]]
        z_i = [df['c'].to_list()[u], df['c'].to_list()[v]]
        color = 'black'
        if plot_cut:
            if (u,v) in cut_edges or (v,u) in cut_edges: color = 'red'
        ax.plot3D(x_i, y_i, z_i, color=color)

    handles, labels = plt.gca().get_legend_handles_labels()
    by_label = dict(zip(labels, handles))
    ax.set_xlabel('x-axis')
    ax.set_ylabel('y-axis')
    ax.set_zlabel('z-axis')
    ax.legend(by_label.values(), by_label.keys(), loc='best')
    plt.show()
    return

#df = pd.read_csv('data.csv', sep=';')
filename = 'MgTiO3_conventional_standard.cif'
df = get_space_distro_matrix(filename)

#Distances
color_set = {'Mg': 'green', 'Ti':'blue', 'O':'yellow'}
color_map = []
distances = []
for i, row in df.iterrows():
    color_map.append(color_set[row['Atom']])
    p1 = np.array([row['a'], row['b'], row['c']])
    dist_i = {}
    for j, row2 in df.iterrows():
        p2 = np.array([row2['a'], row2['b'], row2['c']])
        if j != i:
            dist_i[j] = space_euclidean_distance(p1, p2)
    distances.append(dist_i)

#Create Graph
N = df.shape[0]
G = nx.Graph()
G.add_nodes_from(np.arange(N))
for i in range(N):
    min_dist = min(distances[i].values())
    max_dist = max(distances[i].values())
    d = np.mean(list(distances[i].values()))
    for key, value in distances[i].items():
        if value <= d: #metrica para la conexion de 2 atomos
            G.add_edge(i, key, capacity=value)

#Plot graph in 2D
#plot_2d_graph(G, df, color_map)

#Max-flow min-cut networkx
s, t = 0, 2
cut_value, partition = nx.minimum_cut(G, s, t)
print(cut_value)
cut_edges = get_cut_edges(G, partition)

"""
#near Max-flow min-cut
kwargs={'eps':0.125}
s, t = 0, 2
sol = nmc.near_min_cuts(G,s,t,**kwargs)
cut_edges = sol[0][1]
"""

#Plot graph in 3D
plot_3d_graph(G, df, color_set, False)
plot_3d_graph(G, df, color_set, True, cut_edges)