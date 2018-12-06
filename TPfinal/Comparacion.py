#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Distancia entre gráficos.

Ver: http://www.xavierdupre.fr/app/ensae_teaching_cs/helpsphinx3/specials/graph_distance.html#l-graph-distance
"""
import networkx as nx
import sys
import matplotlib.pyplot as plt

pathHeli = '/home/heli/Documents/Redes/Practicas/TPs/Redes2018/TPfinal/'
pathJuancho = '/home/gossn/Dropbox/Documents/Materias_doctorado/RedesComplejasBiologicos/redes2018/TPfinal/'
pathSanti = '/home/santiago/Documentos/RC/Redes2018/TPfinal/'
pathDocente = '?'

#%%
path = pathSanti

sys.path.append(path)
import graph_distance as GD
Clase = GD.GraphDistance() # importo la clase en la que está la función.



#%%
# Dos grafos chicos como ejemplo:

graph1 = [("a","b"), ("b","c"), ("b","d"), ("d","e"), ("e","f"), ("b","f"), ("b","g"), ("f", "g"), ("a","g"), ("a","g"), ("c","d"), ("d", "g"), ("d","h"), ("aa","h"), ("aa","c"), ("f", "h")]
graph2 = GD.copy.deepcopy(graph1) + [ ("h", "m"), ("m", "l"), ("l", "C"), ("C", "r"), ("a", "k"), ("k", "l"), ("k", "C")]

graph1 = nx.Graph(graph1)
graph2 = nx.Graph(graph2)

Grafos=[graph1,graph2]
ColorGrafos=['red','blue']

# Visualizo los dos grafos:
for i in range(len(Grafos)):
    plt.figure()
    nx.draw(Grafos[i],
            width=0.1,
            edge_color = 'k',
            node_color = ColorGrafos[i], 
            node_size=50,
            font_size=10,
            with_labels=False,
            )
    plt.show()

# Esta línea no anda. "distance_matching_graphs_paths" es una función de la clase "GraphDistance" del archivo "graph_distance".

distance, graph = graph1.distance_matching_graphs_paths(graph2, use_min=False)#, store=store)
