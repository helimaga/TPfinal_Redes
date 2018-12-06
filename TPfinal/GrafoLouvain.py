#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Dec  5 12:29:38 2018

@author: santiago

"""

import community as cm
import networkx as nx
import matplotlib.pyplot as plt

# Lo que está comentado abajo de esta línea es un grafo chico de ejemplo para ver cómo funciona el programa.
'''
plt.close('all')
# Replace this with your networkx graph loading depending on your format !
G = [("a","b"),("b","c"), ("b","d"), ("d","e"), ("e","f"), ("b","f"), ("b","g"), ("f", "g"), ("a","g"), ("a","g"), ("c","d"), ("d", "g"), ("d","h"), ("aa","h"), ("aa","c"), ("f", "h"), ("x", "z"), ("a", "z"), ("x", "y"), ("x", "w"), ("w", "y"), ("z", "w"), ("z", "y")]

G=nx.Graph(G)
#first compute the best partition
partition = cm.best_partition(G)

#drawing
size = float(len(set(partition.values())))
pos = nx.spring_layout(G)
count = 0.
for com in set(partition.values()) :
    count = count + 1.
    list_nodes = [nodes for nodes in partition.keys()
                                if partition[nodes] == com]
    pos=nx.kamada_kawai_layout(G)
    nx.draw_networkx_nodes(G, pos, list_nodes, node_size = 50,
                                node_color = str(count / size))


nx.draw_networkx_edges(G, pos, alpha=0.5,with_labels=True)
plt.suptitle('Grafo original',fontsize=20)
plt.show()
'''

def GrafoLouvain(G):
# Se crea el grafo LouvainFinalGraph, que es el último grafo alcanzado por Louvain (al que no se le puede optimizar la modularidad):

    LouvainFinalGraph=G # Primero, se tiene el grafo original.
    pos=nx.kamada_kawai_layout(G)
    
    particion= cm.best_partition(G)  # diccionario con las comunidades según Louvain.

    # Lo de los colores no supe hacerlo:
    colores = ['aqua',
         'blueviolet',
         'peru',
         'darkgrey',
         'darkseagreen',
         'violet',
         'red',
         'hotpink',
         'moccasin',
         'green',
         'lime',
         'mediumpurple',
         'sandybrown',
         'orangered',
         'peru',
         'royalblue',
         'skyblue',
         'teal',
         'yellow']


    plt.figure()
    nx.draw(LouvainFinalGraph,
        pos,
        width=0.1,
        edge_color = 'k',
        node_color = list(particion.values()),
        node_size=50,
        font_size=20,
        with_labels=True,
        )
    plt.suptitle('Grafo incial de Louvain',fontsize=20)
    plt.show()


    for nodo_i in LouvainFinalGraph:
        for nodo_j in LouvainFinalGraph:
            if particion[nodo_i]==particion[nodo_j]:    # si forman parte de la misma comunidad
                LouvainFinalGraph=nx.contracted_nodes(LouvainFinalGraph,nodo_i,nodo_j)  # se contraen.
    
    
    plt.figure()
    nx.draw(LouvainFinalGraph,
            pos,
            width=0.1,
            edge_color = 'k',
            #node_color = ,# No puedo arreglar esto. No entiendo cómo se hace.
            node_size=50,
            font_size=20,
            with_labels=True,
            )
    plt.suptitle('Grafo final de Louvain',fontsize=20)
    plt.show()


