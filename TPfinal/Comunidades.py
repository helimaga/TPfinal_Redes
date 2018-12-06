#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Nov 11 13:44:26 2018

@author: santiago

Descripción:    Archivo de comunidades, basado en las funciones de Juan para el TC3.
"""

#%%
import networkx as nx
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import community as cm
import sys
import igraph
import itertools as itr
from operator import itemgetter
#%%

#paths

pathHeli = '/home/heli/Documents/Redes/TPfinal/'
pathJuancho = '/home/gossn/Dropbox/Documents/Materias_doctorado/RedesComplejasBiologicos/redes2018/TPfinal/'
pathSanti = '/home/santiago/Documentos/RC/Redes2018/TPfinal/'
pathDocente = '?'

#%%
path = pathHeli

sys.path.append(path)
import modularity_max

#%% Particionamos segun cada uno de los metodos (encapsulamos en funcion...)
def Communities(red,path,metodo):
    if metodo=='l':
        # Metodología Louvain networkx
        red_part= cm.best_partition(red)

    elif metodo=='fg':
    # Metodología Fast Greedy (max Modularity) networkx
    # de: https://networkx.github.io/documentation/latest/reference/algorithms/generated/networkx.algorithms.community.modularity_max.greedy_modularity_communities.html
    
    # para importar modulos provenientes de la version trial de networkx...
        red_part_FGreedy0 = list(modularity_max.greedy_modularity_communities(red))
        red_part_FGreedy0 = [list(x) for x in red_part_FGreedy0]
        
        red_part = {}
        for n in list(red.nodes()):
            for comm in range(len(red_part_FGreedy0)):
                if n in red_part_FGreedy0[comm]:
                    red_part[n] = comm
    elif metodo=='im':
        
        # Metodología Infomap
        #https://www.youtube.com/watch?v=mO0J_H4YLJA
        nx.write_gml(red, path + 'redCommunities.gml', stringizer=None)
        redI = igraph.Graph.Read_GML(path + 'redCommunities.gml') # VER ESTE PATH.
        red_part_IMap0 = list(redI.community_infomap())

        red_part = {}
        vs = igraph.VertexSeq(redI)
        for n in range(len(red.nodes())):
            vStr = vs[n]['label']
            for comm in range(len(red_part_IMap0)):
                if n in red_part_IMap0[comm]:
                    red_part[vStr] = comm
    elif metodo=='ng':
        # Metodología Newman-Girvan (Edge Betweenness) networkx
        # http://materias.df.uba.ar/redesa2018c2/files/2018/10/15_Clusters_2.pdf... 
        # = Ejercicio TC1.2c  ... pero continuar la division un paso mas? (4 grupos?
        # Hay método implementado?
        
        red_part_NewGir0 = list(nx.algorithms.community.centrality.girvan_newman(red))
        
        red_part = [] #lista de diccionarios, cada diccionario es una particion
        for numPart in range(len(red_part_NewGir0)): # len(red_part_NewGir0) = N-1
            red_part_NewGir = {}
            for n in list(red.nodes()):
                for comm in range(len(red_part_NewGir0[numPart])):
                    if n in red_part_NewGir0[numPart][comm]:
                        red_part_NewGir[n] = comm
            red_part.append(red_part_NewGir)
    return red_part
#%% Definimos Silhouette
def silhouetteJuancho(graph,commPartition,outputOpt):#modificada respecto de TC3
    numComm = max(commPartition.values())+1
    silhouette = {}
    silhouetteAvg = []
    nodos = sorted(list(graph.nodes()))
    for nSource in nodos:
        cnk = []
        for comm in range(numComm):
            nNodesComm = 0
            dST = 0
            for nTarget in graph.nodes():
                if (commPartition[nTarget] == comm) & (nTarget != nSource):
                    nNodesComm += 1
                    dST += nx.shortest_path_length(graph,source=nSource, target=nTarget)
            nNodesComm = max(nNodesComm,1) # por si nSource es solitario
            dST = max(dST,1) # por si nSource es solitario
            cnk.append(dST/nNodesComm)
        commSource = commPartition[nSource]
        an = cnk[commSource]
        del cnk[commSource]
        bn = min(cnk)
        s0 = (bn-an)/max(an,bn)
        silhouette[nSource]=s0
    silhouetteAvg = np.mean(list(silhouette.values()))
    if outputOpt == 'all':
        return silhouette
    elif outputOpt == 'mean':
        return silhouetteAvg
#%% Maximize Diagonal... intercambia filas de una matriz para maximizar "Traza(Mat) - NO-Traza(Mat)"
#def MaximizeDiagonal(matrix):
#    import numpy as np
#    import itertools as itr
#    from operator import itemgetter
#    
##    matrix = [[10,56,80],[1321,2,15],[3,540,123],[0,2,3]]
#    
#    N,M = np.shape(matrix)
#    
#    posPerm = list(itr.permutations(range(N)))
#    
#    P = np.shape(posPerm)[0]
#    
#    
#    costF = []
#    for p in range(P):
#        B = []
#        perm = posPerm[p]
#        for n in perm:
#            C = []
#            for m in range(M):
#                C.append(matrix[n][m])
#            B.append(C)
#        costF0=0
#        for n in range(N):
#            for m in range(M):
#                if n==m:
#                    costF0 = costF0+B[n][m]
#                else:
#                    costF0 = costF0-B[n][m]
#        costF.append(costF0)
#    
#    
#    bestPermIdx = max(enumerate(costF), key=itemgetter(1))[0]
#    bestPerm = posPerm[bestPermIdx]
#    
#    MatPerm = []
#    for n in range(len(bestPerm)):
#        MatPerm.append(matrix[bestPerm[n]])
#    return MatPerm,bestPerm
#%% 
def dictsValues2Mat(AD,BD):
    # Readapta diccionario A (o B, el de menor numero de comunidades) para maximizar coincidencias...
    # https://stackoverflow.com/questions/835092/python-dictionary-are-keys-and-values-always-the-same-order
    # If items(), keys(), values(),  iteritems(), iterkeys(), and  itervalues() are called with no intervening modifications to the dictionary, the lists will directly correspond.
    
    keyA = list(AD.keys())
    keyB = list(BD.keys())
    if keyA != keyB:
        raise NotImplementedError('keys son distintassss!')
    
    A = list(AD.values())
    B = list(BD.values())
    
    N = len(A)
    
    permA = list(set(A))
    permB = list(set(B))
    
    permMin = min(permA,permB, key=len)
    
    if permMin == permA: # D es el de menor numero de especies...
        C = B
        D = A
    else:
        C = A
        D = B
    
    
    posPerm = list(itr.permutations(permMin))
    
    costF = []
    for p in range(len(posPerm)):
        DPerms = []
        for d in range(len(D)):
            DPerms.append(posPerm[p][D[d]])
        costF0 = 0
        for n in range(N):
            if C[n]==DPerms[n]:
                costF0+=1
        costF.append(costF0)
    
    bestPermIdx = max(enumerate(costF), key=itemgetter(1))[0]
    DPerm = []
    for d in range(len(D)):
        DPerm.append(posPerm[bestPermIdx][D[d]])
    
    if permMin == permA:
        AD = dict(zip(keyA, DPerm))
        BD = dict(zip(keyB, C))
    else:
        AD = dict(zip(keyA, C))
        BD = dict(zip(keyB, DPerm))
    
    return AD,BD