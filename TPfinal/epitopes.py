#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Nov  2 14:58:22 2018

@author: heli
"""
#%%

#librerias

import pandas as pd
import itertools as itr
import networkx as nx
import matplotlib.pyplot as plt
import numpy as np
from Bio import pairwise2
from Bio.SubsMat.MatrixInfo import blosum62
import matplotlib.patches as mpatches
import csv
import sys



#%%

#paths

pathHeli = '/home/heli/Documents/Redes/Practicas/TPs/Redes2018/TPfinal/'
pathJuancho = '/home/gossn/Dropbox/Documents/Materias_doctorado/RedesComplejasBiologicos/redes2018/TPfinal/'
pathSanti = '/home/santiago/Documentos/RC/Redes2018/TPfinal/'
pathDocente = '?'


path = pathSanti

sys.path.append(path)
import Comunidades as COM


#%%

#cargo las tablas con los TCRs asociados a epitopes conocidos de la base de datos publica VDJdb
#filtro por las secuencias con Score mayor a 0 (con algun grado de confianza para su anotacion), Gene igual a TRB == T-cell Receptor Beta, Species igual a Homo Sapiens y MHC_class (clase del Complejo Mayor de Histocompatibilidad) igual a I. 
#me quedo con las columnas CDR3 (Complementary Determining Region 3) que corresponde a la secuencia de aminoacidos del TRB, Epitope que corresponde a la secuencia del epitope y Epitope_gene que corresponde al gen del cual proviene dicho epitope.


virusStr = ['YFV', 'CMV', 'EBV', 'HIV-1', 'InfAV', 'HCV']
TRBdf = {}

for i in virusStr:    
    df0 = pd.read_csv(path + 'Datos/'+ i, sep='\t')
    df1 = df0[(df0.Score > 0) &(df0.Gene == 'TRB') & (df0.Species == 'HomoSapiens') & (df0.MHC_class == 'MHCI')]
    df2 = df1 [['CDR3', 'Epitope', 'Epitope_gene', 'Epitope_species']]
    
    df3 = df2.drop_duplicates(subset = 'CDR3')
    TRBdf[i] = df3
    print(len(TRBdf[i]))
#%%% 
    
#prueba de pairwise2 de Biopython
#pairwise sequence alignment

#alignments = pairwise2.align.localds("LSPADKTNVKAA", "PEEKAVA", pMHCmatrix, -10, -1)

#local -> uses Smith Waterman local alignment algorithm
#match parameter 'd': A dictionary returns the score of any pair of characters -> BLOSUM62
#gap parameter 's': Same open and extend gap penalties for both sequences -> in this case, gap open penalty of 10 and a gap extension penalty of 1

#score_only=True prints only the best alignment score

#investigar que matriz BLOSUM usar para alignment de los CDR3 y que parametros para gap penalty

#%%%

#tomo 145 TRB de cada categoria (virus) de forma random 
#genero un dataframe conjunto con 870 TRB

TRBdf_rnd0 = []

for i in virusStr:
    TRBdf_rnd0.append(TRBdf[i].sample(n=145))

TRBdf_rnd = pd.concat(TRBdf_rnd0)
TRBdf_rnd = TRBdf_rnd.drop_duplicates(subset = 'CDR3')

#hay un solo CDR3 que esta repetido, ver!
TRBdf_rnd = TRBdf_rnd.reset_index(drop=True)

#%%

#cargo la matriz de similaridad del paper de Kim et al. (2009) (novelmatrixTCR.pdf) 

peptideMHCmatrix = pd.read_csv(path + 'newmatrix.txt', sep=' ')
peptideMHCmatrix.columns = (list(peptideMHCmatrix)[1:21]+[0])
peptideMHCmatrix = peptideMHCmatrix.iloc[:,0:20]

#para reescalear todos los valores de la matriz (dataframe)
#peptideMHCmatrix = peptideMHCmatrix.multiply(1000).astype(int)

pMHCmatrix = {}

pairsmatrix = list(itr.combinations_with_replacement(list(peptideMHCmatrix),2))
for (i,j) in pairsmatrix:
    pMHCmatrix[(i,j)] = peptideMHCmatrix[i][j]


#%%

#falta terminar ojo, no correr este chunk!
#escribo un archivo .csv de los TRB CDR3 (para luego calcular la kernel similarity measure)

#TRB_CDR3 = list(TRBdf_rnd['CDR3'])
#
#path_TRB_CDR3 = '/home/heli/Documents/Lab Immunoinformatics/epitope_CDR3.csv'
#
#with open(path_TRB_CDR3, mode='w') as TRB_file:
#    TRB_writer = csv.writer(TRB_file, delimiter='\n')
#    TRB_writer.writerow(TRB_CDR3)

#%%

#tomo todos los pares unicos de secuencias y hago pairwise alignment de las 
#mismas para obtener una score que refleje la similud de secuencia con las dos matrices
#Blosum 62 y pMHCmatrix

CDR3_pairs = list(itr.combinations(TRBdf_rnd["CDR3"], 2))

seq_scores = {}
selected_CDR3_pairs = {}

matrices = {'blosum62':blosum62, 'pMHCmatrix':pMHCmatrix}
gap_open_penalty = {'blosum62':-10, 'pMHCmatrix':-1}
gap_extension_penalty = {'blosum62':-1, 'pMHCmatrix':-0.1}
threshold = {'blosum62':45, 'pMHCmatrix':2}

for i in matrices.keys():
    scorelist = []
    selected_CDR3_list = []
    for (seq1, seq2) in CDR3_pairs:
        score = pairwise2.align.localds(seq1, seq2, matrices[i], gap_open_penalty[i], gap_extension_penalty[i], score_only=True)
        scorelist.append(score)
        if score > threshold[i]:
            selected_CDR3_list.append((seq1, seq2, {'weight': score}))
    selected_CDR3_pairs[i] = selected_CDR3_list
    seq_scores[i] = scorelist
    print(i)

#para encontrar el cutoff o threshold adecuado, ver como estan distribuidos los scores


#hay que evaluar adonde poner el cutoff para el score calculado
#por el momento lo fije en 45 que es mas que media + 1 desvio std (Blosum62)
#%%
#hago un histograma con los scores de las secuencias
plt.figure()
#figManager = plt.get_current_fig_manager()
#figManager.window.showMaximized()
sp=0
for i in matrices.keys():
    sp+=1
    plt.subplot(210+sp)
    plt.hist(seq_scores[i], bins=500)
    plt.axvline(x=threshold[i], color='k',LineWidth=0.5)
    #plt.xticks(np.arange(0, 80, step=5))
    plt.xlim((0,1.5*threshold[i]))
#    plt.ylim((0,60000))
    plt.subplots_adjust(left=None, bottom=None, right=None, top=None, wspace=None, hspace=0.3)
    plt.tick_params(axis='both', which='major', labelsize=10)
    plt.xlabel('Sequence score', fontsize=10)
    plt.ylabel('Occurrence (' + i + ')', fontsize=10)
    plt.show()
plt.savefig(path+'/Figuras/SequenceScore.pdf')
#%%
#armo el grafo con el output del pairwise alignment de secuencias
#le agrego el atrobuto a los nodos de epitope_sp es decir a que virus corresponde el epitopes al cual se une la secuencia CDR3

TRB = {}

for i in matrices:
    TRB[i] = nx.Graph()
    TRB[i].add_edges_from(selected_CDR3_pairs[i])
    for n in TRBdf_rnd['CDR3']:
        if n in list(TRB[i].nodes()):
            TRB[i].nodes[n]['epitopeSp'] = TRBdf_rnd['Epitope_species'][TRBdf_rnd['CDR3'] == n].to_string(index=False)


TRB[i].number_of_nodes()
TRB[i].number_of_edges()

#los nodos que quedaron no enlazados los descartamos

#    else:
#        TRBdf.add_node(n, epitope_sp = TRBdf_rnd['Epitope_species'][TRBdf_rnd['CDR3'] == n].to_string(index=False))


#tambien es interesante evaluar la distribucion de grado de los nodos para distintos thresholds o cutoffs
netDeg = np.array(list(TRB[i].degree()))
CDR3_netDeg = [int(x) for x in netDeg[:,1]]

np.median(CDR3_netDeg)
np.mean(CDR3_netDeg)
np.std(CDR3_netDeg)

plt.figure()
plt.hist(CDR3_netDeg, bins=60)
plt.xticks(np.arange(0, 50, step=50))
plt.tick_params(axis='both', which='major', labelsize=10)
plt.xlabel('Selected sequence degree', fontsize=15)
plt.show()
#plt.savefig(path+'/Figuras/SelectedSequenceDegree.pdf')

'''
¿Qué es este gráfico? Agregar título, nombre del eje y, etc.
'''
#%%

#plots de los grafos TRB constriudos

# use one of the edge properties to control line thickness
edgewidth = [ d['weight']/20 for (u,v,d) in TRB[i].edges(data=True)]

plt.figure()
nx.draw(TRB[i], width = edgewidth, node_size=80, font_size=20)
plt.suptitle('TRB sequence similarity network', fontsize=20)
plt.show()

color_code = {'YellowFeverVirus':'purple', 'CMV':'green', 'HIV-1':'black',  'EBV':'blue', 'InfluenzaA':'red', 'HCV':'orange'}

TRB_node_color = {}
pos = {}

for i in matrices:
    TRB_epitopes = list(nx.get_node_attributes(TRB[i], "epitopeSp").values())
    node_color = []
    for k in TRB_epitopes:
        node_color.append(color_code[k])
    TRB_node_color[i] = node_color
    pos[i] = nx.kamada_kawai_layout(TRB[i])

plt.figure()
nx.draw(TRB[i],
        width=0.1,
        edge_color = 'k',
        node_color = TRB_node_color[i], 
        node_size=50,
        font_size=10,
        with_labels=False,
        )
plt.show()

plt.figure()
nx.draw(TRB[i],
        pos[i],
        width=0.1,
        edge_color = 'k',
        node_color = TRB_node_color[i], 
        node_size=50,
        font_size=10,
        with_labels=False,
        )
YFV = mpatches.Patch(color='purple', label='Yellow Fever Virus')
CMV = mpatches.Patch(color='green', label='Citomegalovirus')
HIV = mpatches.Patch(color='black', label='HIV-1')
EBV = mpatches.Patch(color='blue', label='Epstein-Barr Virus')
InfluenzaA = mpatches.Patch(color='red', label='Influenza A Virus')
HCV = mpatches.Patch(color='orange', label='Hepatitis C Virus')
plt.legend(handles=[YFV, CMV, HIV,  EBV, InfluenzaA, HCV], fontsize=20)
plt.show()
plt.savefig(path+'/Figuras/Red.pdf')

#%%%

# Comunidades:

methods=['l','fg','im']#,'ng']
methodsStr = ['Louvain','Fast Greedy','Infomap']#,'Newman-Girvan']

Particion={}    # diccionario que contiene las particiones de la red para los cuatro métodos.

for m in methods:
    Particion[m]=COM.Communities(TRB[i],path,m)

'''
TRB[i] --> Hay que iterar i.
'''
#%%



