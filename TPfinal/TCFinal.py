#!/usr/bin/env python3
# -*- coding: utf-8 -*-

#%% Importar modulos

import pandas as pd
import itertools as itr
import networkx as nx
import matplotlib.pyplot as plt
import numpy as np
from Bio import pairwise2
import community
import scipy.stats as stats
import sys
import os
import matplotlib.patches as mpatches
import pickle
from collections import Counter


#%% Determinar path en el cual esta TCFinal
#por primera y unica vez correr con F5 (y luego frenar corrida) (de lo contrario no funciona!)

path = os.path.dirname(os.path.realpath('__file__'))

#deberia arrojar el directorio del archivo TCFinal 

#como chequeo correr en la consola:
#os.getcwd()
#esto deberia arrojar a path como el current working directory

sys.path.append(path)

#%% Importar modulos relevantes para la deteccion de comunidades 

import Comunidades as COM
import GrafoLouvain as GL

#%% OPCIONES DE ACHICAMIENTO DE RED

nSample = -1 # si -1, no quita ninguno...
virusSelect = [0,1,2,3,4,5]

#%% CARGAR DATOS VDJdb

virusStr = ['YFV', 'CMV', 'EBV', 'HIV-1', 'InfAV', 'HCV']
virusStr = list( virusStr[i] for i in virusSelect)

TRBVirDf = {} # Diccionario de DataFrames

for i in virusStr:

    # Cargo las tablas con los TCRs asociados a epitopes asociados a cada especie de virus en la lista
    df0 = pd.read_csv(path + '/Datos/'+ i, sep='\t')

    #filtro por las secuencias con Score mayor a 0 (con algun grado de confianza para su anotacion), Gene igual a TRB == T-cell Receptor Beta, Species igual a Homo Sapiens y MHC_class (clase del Complejo Mayor de Histocompatibilidad) igual a I. 
    df1 = df0[(df0.Score > 0) &(df0.Gene == 'TRB') & (df0.Species == 'HomoSapiens') & (df0.MHC_class == 'MHCI')] # 
    df1 = df1.drop_duplicates('CDR3')
    #me quedo con las columnas CDR3 (Complementary Determining Region 3) que corresponde a la secuencia de aminoacidos del TRB, Epitope que corresponde a la secuencia del epitope y Epitope_gene que corresponde al gen del cual proviene dicho epitope.
    df2 = df1 [['CDR3', 'V', 'J', 'Epitope', 'Epitope_gene', 'Epitope_species']]
    if nSample!=-1:
            df2 = df2.sample(n=nSample, random_state=457)
    print(len(df2))

#    # Elimino secuencias CDR3 duplicadas
#    df3 = df2.drop_duplicates(subset = 'CDR3')

    TRBVirDf[i] = df2
#%% Armar un unico DataFrame con todxs lxs virus :P

TRBdf0 = []

for i in virusStr:
    TRBdf0.append(TRBVirDf[i]) # Lista de dataframes

# Dataframe
TRBdf = pd.concat(TRBdf0)

#veo cuales aun estan duplicados luego de concatenar todos los dataframes
TRBdf[TRBdf.duplicated(subset='CDR3', keep=False)].sort_values('CDR3')
#son 10!

#arbitrariamente me quedo con el que aparece primero en la tabla de TRBdf
TRBdf = TRBdf.drop_duplicates('CDR3', keep='first')

# Reseteo de indices, nombrarlos de 0 a N...
TRBdf = TRBdf.reset_index(drop=True)

#analizo lo que quedo en el TRBdf final

TRBdf['Epitope'].nunique()
TRBdf['Epitope_gene'].nunique()
TRBdf['Epitope_species'].nunique()

table_s = TRBdf.groupby(['Epitope_species']).count()
table_s['Epitope'].sort_index()

table_gs = TRBdf.groupby(['Epitope_species','Epitope_gene']).count()
table_gs['Epitope'].sort_index()

table_e=TRBdf.groupby(['Epitope_species', 'Epitope_gene','Epitope']).count()

table_e_aux=TRBdf.groupby(['Epitope', 'Epitope_gene']).count()
table_e_aux.sort_index()

table_e.to_csv(path_or_buf=path+'/table_epitopes')

#%% Cargar METRICAS

# peptide-MHC matrix (para epitopes o peptides unicamente!)

# Cargo la matriz de similaridad del paper de Kim et al. (2009) (novelmatrixTCR.pdf) 
peptideMHCmatrix = pd.read_csv(path + '/newmatrix.txt', sep=' ')
peptideMHCmatrix.columns = (list(peptideMHCmatrix)[1:21]+[0])
peptideMHCmatrix = peptideMHCmatrix.iloc[:,0:20]

#para reescalear todos los valores de la matriz (dataframe)
#peptideMHCmatrix = peptideMHCmatrix.multiply(1000).astype(int)

pMHCmatrix = {}

pairsmatrix = list(itr.combinations_with_replacement(list(peptideMHCmatrix),2))
for (i,j) in pairsmatrix:
    pMHCmatrix[(i,j)] = peptideMHCmatrix[i][j]

#metricas para los CDR3 de los TCRbetas

from Bio.SubsMat.MatrixInfo import blosum62
from pt_tcr_distances import *
from vj_distances import vdist,jdist,levenshteinDistance
from trm_distances import compareTrimer
from bm_distances import compareDimer
from profile_distances import profile_distance_allprop

#%% ELIMINO DATAFRAMES AL PEDO... ASI al apretar en la columna de "Type" en el 
# variable explorer los tienen todos arriba y pipí cucú
del(df0,df1,df2,TRBVirDf,TRBdf0, table_e, table_e_aux, table_gs, table_s, peptideMHCmatrix)
#%% DEFINIR REDES

columns = ['nodeType', 'matrixName', 'matrix', 'gop', 'gep', 'computOpt']

dfNetworks = pd.DataFrame([\
                           ['CDR3','blosum62',blosum62,-10,-1,'pairwiseAlign'],\
                           ['CDR3', 'distancePT' , 'none' ,'-','-','distancePT'],\
                           ['CDR3', 'distanceProfile' , 'none' ,'-','-','distanceProfile'],\
#                           ['CDR3', 'compareLength' , 'none' ,'-','-','compareLength'],\
                           ['CDR3', 'distanceTrimer' , 'none' ,'-','-','distanceTrimer'],\
                           ['CDR3', 'distanceDimer' , 'none' ,'-','-','distanceDimer'],\
#                           ['CDR3', 'none' , 'none' ,'-','-','distanceVJ'],\
                           ['CDR3', 'distanceLvsh' , 'none' ,'-','-','distanceLvsh'],\
                           ['CDR3', 'distanceKernel' , 'none' ,'-','-','distanceKernel'],\
                           ['Epitope','blosum62',blosum62,-10,-1,'pairwiseAlign'],\
                           ['Epitope','pMHCmatrix',pMHCmatrix,-1,-0.1,'pairwiseAlign'],\
                           ], columns=columns)

nNets = dfNetworks.index

#%%

#cargar dataframe con los resultados de la metrica Kernel SeqSimilarity ('distanceKernel')

KD_scores = pd.read_csv(path + '/KD_scores', sep=' ', header=None)
KD_scores = KD_scores.rename({0: 'score_cat' , 1: 'seq1', 2: 'seq2', 3:'score'}, axis='columns')

KD_scores = KD_scores.set_index(['seq1', 'seq2'])
KD_scores = KD_scores.sort_index()

#%% COMPUTAR SCORES!
#tomo todos los pares unicos de secuencias y hago pairwise alignment de las 
#mismas u otro metodo para obtener una score

nodesNets = {}
scoresNets = {}
nodePairsNet = {}
for net in nNets:
    print(net)
    # 1: Elijo los nodos
    nodos = list(TRBdf[dfNetworks.loc[net,'nodeType']])
    # 2: Elimino nodos repetidos
    nodos = list(set(nodos))
    # 3: Los ordeno alfabeticamente
    nodos = sorted(nodos)
    
    nodesNets[net] = nodos
    
    nodePairs = list(itr.combinations(nodos, 2))
    nodePairsNet[net] = nodePairs 
    scorelist = []
    
    matrix = dfNetworks.loc[net,'matrix']
    gop = dfNetworks.loc[net,'gop']
    gep = dfNetworks.loc[net,'gep'] 
    computOpt = dfNetworks.loc[net,'computOpt'] 
    
    TRBdf = TRBdf.set_index('CDR3')
    
    for (seq1, seq2) in nodePairs:
        if computOpt == 'pairwiseAlign':
            score = pairwise2.align.localds(seq1, seq2, matrix, gop, gep, score_only=True)
        elif computOpt == 'distancePT':
            score = -weighted_cdr3_distance(seq1, seq2, default_distance_params)
        elif computOpt == 'distanceProfile':
            score = -profile_distance_allprop(seq1, seq2)
#        elif computOpt == 'compareLength':
#            score = -abs(len(seq1) - len(seq2))
        elif computOpt == 'distanceTrimer':
            score = -compareTrimer(seq1, seq2)
        elif computOpt == 'distanceDimer':
            score = -compareDimer(seq1, seq2)
#        elif computOpt == 'distanceVJ':
#            v1, j1  = TRBdf.loc[seq1, ['V', 'J']]
#            v2, j2  = TRBdf.loc[seq2, ['V', 'J']]
#           score = vdist(v1, v2) + jdist(j1, j2)
        elif computOpt == 'distanceLvsh':
            score = -levenshteinDistance(seq1, seq2)
        elif computOpt == 'distanceKernel':
            score = KD_scores.loc[(seq1, seq2), 'score'][0]
        scorelist.append(score)
    scoresNets[net] = scorelist
    TRBdf = TRBdf.reset_index()

#%%

#guardo los diccionarios con los scores de la red completa (todos los nodos posibles)
#check os.getcwd() before saving object


def save_obj(obj, name):
    with open(name + '.pkl', 'wb') as f:
        pickle.dump(obj, f, pickle.HIGHEST_PROTOCOL)

save_obj(nodesNets, "nodesNets")
save_obj(scoresNets, "scoresNets")
save_obj(nodePairsNet, "nodePairsNet")

def load_obj(name):
    with open(name + '.pkl', 'rb') as f:
        return pickle.load(f)


#%% DEFINIR ENLACES A PARTIR DE UMBRAL

perScore = {'CDR3': 0.9995, 'Epitope': 0.95}

cdr3Flag=True
epiFlag=True
posCDR3={}
posEpi={}

edgesNets = {}
threshNets = {}
graphsNets = {}

for net in nNets:
    print(net)
    threshNets[net] = np.percentile(scoresNets[net],100*perScore[dfNetworks.loc[net,'nodeType']])
    selectedPairs = []
    pair = -1
    for (seq1, seq2) in nodePairsNet[net]:
        pair+=1
        if scoresNets[net][pair] > threshNets[net]:
            selectedPairs.append((seq1, seq2))
    edgesNets[net] = selectedPairs
    
    graphsNets[net] = nx.Graph()
    graphsNets[net].add_edges_from(edgesNets[net])
    graphsNets[net].add_nodes_from(nodesNets[net])
    if cdr3Flag and dfNetworks.loc[net,'nodeType']=='CDR3' and net==6:
        posCDR3=nx.kamada_kawai_layout(graphsNets[net])
        #net=distanceKernel
        cdr3Flag=False
    if epiFlag and dfNetworks.loc[net,'nodeType']=='Epitope' and net==8:
        posEpi=nx.kamada_kawai_layout(graphsNets[net])
        #net=pMHCmatrix
        epiFlag=False

for net in nNets:
    print(net)
    #cobertura
    totalNodes = len(graphsNets[net].nodes)
    giantNodes = len(max(nx.connected_component_subgraphs(graphsNets[net]), key=len))
    coverage =  (giantNodes/totalNodes)*100
    print(coverage)
    #retencion
    isolateNodes = len(list(nx.isolates(graphsNets[net])))
    retention =  ((totalNodes-isolateNodes)/totalNodes)*100
    print(retention)
    print(giantNodes)
    
#tomo la componente gigante!   
for net in nNets:
    graphsNets[net] = max(nx.connected_component_subgraphs(graphsNets[net]), key=len)    
    
#%% HISTOGRAMAS SCORES + UMBRAL

#hago un histograma con los scores de las secuencias
f = plt.figure()
#figManager = plt.get_current_fig_manager()
#figManager.window.showMaximized()

sp=0
for net in nNets:
    sp+=1
    plt.figure(net)
    plt.hist(scoresNets[net], bins=int(np.floor(np.sqrt(len(scoresNets[net])))))
    plt.axvline(x=threshNets[net], color='k',LineWidth=0.5)
    #plt.xticks(np.arange(0, 80, step=5))
    #plt.xlim((0,1.5*threshNets[net]))
    titStr = dfNetworks.loc[net,'nodeType'] + '; ' + dfNetworks.loc[net,'matrixName']
    plt.title(titStr, fontsize=20)
    plt.tick_params(axis='both', which='major', labelsize=10)
    if net>1:
        plt.xlabel('Sequence score', fontsize=20)
    plt.ylabel('Occurrence', fontsize=20)
    plt.show()

#%%
#Grafico GRAFOS con distintas metricas
#usar layout Kamada Kawai cuando se toma giant component

color_code = {'YellowFeverVirus':'purple', 'CMV':'green', 'HIV-1':'black',  'EBV':'blue', 'InfluenzaA':'red', 'HCV':'orange'}

label_flag=False

for net in nNets:
    print(net)
    
    colNodeList = []
    for n in graphsNets[net].nodes():
        colNodeList.append(color_code[list(TRBdf.loc[TRBdf[dfNetworks.loc[net,'nodeType']]==n,'Epitope_species'])[0]])
        
    labels_gene = {}
    if dfNetworks.loc[net,'nodeType']=='Epitope':     
        for n in graphsNets[net].nodes():
           labels_gene[n]=list(TRBdf.loc[TRBdf[dfNetworks.loc[net,'nodeType']]==n,'Epitope_gene'])[0]
           label_flag=True

    titStr = dfNetworks.loc[net,'nodeType'] + '; ' + dfNetworks.loc[net,'matrixName']
    pos = nx.kamada_kawai_layout(graphsNets[net])
    
    plt.figure()
    nx.draw(graphsNets[net],
            pos,
            width=0.2,
            edge_color = 'k',
            labels = labels_gene,
            node_color= colNodeList,
            node_size=200,
            font_size=20,
            with_labels=label_flag,
           )
   
    plt.suptitle(titStr, fontsize=30)

    YFV = mpatches.Patch(color='purple', label='Yellow Fever Virus')
    CMV = mpatches.Patch(color='green', label='Citomegalovirus')
    HIV = mpatches.Patch(color='black', label='HIV-1')
    EBV = mpatches.Patch(color='blue', label='Epstein-Barr Virus')
    InfluenzaA = mpatches.Patch(color='red', label='Influenza A Virus')
    HCV = mpatches.Patch(color='orange', label='Hepatitis C Virus')
    plt.legend(handles=[YFV, CMV, HIV,  EBV, InfluenzaA, HCV], fontsize=20)
   
#%% GRAFOS a posicion fija

color_code = {'YellowFeverVirus':'purple', 'CMV':'green', 'HIV-1':'black',  'EBV':'blue', 'InfluenzaA':'red', 'HCV':'orange'}

posCDR3Flag=True
posEpiFlag=True
del(ax)
for net in nNets:
    print(net)
    
    colNodeList = []
    for n in graphsNets[net].nodes():
        colNodeList.append(color_code[list(TRBdf.loc[TRBdf[dfNetworks.loc[net,'nodeType']]==n,'Epitope_species'])[0]])

    
    if dfNetworks.loc[net,'nodeType']=='CDR3':
        if posCDR3Flag:
            plt.figure(0)
            posCDR3Flag=False
            sp=0
        pos0 = posCDR3
        sp+=1
        ax = plt.subplot(240+sp)
    elif dfNetworks.loc[net,'nodeType']=='Epitope':
        if posEpiFlag:
            plt.figure(1)
            posEpiFlag=False
            sp=0
        pos0 = posEpi
        sp+=1
        ax = plt.subplot(120+sp)
    pos={}
    for n in graphsNets[net].nodes():
            pos[n]=pos0[n]


    metrica = dfNetworks.loc[net,'matrixName']
    
    nx.draw(graphsNets[net],
            pos,
            width=0.2,
            edge_color = 'k',
            node_color= colNodeList, 
            node_size=50,
            font_size=10,
            with_labels=False,
           )
    plt.xlim(-1,1)
    plt.ylim(-1,1)
    ax.set_title('Net: ' + metrica)
    if sp==1:
        YFV = mpatches.Patch(color='purple', label='Yellow Fever Virus')
        CMV = mpatches.Patch(color='green', label='Citomegalovirus')
        HIV = mpatches.Patch(color='black', label='HIV-1')
        EBV = mpatches.Patch(color='blue', label='Epstein-Barr Virus')
        InfluenzaA = mpatches.Patch(color='red', label='Influenza A Virus')
        HCV = mpatches.Patch(color='orange', label='Hepatitis C Virus')
        plt.legend(handles=[YFV, CMV, HIV,  EBV, InfluenzaA, HCV], fontsize=20)
#%% CLUSTERIZACION: INFOMAP, LOUVAIN, ...

clusterType = {'Infomap':{'method':'TC3','acr':'im'},\
               'Louvain':{'method':'TC3','acr':'l'},}

cMethods = list(clusterType.keys())

clusterNets = pd.DataFrame()
clusterNets = clusterNets.astype('object')
idx = -1
for net in nNets:
    for clust in cMethods:
        print(str(net) + ' ' + clust)
        idx+=1
        clusterNets.loc[idx,'IDNet'] = int(net)
        clusterNets.loc[idx,'clusterType'] = clust
        if clusterType[clust]['method'] == 'TC3':
            methodStr = clusterType[clust]['acr']
            clusterNets.loc[idx,'clusterDict'] = [COM.Communities(graphsNets[net],path,methodStr)] # Va como lista, no como diccionario...
        commPartition = clusterNets.loc[idx,'clusterDict']
        clusterNets.loc[idx,'silhouette'] = [COM.silhouetteJuancho(graphsNets[net],commPartition,'all')]
        clusterNets.loc[idx,'silhouetteAvg'] = COM.silhouetteJuancho(graphsNets[net],commPartition,'mean')
        clusterNets.loc[idx,'modularity'] = community.modularity(commPartition,graphsNets[net])
clusterNets['IDNet'] = clusterNets['IDNet'].astype(int)

#%% GRAFOS para InfoMap y Louvain

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

markerComm = ['*','.','v','^','3','4','8','+','X','H']
#colors = list(dict(mcolors.BASE_COLORS, **mcolors.CSS4_COLORS).keys())
#colores = list(range(len(colors)))
#random.shuffle(colores)

for net in nNets:
    print(net)
    commPart = [\
    list(clusterNets.loc[(clusterNets['IDNet']==net) & (clusterNets['clusterType']==cMethods[0]),'clusterDict'])[0],\
    list(clusterNets.loc[(clusterNets['IDNet']==net) & (clusterNets['clusterType']==cMethods[1]),'clusterDict'])[0]]
    A,B = COM.dictsValues2Mat(commPart[0],commPart[1])
    commPart[0] = A
    commPart[1] = B
    
    pos = nx.kamada_kawai_layout(graphsNets[net])
        
    plt.figure(net)
    sp = -1
    
    netType = dfNetworks.loc[net,'nodeType']
    metrica = dfNetworks.loc[net,'matrixName']
    
    for clust in cMethods:
        sp+=1
        colNodeList = []
        for n in graphsNets[net].nodes():
            colNodeList.append(colores[commPart[sp][n]])
        ax = plt.subplot(121+sp)
        nx.draw(graphsNets[net],
                pos,
                width=0.2,
                edge_color = 'k',
                node_color= colNodeList, 
                node_size=50,
                font_size=10,
                with_labels=False,
               )
        ax.set_title('Net: ' + netType + '-' + metrica + '\n' + 'Partition: ' + clust)

#%% SILUETAS para InfoMap y Louvain
for net in nNets:
    commPart = [\
    list(clusterNets.loc[(clusterNets['IDNet']==net) & (clusterNets['clusterType']==cMethods[0]),'clusterDict'])[0],\
    list(clusterNets.loc[(clusterNets['IDNet']==net) & (clusterNets['clusterType']==cMethods[1]),'clusterDict'])[0]]
    A,B = COM.dictsValues2Mat(commPart[0],commPart[1])
    commPart[0] = A
    commPart[1] = B

    plt.figure(net)
    sp = -1
    
    netType = dfNetworks.loc[net,'nodeType']
    metrica = dfNetworks.loc[net,'matrixName']
    for clust in cMethods:
        sp+=1

        dfsillpart = pd.DataFrame()
        silueta0 = list(clusterNets.loc[(clusterNets['IDNet']==net)&(clusterNets['clusterType']==clust),'silhouette'])[0]
        silmed = float(clusterNets.loc[(clusterNets['IDNet']==net)&(clusterNets['clusterType']==clust),'silhouetteAvg'])
        print(silmed)
        for n in graphsNets[net].nodes():
            dfsillpart.loc[n,'community'] = commPart[sp][n]
            dfsillpart.loc[n,'color'] = colores[commPart[sp][n]]
            dfsillpart.loc[n,'silhouette'] = silueta0[n]
        dfsillpart = dfsillpart.sort_values(by=['community', 'silhouette'], ascending=[True, False])
        
        silueta1 = list(reversed(dfsillpart['silhouette']))
        color1 = list(reversed(dfsillpart['color']))
        nodos1 = list(reversed(dfsillpart.index))
        
        
        color2 = []
        for n in graphsNets[net].nodes():
            color2.append(dfsillpart.loc[n,'color'])
        
        ax = plt.subplot(121+sp)
        
        plt.barh(range(len(nodos1)), silueta1, align='center',color=color1)
        plt.axvline(x=silmed,color='k',linestyle='--')
        plt.axvline(x=0,color='k',linestyle='-',linewidth=0.5)
        plt.text(silmed+0.01,1,'mean= ' + str(round(silmed*100)/100),fontsize=12,color='k')
        plt.yticks(range(len(nodos1)),'')
        plt.xlabel('silueta')
        plt.xlim(-1,1)
        ax.set_title('Net: ' + netType + '-' + metrica + '\n' + 'Partition: ' + clust)
        plt.show()
        plt.draw()
del(dfsillpart)
#%% TESTS de FISHER 2x2 (SIGUE los IDS del DF clusterNets)

# Pueden elegir el criterio separador que les parezca acá (columnas de TRBdf)
sepCrit = 'Epitope_species'
pThresh = 0.025 # p Muy exigente

contMatrices = []
fisherTests = []

idx=-1
for net in nNets:
    sepDict = {}
    comPart = {}
    N = len(graphsNets[net].nodes())
    for n in graphsNets[net].nodes():
        sepDict[n] = list(TRBdf.loc[TRBdf[dfNetworks.loc[net,'nodeType']]==n,sepCrit])[0]
    comPart = list(clusterNets.loc[(clusterNets['IDNet']==net) & (clusterNets['clusterType']==cMethods[0]),'clusterDict'])[0]
    
    sepCritVal = list(sepDict.values())
    comPartVal = list(comPart.values())
    sepCritList = list(set(sepCritVal))
    comPartList = list(set(comPartVal))

    for clust in cMethods:
        idx=idx+1
        contMatrix = pd.DataFrame(columns=sepCritList)
        fisherTest = pd.DataFrame(columns=sepCritList)
        for com in comPartList:
            A = comPartVal.count(com)
            for sepval in sepCritList:
                Nsepval = sepCritVal.count(sepval)
                Asepval = 0
                for n in graphsNets[net].nodes():
                    if (sepDict[n] == sepval) & (comPart[n] == com):
                        Asepval+=1
                contMatrix.loc[com,sepval] = Asepval
                oddsratio, pvalue = stats.fisher_exact([[Asepval,A-Asepval], [Nsepval-Asepval,N-A-(Nsepval-Asepval)]],'less')
                if pvalue > 1-pThresh:
                    fisherTest.loc[com,sepval] = '+' + sepval
                else:
                    fisherTest.loc[com,sepval] = 'Neutro'
        contMatrices.append(contMatrix)
        fisherTests.append(fisherTest)

        print(dfNetworks.loc[net,'matrixName'])
        print(dfNetworks.loc[net,'nodeType'])
        print(clusterNets.loc[idx,'clusterType'])
        print(contMatrix)
        print(fisherTest)
        print('\n\n\n')
        
#%%

#calculo la pureza y la consistencia de los clusters definidos por cada metrica y metodo de clustering
        
clusterNets['Purity_clusters'] = clusterNets['Purity_clusters'].astype(object)
clusterNets['Consistency_clusters'] = clusterNets['Consistency_clusters'].astype(object)

for n in range(len(clusterNets)):
    print(n)
    net = clusterNets.loc[n, 'IDNet']
    if dfNetworks.loc[net,'nodeType']=='CDR3':
        print('CDR3')
        dCDR3=clusterNets['clusterDict'][n]
        mx=max(dCDR3.values())+1
        x = [[] for i in range(mx)]
        x_epi = []
        x_cluster = []
        consistency=[]
        w_consistency = []
        purity=[]
        w_epi = 0
        
        for key, value in dCDR3.items():
            epitope = TRBdf.loc[TRBdf['CDR3']==key, 'Epitope'].values[0]
            x[value].append(epitope)
            x_epi.append(epitope)
            x_cluster.append(value)
            
        for j in range(mx):
            count_elements=Counter(x[j])
            max_element = count_elements.most_common()[0][1]
            purity.append((max_element/len(x[j]))*100)

        pur_avg = sum(purity)/len(purity)
        
        for k in set(x_epi):
            zcluster = []
            for iz, z in enumerate(x_epi):
                if z == k:
                    zcluster.append(x_cluster[iz])
            c_zcluster=Counter(zcluster)
            max_zcluster=c_zcluster.most_common()[0][1]
            if len(zcluster) > 1:
                consistency.append((max_zcluster/len(zcluster))*100)
                w_consistency.append(max_zcluster*100)
                w_epi += len(zcluster)
        
        con_avg = sum(consistency)/len(consistency)
        w_con_avg = sum(w_consistency)/w_epi
        
    else:
        purity = 'Nan'
        pur_avg = 'Nan'
        consistency = 'Nan'
        con_avg = 'Nan'
    
    clusterNets.at[n, 'Purity_clusters']=purity
    clusterNets.at[n, 'Purity_avg']=pur_avg
    clusterNets.at[n, 'Consistency_clusters']=consistency
    clusterNets.at[n, 'Consistency_avg']=con_avg
    clusterNets.at[n, 'Consistency_w_avg']=w_con_avg
  
#%%

# Comparación entre las redes.

# Primero, se reduce la red grande:
GrafoResumido=GL.GrafoLouvain(X)    # X es el grafo a reducir

distancia=nx.algorithms.similarity.graph_edit_distance(GrafoResumido,Y) # Y es el  grafo con el que lo vamos a comparar una vez reducido.

# FALTA REEMPLAZAR 'X' E 'Y' ARRIBA.










