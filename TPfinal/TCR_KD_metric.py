#!/usr/bin/env python3
# -*- coding: utf-8 -*

#%% Importar modulos

import pandas as pd
import os 

#%%

path = os.path.dirname(os.path.realpath('__file__'))

#%%

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

    #me quedo con las columnas CDR3 (Complementary Determining Region 3) que corresponde a la secuencia de aminoacidos del TRB, Epitope que corresponde a la secuencia del epitope y Epitope_gene que corresponde al gen del cual proviene dicho epitope.
    df2 = df1 [['CDR3']]
    if nSample!=-1:
            df2 = df2.sample(n=nSample)
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

# Reseteo de indices, nombrarlos de 0 a N...
TRBdf = TRBdf.reset_index(drop=True)

#%%

#escribo TRBdf en un .csv 

TRBdf.to_csv(path_or_buf = path + '/KD_epitopes', columns=['CDR3'], header=False, index=False)
    
###corro en la consola seq2score_db_kernel y obtengo KD_scores

#%%

#recupero los scores

KD_scores = pd.read_csv(path + '/KD_scores', sep=' ', header=None)
KD_scores = KD_scores.rename({0: 'score_cat' , 1: 'seq1', 2: 'seq2', 3:'score'}, axis='columns')

KD_scores = KD_scores.set_index(['seq1', 'seq2'])
KD_scores = KD_scores.sort_index()

seq1 = 'CASSEATGASYEQYF' 
seq2 = 'CASSEGGQAYNEQFF'

score = KD_scores.loc[(seq1, seq2), 'score'][0]








