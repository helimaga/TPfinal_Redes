#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Dec  3 13:44:44 2018

@author: heli
"""

#running Python 2.7 in virtual environment ENV, OLGA is a python 2.7 software!

#open terminal
#run: source /home/heli/ENV/bin/activate --> ENV activated
#run: python -c "import sys; print(sys.executable)" --> copy path
#run: deactivate --> VEnvironment ENV deactivated
#paste path in Spyder Tools --> Preferences--> Python Interpreter, save changes
#reinitiate Spyder

#%%

import numpy as np
import pandas as pd
import olga.load_model as load_model
import olga.generation_probability as pgen
import olga.sequence_generation as seq_gen

#%%

path = '/home/heli/ENV/lib/python2.7/site-packages/olga/'

#Define the files for loading in generative model/data
params_file_name = path + 'default_models/human_T_beta/model_params.txt'
marginals_file_name = path + 'default_models/human_T_beta/model_marginals.txt'
V_anchor_pos_file = path + 'default_models/human_T_beta/V_gene_CDR3_anchors.csv'
J_anchor_pos_file = path + 'default_models/human_T_beta/J_gene_CDR3_anchors.csv'

#Load data
genomic_data = load_model.GenomicDataVDJ()
genomic_data.load_igor_genomic_data(params_file_name, V_anchor_pos_file, J_anchor_pos_file)
#Load model
generative_model = load_model.GenerativeModelVDJ()
generative_model.load_and_process_igor_model(marginals_file_name)

#Process model/data for pgen computation by instantiating GenerationProbabilityVDJ
pgen_model = pgen.GenerationProbabilityVDJ(generative_model, genomic_data)

#example
#calculating pgen with restriction to V, J gene usage 
pgen_model.compute_aa_CDR3_pgen('CAWSVAPDRGGYTF', 'TRBV30*01', 'TRBJ1-2*01')
#calculating pgen without restriction to V, J gene usage 
pgen_model.compute_aa_CDR3_pgen('CAWSVAPDRGGYTF')

#Process model/data for sequence generation by instantiating SequenceGenerationVDJ
seq_gen_model = seq_gen.SequenceGenerationVDJ(generative_model, genomic_data)
#Generate some random sequences
seq_gen_model.gen_rnd_prod_CDR3()
#('TGTGCCAGCAGTGAAAAAAGGCAATGGGAAAGCGGGGAGCTGTTTTTT', 'CASSEKRQWESGELFF', 27, 8)
seq_gen_model.gen_rnd_prod_CDR3()
#('TGTGCCAGCAGTTTAGTGGGAAGGGCGGGGCCCTATGGCTACACCTTC', 'CASSLVGRAGPYGYTF', 14, 1)
seq_gen_model.gen_rnd_prod_CDR3()
#('TGTGCCAGCTGGACAGGGGGCAACTACGAGCAGTACTTC', 'CASWTGGNYEQYF', 55, 13)

#%%

#genero 5000 secuencias y las guardo 

path_redes = '/home/heli/Documents/Redes/TPfinal/'

rnd_seq = []
n=5000

for i in range(n):
    seq = seq_gen_model.gen_rnd_prod_CDR3()
    rnd_seq.append(seq)

np.savetxt(path_redes + 'rnd_seq', rnd_seq, delimiter=',', fmt='%s')




