#!/usr/bin/env python3
import numpy as np
import glob
import re
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
import seaborn as sns
sns.set_style('white')
sns.set_style('ticks')
sns.set_context('notebook')


#import file of sample names and 
SampleData = pd.read_csv("SampleLabels.txt",sep='\t')

path = "K3_final_simQfiles/" + "*_theta"
files = glob.glob(path)

df_list = []
for file in files:
	SplitFile = file.split('_')
	Gen = SplitFile[4]
	Mig = SplitFile[5]
	K = SplitFile[7]
	Ne = SplitFile[8]
	
	Gen2 = re.sub(r"Gen","", Gen)
	#Mig0.00015.geno
	Mig2 = re.sub(r"(Mig)|.geno","",Mig)
	#print(Mig2) 
	df1 = pd.DataFrame()
	df1['Population'] = SampleData['Population']
	df1['Epoch'] = SampleData['Epoch']
	df1['GenTime'] = Gen2
	df1['MigRate'] = Mig2
	df1['K'] = K
	df1['Ne'] = Ne
	

	Qfile = pd.read_csv(file, sep='\t', header=None)
	#print(Qfile.tail())
	Q1 = Qfile.idxmax(1).ix[99]	
	Q2 = Qfile.idxmax(1).ix[124]
	Q3 = Qfile.idxmax(1).ix[149]
	#print(Q1,Q2,Q3, sep=",")
	df1['Q1'] = Qfile[Q1]
	df1['Q2'] = Qfile[Q2]
	df1['Q3'] = Qfile[Q3]
	df_list.append(df1)


Sim_Dystruct_Out = pd.concat(df_list, axis=0,ignore_index=True)


#subset to include only BayArea Mod
#plot Q3 as a function of mig

BA_mod_df = Sim_Dystruct_Out[(Sim_Dystruct_Out['Population']=='BayArea') & (Sim_Dystruct_Out['Epoch']=='Mod')]
BA_mod_small = BA_mod_df[(BA_mod_df['Ne']=='small')] 
BA_mod_mid = BA_mod_df[(BA_mod_df['Ne']=='mid')] 
BA_mod_large = BA_mod_df[(BA_mod_df['Ne']=='large')] 

fig,ax = plt.subplots(figsize=(6, 5))
sns.stripplot(x='GenTime',y='Q3',hue='MigRate', edgecolor='black', order=['NA','5','10','15','20','25'], data=BA_mod_large)
ax.legend(bbox_to_anchor=(1, 1), loc='upper left', fontsize='medium')
ax.axhline(0.44, linestyle='--')
ax.set_xlabel('Generations before present of migration rate change')
ax.set_ylabel('Proportion eastern Cal ancestry in Bay Area')

fig.savefig("Eastcal_ancestry_bayarea_mid.pdf", dpi=300, bbox_inches='tight')