import numpy as np
import os
import scipy
import scipy.stats as stats
import copy
import random
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
import seaborn as sns
sns.set_style('white')
sns.set_style('ticks')
sns.set_context('notebook')

Final_dat = pd.read_csv("AllPopulations_WindowedDeltaFst_final_7Dec2021.csv")

#plot figure comparing delta Fst between historical and modern datasets
fig, ax = plt.subplots(figsize=(12, 6))
sns.despine(ax=ax, offset=5)
ax = sns.violinplot(x='Population', y='delta_Fst', hue='SNPtype', data=Final_dat, order=["Orange", "Morro", "BayArea", "Humboldt"])	
ax.legend(bbox_to_anchor=(1, 1), loc='upper left', fontsize='medium')
ax.axhline(0, color="black", lw=0.5)
#ax.set_ylabel(ylab)
fig.tight_layout()
fig.savefig("HistoricModern_deltaFst_7Dec2021.pdf", dpi=300)

#perform permutation test to determine whether outliers declined in frequency more than expected by chance. 
populations = ["Orange","Morro", "BayArea", "Humboldt"]
OutFileName = "PermutationAnalysis_results_7Dec21.txt"
OutFile = open(OutFileName, 'a') 
for pop in populations:
	Pop_df = Final_dat[Final_dat["Population"]==pop]
	Pop_df.reset_index(drop=True, inplace=True)
	#print(Pop_df.head())
	
	SNPtype_mean = Pop_df.groupby(["SNPtype"])["delta_Fst"].mean()
	#gT = np.abs(np.average(feat_vir[:,0]) - np.average(feat_ver[:,0]))
	gT = np.abs(SNPtype_mean[0]-SNPtype_mean[1])
	deltaFst = Pop_df["delta_Fst"]

	#Copy pooled distribution:
	pS = copy.copy(deltaFst)
	#Initialize permutation:
	pD = []
	#Define p (number of permutations):
	p=10000
	# Permutation loop:
	for i in range(0,p):
		# Shuffle the data:
		random.shuffle(pS)
		# Compute permuted absolute difference of your two sampled distributions and store it in pD
		pD.append(np.abs(np.average(pS[0:int(len(pS)/2)]) - np.average(pS[int(len(pS)/2):])))
	
	p_val = len(np.where(pD>=gT)[0])/p
	pd_df = pd.DataFrame(pD, columns=['pD'])

	pd_mean = np.mean(pD)
	pd_max = np.max(pD)
	print("Population:", pop, "Difference:", gT,"p-value:", p_val,"perm. mean:", pd_mean,"perm. max", pd_max, sep='\t', file=OutFile)	
	
	#make histogram of output data with actual mean difference plotted.
	fig, ax = plt.subplots(figsize=(7, 5))
	sns.despine(ax=ax, offset=5)
	ax = sns.histplot(data=pd_df, x="pD", kde=True)
	ax.set_xlabel("Permuted difference")
	ax.set_ylabel("Proportion")
	ax.axvline(gT, color = "red")
	fig.tight_layout()
	FigName = "%s_PermTestResults_7Dec21.pdf" % pop
	fig.savefig(FigName, dpi=300)	