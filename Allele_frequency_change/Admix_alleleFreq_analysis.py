#!/usr/bin/env python3
import os
import sys
import pandas as pd
import numpy as np
from scipy import stats
import matplotlib as mpl
import matplotlib.pyplot as plt
import seaborn as sns
sns.set_style('white')
sns.set_style('ticks')
sns.set_context('notebook')


Usage="""This script parses data on admixture proportions from DyStruct (K=3) and gets mean
ancestry proportions for each time period and population. 
Second, creates text file from each population with individuals in order to use vcftools to estimate
allele frequencies for each population. 
Get outlier non-outlier Fst estimates too?
"""

#import csv file with sample data and admixture proportions.
VCF_path = "/Users/phbenham/Dropbox/CA_SAVS_temporalGenomics_analyses_ms/Bioinformatics_analyses/Data/SAVS_CA_downsampled_CT-GAremoved_13May21.vcf"
AdmixData = pd.read_csv("Admixture_Prop_forAF_analysis.csv")
print(AdmixData.head())

#if directory for individual files does not exist make file.
if not os.path.exists("Population_individual_datafiles"):
		os.mkdir("Population_individual_datafiles")

if not os.path.exists("VCFtools_AF_output"):
		os.mkdir("VCFtools_AF_output")

if not os.path.exists("VCFtools_Fst_outlier_output"):
		os.mkdir("VCFtools_Fst_outlier_output")					

#for loop to go over each population, produce text file, means for admixture proportion.
Pop_times = AdmixData.Region_time.unique()

Admix_mean = AdmixData.groupby(["Region_time"])["Nev_anc"].mean()

Admix_groups = AdmixData.groupby("Region_time")["SampleID"].apply(list).to_dict()




for key in Admix_groups:
	OutFileName = "Population_individual_datafiles/%s_individuals.txt" % key
	OutFile = open(OutFileName, 'w')
	for item in Admix_groups[key]:
		print(item, file=OutFile)
	OutFile.close()	

	#AlleleFreq_Out = "VCFtools_AF_output/%s_AlleleFreq_out" % key
	#VCFtools_cmd = "vcftools --vcf %s --keep %s --freq --out %s" % (VCF_path, OutFileName, AlleleFreq_Out)
	#os.system(VCFtools_cmd)

populations = ['BayArea', 'Cbeldingi', 'Humboldt', 'Morro']
for pop in populations:
	Hist_out = "VCFtools_Fst_outlier_output/%s_Fst_historical" % pop 
	Mod_out =  "VCFtools_Fst_outlier_output/%s_Fst_modern" % pop
	VCFtools_fst_hist_cmd = "vcftools --vcf %s --weir-fst-pop Population_individual_datafiles/Nevadensis_Historic_individuals.txt --weir-fst-pop Population_individual_datafiles/%s_Historic_individuals.txt --out %s" % (VCF_path,pop,Hist_out)
	VCFtools_fst_mod_cmd = "vcftools --vcf %s --weir-fst-pop Population_individual_datafiles/Nevadensis_Modern_individuals.txt --weir-fst-pop Population_individual_datafiles/%s_Modern_individuals.txt --out %s" % (VCF_path,pop,Mod_out)
	os.system(VCFtools_fst_hist_cmd)
	os.system(VCFtools_fst_mod_cmd)


AF_data = pd.read_csv("Alt_alleleFreq_Fst_Data.csv")

AF_data_new = AF_data.dropna(thresh=13)
AF_data_new.reset_index(drop=True, inplace=True)


#Nev_hist_freq	Nev_mod_freq	Hum_hist_freq	Hum_mod_freq	BA_hist_freq	BA_mod_freq	Morro_hist_freq	Morro_mod_freq	Bel_hist_freq	Bel_mod_freq
populations = ['Hum','BA', 'Morro', 'Bel']

for pop in populations:
	#calculate mean freq hist
	MeanFreqCol = "%s_meanAF_hist" % pop
	Salt_pop = "%s_hist_freq" % pop
	HistDiff = "%s_hist_diff" % pop
	AF_data_new[MeanFreqCol] = (AF_data_new[Salt_pop] + AF_data_new["Nev_hist_freq"])/2
	AF_data_new[HistDiff] = abs(AF_data_new["Nev_hist_freq"] - AF_data_new[Salt_pop])
	
	#calculate historical difference between historic mean AF and historic salt AF (historic deviation)
	HistDev = "%s_hist_dev" % pop
	AF_data_new[HistDev] = AF_data_new[MeanFreqCol] - AF_data_new[Salt_pop]
	
	#calculate modern difference between historic mean AF and salt modern salt AF (modern deviation)
	ModDev = "%s_Mod_dev" % pop
	Salt_pop_mod = "%s_mod_freq" % pop
	ModDiff = "%s_mod_diff" % pop
	Salt_Pop_Change = "%s_temp_change" % pop
	AF_data_new[ModDev] = abs(AF_data_new[MeanFreqCol] - AF_data_new[Salt_pop_mod])
	AF_data_new[ModDiff] = abs(AF_data_new["Nev_hist_freq"] - AF_data_new[Salt_pop_mod])
	AF_data_new[Salt_Pop_Change] = AF_data_new[Salt_pop_mod] - AF_data_new[Salt_pop]
	
	#calculate difference between modern and historic deviation
	delta_Dev = "%s_delta_Dev" % pop
	delta_Diff = "%s_delta_Diff" % pop
	AF_data_new[delta_Dev] = AF_data_new[ModDev] - AF_data_new[HistDev]
	AF_data_new[delta_Diff] = AF_data_new[ModDiff] - AF_data_new[HistDev]
	
	
	#estimate delta_fst for each snp
	Hist_Fst = "%s_hist_fst" % pop
	Mod_Fst = "%s_mod_fst" % pop
	delta_Fst = "%s_delta_fst" % pop
	AF_data_new[delta_Fst] = AF_data_new[Mod_Fst] - AF_data_new[Hist_Fst]


#AF calculations for Nev population.
AF_data_new["Nev_hist_dev"] = abs(AF_data_new["BA_meanAF_hist"] - AF_data_new["Nev_hist_freq"])
AF_data_new["Nev_Mod_dev"] = abs(AF_data_new["BA_meanAF_hist"] - AF_data_new["Nev_mod_freq"])
AF_data_new["Nev_delta_Dev"] = AF_data_new["Nev_Mod_dev"] - AF_data_new["Nev_hist_dev"]
	
print(AF_data_new)
AF_data_new.to_csv("AlleleFreqData_final_14Dec2021.csv")



AF_FST_DATA = pd.read_csv("AlleleFreqData_final_15Dec2021.csv",header=[0,1])

Af_fst_fin = AF_FST_DATA.stack(level=[0])
print(Af_fst_fin.head())

Af_fst_fin['Pop1'] = Af_fst_fin.index


temp_df= pd.DataFrame(Af_fst_fin['Pop1'].to_list(), columns=['num','Population'])

Af_fst_fin.insert(0,'Population', np.array(temp_df['Population']))
Af_fst_fin = Af_fst_fin.drop('Pop1',1)
Af_fst_fin = Af_fst_fin.sort_values(by='Population')
Af_fst_fin.reset_index(drop=True, inplace=True)
print(Af_fst_fin.head())

Pop_fst = Af_fst_fin.groupby('Population').hist_fst.quantile(0.99)
print(Pop_fst)

#split dataframe by population
#for each population estimate Fst 0.95 threshold
populations = ['Hum','BA', 'Morro', 'Bel']

new_df = []
for pop in populations:
	pop_df = Af_fst_fin[Af_fst_fin['Population'] == pop]
	Fst_99 = pop_df.hist_fst.quantile(0.99)
	#label SNPs as outliers or not outliers
	pop_df["SNPtype"] = np.where((pop_df["hist_fst"]>Fst_99), "hist_outlier", "not_outlier")
	
	new_df.append(pop_df)

#remerge dataframes and create violin plot.
AF_FST_fin_df = pd.concat(new_df)
AF_FST_fin_df.reset_index(drop=True, inplace=True)

AF_FST_fin_df.to_csv("Outlier_Fst_data.csv")

fig, ax = plt.subplots(figsize=(12, 6))
sns.despine(ax=ax, offset=5)
ax = sns.violinplot(x='Population', y='delta_fst', hue='SNPtype', data=AF_FST_fin_df)
ax.legend(bbox_to_anchor=(1, 1), loc='upper left', fontsize='medium')
ax.axhline(0, color="black", lw=0.5)
#ax.set_ylabel(ylab)
fig.tight_layout()
fig.savefig("Fstchange_saltmarshpops_099_16Dec.pdf", dpi=300)

Hum_perm = []
BA_perm = []
Morro_perm = []
Bel_perm = []
Outlier_corr = []
Not_Outlier_corr = []
for pop in populations:
	Pop_outlier_df = AF_FST_fin_df[(AF_FST_fin_df["Population"] == pop) & (AF_FST_fin_df["SNPtype"] == "hist_outlier")]

	data1 = Pop_outlier_df['hist_dev']
	data2 = Pop_outlier_df['temp_change']
	Outlier_corrP,_ = stats.pearsonr(data1, data2)
	#corrS,_ = stats.spearmanr(data1, data2)
	
	print(pop, len(data1), Outlier_corrP, sep="\t")
	Outlier_corr.append(Outlier_corrP)
	#subset non-outliers 1000 times that are the same length as hist_outliers
	#for each subset estimate Pearson's correlation coefficient and add to list
	Pop_notOutlier_df = AF_FST_fin_df[(AF_FST_fin_df["Population"] == pop) & (AF_FST_fin_df["SNPtype"] == "not_outlier")]
	data1_not = Pop_notOutlier_df['hist_dev']
	data2_not = Pop_notOutlier_df['temp_change']
	NotOutlier_corrP,_ = stats.pearsonr(data1_not, data2_not)
	Not_Outlier_corr.append(NotOutlier_corrP)
	
	# Permutation loop:
	p=1000
	for i in range(0,p):
		# Shuffle the data:
		NumRows = len(data1)
		Random_df = Pop_notOutlier_df.sample(n=NumRows)
		
		# Compute correlation coefficient
		data1_rand = Random_df['hist_dev']
		data2_rand = Random_df['temp_change']
		NonOutlier_corrP,_ = stats.pearsonr(data1_rand, data2_rand)
		
		if pop == 'Hum':
			Hum_perm.append(NonOutlier_corrP)
		elif pop == "BA":
			BA_perm.append(NonOutlier_corrP)
		elif pop == "Morro":
			Morro_perm.append(NonOutlier_corrP)	
		elif pop == "Bel":
			Bel_perm.append(NonOutlier_corrP)				
		
	zero = np.float64(0.0)
	if pop == 'Hum':
		p_val_0 = 1 - (len(np.where(Hum_perm>=zero)[0])/p)
		p_val_out = 1 - (len(np.where(Hum_perm<=Outlier_corrP)[0])/p)
		print(pop,p_val_0, p_val_out)
	elif pop == "BA":
		p_val_0 = 1 - (len(np.where(BA_perm>=zero)[0])/p)
		p_val_out = 1 - (len(np.where(BA_perm<=Outlier_corrP)[0])/p)
		print(pop, p_val_0,p_val_out)
	elif pop == "Morro":
		p_val_0 = 1 - (len(np.where(Morro_perm>=zero)[0])/p)
		p_val_out = 1 - (len(np.where(Morro_perm<=Outlier_corrP)[0])/p)
		print(pop, p_val_0,p_val_out)
	elif pop == "Bel":
		p_val_0 = 1 - (len(np.where(Bel_perm>=zero)[0])/p)
		p_val_out = 1 -  (len(np.where(Bel_perm<=Outlier_corrP)[0])/p)
		print(pop, p_val_0,p_val_out)	

CorrResults = pd.DataFrame(list(zip(Hum_perm, BA_perm,Morro_perm,Bel_perm)),columns =populations)

CorrResults2 = CorrResults.melt(var_name="Population", value_name="Correlation")
#print(CorrResults2)
fig, ax = plt.subplots(figsize=(7, 5))
sns.despine(ax=ax, offset=5)
sns.boxplot(x="Population", y="Correlation", color="gray", order= ['Bel','Morro','BA','Hum'], data=CorrResults2)
ax.plot(3, Outlier_corr[0], color="red", marker='o')
ax.plot(2, Outlier_corr[1], color="red", marker='o')
ax.plot(1, Outlier_corr[2], color="red", marker='o')
ax.plot(0, Outlier_corr[3], color="red", marker='o')
ax.axhline(0, color="black", lw=0.5)
fig.tight_layout()
fig.savefig("Correlation_099_plot.pdf", dpi=300)

Admix_change = [0.64,0.48,0.17,0]

admix_df = pd.DataFrame(list(zip(populations, Outlier_corr, Admix_change, Not_Outlier_corr)), columns = ['Population', 'Corr', 'Admix','Not_corr'])

print(admix_df)
fig, ax = plt.subplots(figsize=(7, 5))
fig = sns.lmplot(x='Admix', y='Corr', data=admix_df)
fig.tight_layout()
fig.savefig("Admxiture_correlation_plot.pdf", dpi=300)

admix_lm = stats.linregress(Admix_change, Outlier_corr)
print(admix_lm)
'''
#fig, ax = plt.subplots(figsize=(7, 5))
fig = sns.lmplot(x='hist_dev', y='temp_change', hue='Population', data=NoNev_df)
#fig.axhline(0, color="black", lw=0.5, ls='--')
#fig.axvline(0, color="black", lw=0.5, ls='--')
fig.tight_layout()
fig.savefig("AFchange_ALL_notout_saltmarshpops_14Dec.pdf", dpi=300)


	

AF_outlier_df = NoNev_df[(NoNev_df["SNPtype"] == "not_outlier")] #(AF_FST_fin_df["SNPtype"] == "not_outlier") & 
#AF_NEV_fin_df = AF_FST_fin_df[(AF_FST_fin_df['Population'] == 'Nev') ]
Median_dev = AF_outlier_df.groupby('Population')['delta_Diff'].median()
print(Median_dev)

pop_dev = AF_outlier_df.groupby('Population')['delta_Diff'].apply(np.array)


F_value,p_value = stats.f_oneway(pop_dev[0],pop_dev[1],pop_dev[2],pop_dev[3])
print(p_value)
#ttest = scipy.stats.ttest_ind(AF_BA_fin_df['delta_Dev'],AF_NEV_fin_df['delta_Dev'])
#print(ttest)
'''
#plot delta deviation as a function of historic Fst 
Pop_outlier_df = AF_FST_fin_df[(AF_FST_fin_df["Population"] == "Bel")]
Pop_outlier_df[Pop_outlier_df['hist_fst']<0] = 0
fig, ax = plt.subplots(figsize=(7, 5))
sns.despine(ax=ax, offset=5)
sns.histplot(x='hist_fst', bins=25, data=Pop_outlier_df)
ax.axvline(0.37, color="red", lw=0.5, ls="--")
ax.axvline(0.73, color="black", lw=0.5, ls="--")
ax.axvline(0.7, color="black", lw=0.5, ls="--")
#ax.axvline(Median_dev[3], color="green", lw=0.5)
fig.savefig("Bel_perSNP_FstDivergence_histogram.pdf",  dpi=300)
