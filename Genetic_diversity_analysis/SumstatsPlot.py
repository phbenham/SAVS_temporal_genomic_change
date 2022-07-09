#! /usr/bin/env python3

import numpy as np
import scipy
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
import seaborn as sns
sns.set_style('white')
sns.set_style('ticks')
sns.set_context('notebook')
import h5py
import allel; print('scikit-allel', allel.__version__)

##########################################################################################
def SumStat_barplot(data,ylab,filename):
	fig, ax = plt.subplots(figsize=(14, 6))
	sns.despine(ax=ax, offset=5)
	#ax = sns.violinplot(x="Population", y="F",hue='DataType', data=data)#, errwidth=0.25, capsize=0.1)
	#x = sns.catplot(x="Population", y='F', hue='SamplePeriod', col="DataType", kind="violin",data=data)
	ax = sns.barplot(x='Population', y='delta_stat', hue='DataType', data=data, errwidth=0.25, capsize=0.1)	
	ax.legend(bbox_to_anchor=(1, 1), loc='upper left', fontsize='medium')
	ax.axhline(0, color="black", lw=0.5)
	ax.set_ylabel(ylab)
	fig.tight_layout()
	fig.savefig(filename, dpi=300)
##########################################################################################
'''
HetVcf = pd.read_csv("Het_data_may2021/Het_DownsampledCTGAremoved.csv")
print(HetVcf.head())

grouped = HetVcf.groupby(['Population'])
dfs = [grouped.get_group('Beldingi'),grouped.get_group('MorroBay'),grouped.get_group('BayArea'),grouped.get_group('NWCal'),grouped.get_group('EastCal')] 
HetVcf1 = pd.concat(dfs)

#HetVcf1["het_obs"] = 1 - (HetVcf1["Obs_HOM"]/HetVcf1["N_SITES"])

grouped2 = HetVcf1.groupby(['SamplePeriod'])
dfs = [grouped2.get_group('Historic'),grouped2.get_group('Twenties'),grouped2.get_group('Middle'),grouped2.get_group('Modern')]#,grouped2.get_group('1930s'),grouped2.get_group('1950s'),grouped2.get_group('1960s'),grouped2.get_group('2010s')] 
HetVcf2 = pd.concat(dfs)

SumStat_barplot(HetVcf2, "Inbreeding_coefficient", "Het_data_may2021/Cal_het_exon_nongen.pdf")

'''
SumStats = pd.read_csv('CaptureRegions_SumStats_CTGAremoved_downsampled_final_27Aug2021.csv', header=[0,1,2,3], index_col=0)
print(SumStats.head())

Sumstats_melt = SumStats.melt()
print(Sumstats_melt.head())
Sumstats_melt.to_csv("CA_diversityStats_5pop_1Sept2021.csv", na_rep='NaN')
 
SumStats1 = SumStats.unstack().unstack(level=1).reset_index(level=3,drop=True).reset_index() #.reset_index(level=1,drop=True).rename_axis('Population').reset_index()
print(SumStats1.head())

SumStats1['delta_stat'] = SumStats1['Modern'] - SumStats1['Historic']


grouped = SumStats1.groupby(['Population'])
dfs = [grouped.get_group('NewportBay'),grouped.get_group('MorroBay'),grouped.get_group('BayArea'),grouped.get_group('NWCal'),grouped.get_group('EastCal')]#,grouped.get_group('Humboldt'),grouped.get_group('brooksi'),grouped.get_group('NECal')]
SumStats3 = pd.concat(dfs)



#grouped = SumStats2.groupby(['DataType'])
#dfs = [grouped.get_group('FullData'),grouped.get_group('FullDownsampled'),grouped.get_group('CT_GAmut_removed'),grouped.get_group('Downsampled_CTGA_removed')]
#SumStats3 = pd.concat(dfs)

#SegSites = SumStats1[SumStats1['Statistic'].isin(['S'])]
NucDiv = SumStats3[SumStats3['Statistic'].isin(['pi'])]
WattTheta = SumStats3[SumStats3['Statistic'].isin(['wattersonTheta'])]
TajD = SumStats3[SumStats3['Statistic'].isin(['tajimaD'])]

#SegSites['value'] = SegSites['value'].astype('float64')
NucDiv['delta_stat'] = NucDiv['delta_stat'].astype('float64')
WattTheta['delta_stat'] = WattTheta['delta_stat'].astype('float64')
TajD['delta_stat'] = TajD['delta_stat'].astype('float64')

print(TajD.head())

data = [NucDiv, WattTheta, TajD]
ylab = ["Change in nucleotide diversity of capture regions", "Change in Watterson's theta of capture regions", "Change in Tajima's D of capture regions"]
filename = ["NucDiv_deltaComparison_plot_1Sept2021.pdf", "WattTheta_deltaComparison_plot_1Sept2021.pdf", "TajD_deltaComparison_plot_1Sept2021.pdf"]

for i in range(3):
	SumStat_barplot(data[i],ylab[i],filename[i])
