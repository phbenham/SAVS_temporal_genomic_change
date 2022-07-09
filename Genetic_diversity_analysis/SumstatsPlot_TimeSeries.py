#! /usr/bin/env python3

import numpy as np
import scipy
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
import seaborn as sns
#sns.set_style('white')
#sns.set_style('ticks')
#sns.set_context('notebook')
import h5py
import allel; print('scikit-allel', allel.__version__)

##########################################################################################
def SumStat_barplot(data,ylab,filename):
	#fig, ax = plt.subplots(figsize=(12, 4))
	#sns.despine(ax=ax, offset=5)
	#x = sns.displot(data, y='value', hue='Time', kind="kde",multiple="stack")#data=data, errwidth=0.5, capsize=0.1)
	#ax = sns.barplot(x='TimePoint',y='value',hue='DataType', data=data, errwidth=0.5, capsize=0.1)	
	#ax.legend(bbox_to_anchor=(1, 1), loc='upper left', fontsize='medium')
	fig = sns.catplot(x='Decade',y='value',col='DataType', color="gray",data=data, kind="bar", ci=95, n_boot=1000, errwidth=0.5, capsize=0.1)		
	fig.set_axis_labels("Decade", ylab)
	#ax.set_ylabel(ylab)
	#fig.tight_layout()
	fig.savefig(filename, dpi=300)
##########################################################################################

SumStats = pd.read_csv('CaptureRegions_SumStats_CTGAremoved_downsampled_SFBay_TimeSeries_1Sept2021.csv', header=[0,1,2,3], index_col=0)
#print(SumStats.head())

SumStats2 = SumStats.melt()

#grouped = SumStats1.groupby(['Population'])
#dfs = [grouped.get_group('Beldingi'),grouped.get_group('Morro'),grouped.get_group('SFBay'),grouped.get_group('SanPablo'),grouped.get_group('Suisun'),grouped.get_group('Humboldt'),grouped.get_group('brooksi'),grouped.get_group('nevadensis')]
#SumStats2 = pd.concat(dfs)

#SegSites = SumStats1[SumStats1['Statistic'].isin(['S'])]
NucDiv = SumStats2[SumStats2['Statistic'].isin(['pi'])]
WattTheta = SumStats2[SumStats2['Statistic'].isin(['wattersonTheta'])]
TajD = SumStats2[SumStats2['Statistic'].isin(['tajimaD'])]

#SegSites['value'] = SegSites['value'].astype('float64')
NucDiv['value'] = NucDiv['value'].astype('float64')
WattTheta['value'] = WattTheta['value'].astype('float64')
TajD['value'] = TajD['value'].astype('float64')

data = [NucDiv, WattTheta, TajD]
ylab = ["Mean nucleotide diversity of capture regions", "Mean Watterson's theta of capture regions", "Mean Tajima's D of capture regions"]
filename = ["NucDiv_plot_SFBay_TimeSeries_1Sept2021.pdf", "WattTheta_plot_SFBay_TimeSeries_1Sept2021.pdf", "TajD_plot_SFBay_TimeSeries_1Sept2021.pdf"]

for i in range(3):
	SumStat_barplot(data[i],ylab[i],filename[i])

SumStats2.to_csv("SFBay_timeseries_data_melt.csv")