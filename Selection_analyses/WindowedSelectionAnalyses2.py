#! /usr/bin/env python3

import numpy as np
import scipy
import scipy.stats as stats
import pandas
import random
import copy
import matplotlib as mpl
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib_venn import venn2, venn2_circles, venn2_unweighted
from matplotlib_venn import venn3, venn3_circles
sns.set_style('white')
sns.set_style('ticks')
sns.set_context('notebook')
import h5py
import os
import allel; print('scikit-allel', allel.__version__)

#convert vcf to h5 file, comment out after first use.
#VCF_path = "/Users/phbenham/Dropbox/CA_SAVS_temporalGenomics_analyses_ms/Bioinformatics_analyses/Data/SAVS_CA_Exons_21Oct2021_sorted.vcf"
#allel.vcf_to_hdf5(VCF_path, '/Users/phbenham/Dropbox/CA_SAVS_temporalGenomics_analyses_ms/Bioinformatics_analyses/Data/SAVS_CA_Exons_21Oct2021_sorted.h5', fields='*', overwrite=True)

#Import hd5 file
SAVScapture_callset = h5py.File('/Users/phbenham/Dropbox/CA_SAVS_temporalGenomics_analyses_ms/Bioinformatics_analyses/Data/SAVS_CA_Exons_21Oct2021_sorted.h5', mode='r')


#import list of samples and other data
samples_fn = '/Users/phbenham/Dropbox/CA_SAVS_temporalGenomics_analyses_ms/Bioinformatics_analyses/Data/CA_sample_data_forGenDiv2.csv'
samples = pandas.read_csv(samples_fn)
print(samples.head())

#import bed file with capture regions that overlap SNPs in VCFfile
#bed_regions = pandas.read_csv('/Users/phbenham/Dropbox/CA_SAVS_temporalGenomics_analyses_ms/Bioinformatics_analyses/Data/CA_SAVS_CaptureRegions_27May2021.bed.csv')

WTSPcontigs = open('/Users/phbenham/Dropbox/CA_SAVS_temporalGenomics_analyses_ms/Bioinformatics_analyses/Data/WTSPcontigs.txt', 'r')
Contig_dict = {}
for line in WTSPcontigs:
	line = line.strip('\n')
	splits = line.split('\t')
	Contig_dict[splits[0]] = splits[1]
	

#create index for each contig and position
chrom = SAVScapture_callset['variants/CHROM']
pos = SAVScapture_callset['variants/POS']
idx = allel.SortedMultiIndex(chrom,pos)

#make chunked variant table
vtbl = allel.VariantChunkedTable(SAVScapture_callset['variants'], names=['CHROM', 'POS', 'REF', 'ALT', 'QUAL'])

#make genotype table
genotypes = allel.GenotypeChunkedArray(SAVScapture_callset['calldata/GT'])
#get array of populations/time period
populations = ["Orange","Morro", "BayArea", "Humboldt"]
times = samples.SamplePeriod.unique()


pos = vtbl['POS']

ac = genotypes.count_alleles()

chrom = np.unique(vtbl['CHROM'])
counter = 0

df_pop_append = []
for pop in populations:
	df_append = []
	for item in chrom:
		if int(Contig_dict[item]) > 5000:
			index_length = int(Contig_dict[item])/5000
			print(item)
			#print(int(index_length))
			#subset vtbl and genotypes based on chromosome
			loc_chrom = idx.locate_range(item)
			vtbl_chrom = vtbl[loc_chrom]
			vtbl_chrom_pos = vtbl_chrom['POS']
			#print(vtbl_chrom_pos)
			gt_chrom = genotypes[loc_chrom]
			counter +=1
		
			if len(vtbl_chrom_pos) > 1:
				#set up empty df that includes time period and population 
				#population
				#historic modern
				#then underneath each column list stats

	
				fresh_hist = "Fresh_Historic"
				salt_hist = "%s_Historic" % pop

				fresh_mod = "Fresh_Modern"
				salt_mod = "%s_Modern" % pop

				subpops = {
					"Fresh_Hist_pop": samples[samples.Time_Pop == fresh_hist].index,
					"Salt_Hist_pop": samples[samples.Time_Pop == salt_hist].index,
					"Fresh_Mod_pop": samples[samples.Time_Pop == fresh_mod].index,
					"Salt_Mod_pop": samples[samples.Time_Pop == salt_mod].index,
				}

				acs = gt_chrom.count_alleles_subpops(subpops)
		
				ac_freshHist = allel.AlleleCountsArray(acs["Fresh_Hist_pop"])
				ac_saltHist = allel.AlleleCountsArray(acs["Salt_Hist_pop"])
				ac_freshMod = allel.AlleleCountsArray(acs["Fresh_Mod_pop"])
				ac_saltMod = allel.AlleleCountsArray(acs["Salt_Mod_pop"])

				#estimate Fst and Dxy use bed file coordinates to set start and end position of window.	
				Fst_hist,windows,_ = allel.windowed_hudson_fst(vtbl_chrom_pos, ac_freshHist, ac_saltHist, size=5000, step=2500)
				Fst_mod,_,_ = allel.windowed_hudson_fst(vtbl_chrom_pos, ac_freshMod, ac_saltMod, size=5000, step=2500)
		
				Dxy_hist,_,_,_ = allel.windowed_divergence(vtbl_chrom_pos, ac_freshHist, ac_saltHist, size=5000, step=2500)
				Dxy_mod,_,_,_ = allel.windowed_divergence(vtbl_chrom_pos, ac_freshMod, ac_saltMod, size=5000, step=2500)
			
				#estimate delta pi
				pi_hist_salt,_,_,_ = allel.windowed_diversity(vtbl_chrom_pos, ac_saltHist, size=5000, step=2500)
				pi_mod_salt,_,_,_ = allel.windowed_diversity(vtbl_chrom_pos, ac_saltMod, size=5000, step=2500)
			
				pi_hist_fresh,_,_,_ = allel.windowed_diversity(vtbl_chrom_pos, ac_freshHist, size=5000, step=2500)
				pi_mod_fresh,_,_,_ = allel.windowed_diversity(vtbl_chrom_pos, ac_freshMod, size=5000, step=2500)
			
			
				#estimate delta Tajimas D
				TajD_salt_hist,_,_ = allel.windowed_tajima_d(vtbl_chrom_pos, ac_saltHist, size=5000, step=2500, min_sites=1)
				TajD_fresh_hist,_,_ = allel.windowed_tajima_d(vtbl_chrom_pos, ac_freshHist, size=5000, step=2500, min_sites=1)
			
				TajD_salt_mod,_,_ = allel.windowed_tajima_d(vtbl_chrom_pos, ac_saltMod,  size=5000, step=2500, min_sites=1)
				TajD_fresh_mod,_,_ = allel.windowed_tajima_d(vtbl_chrom_pos, ac_freshMod, size=5000, step=2500, min_sites=1)

				delta_hist = TajD_salt_hist - TajD_fresh_hist
				delta_mod = TajD_salt_mod - TajD_fresh_mod 
				dTajD_hist = (delta_hist - np.mean(delta_hist)) / np.std(delta_hist)
				dTajD_mod = (delta_mod - np.mean(delta_mod)) / np.std(delta_mod)
				
				#append to pandas dataframe
				index_range = range(len(Fst_hist))
				
				#initialize data frame with chrom, start, end position for each window
				Chrom = np.full((len(Fst_hist)), item)
				df = pandas.DataFrame(windows, columns=["Start","End"])
				df.insert(loc=0, column="Chrom", value=Chrom)
				
				#add stats to dataframe
				df["Fst_mod"] = Fst_mod
				df["Dxy_mod"] = Dxy_mod
				df["TajD_mod_salt"] = TajD_salt_mod
				df["TajD_mod_fresh"] = TajD_fresh_mod
				df["delta_TajD_mod"] = dTajD_mod
				df["Pi_mod_fresh"] = pi_mod_fresh
				df["Pi_mod_salt"] = pi_mod_salt
				df["Fst_hist"] = Fst_hist
				df["Dxy_hist"] = Dxy_hist
				df["TajD_hist_salt"] = TajD_salt_hist
				df["TajD_hist_fresh"] = TajD_fresh_hist
				df["Pi_hist_fresh"] = pi_hist_fresh
				df["Pi_hist_salt"] = pi_hist_salt
				df["delta_TajD_hist"] = dTajD_hist
				df["delta_pi_hist"] = df["Pi_hist_salt"] - df["Pi_hist_fresh"]
				df["delta_pi_mod"] = df["Pi_mod_salt"] - df["Pi_mod_fresh"]

			df_append.append(df)
		
	df_pop = pandas.concat(df_append)
	df_pop_fin = df_pop.dropna(subset=["Fst_mod","Fst_hist"])
	df_pop_fin.reset_index(drop=True, inplace=True)	
	print(df_pop_fin.head())
	OutFileName = "%s_WindowedDivergenceStats_5kb_22October2021_exons.csv" % pop
	df_pop_fin.to_csv(OutFileName , na_rep='NaN')


populations = ["Orange","Morro", "BayArea", "Humboldt"]
df_pop_append = []
df_Fst_append = []
df_Dxy_append = []
for pop in populations:
	InFileName = "%s_WindowedDivergenceStats_5kb_22October2021_exons.csv" % pop
	DivData = pandas.read_csv(InFileName) #df_pop_fin
	DivData["window_name"] = DivData["Chrom"]+':'+DivData["Start"].astype(str)+'-'+DivData["End"].astype(str)
	print(DivData.head())

	#get 95% thresholds for Dxy, Fst, deltaPi[lowest 5%] in historic group.

	Fst_95 = DivData.Fst_hist.quantile(0.99)
	Dxy_95 = DivData.Dxy_hist.quantile(0.95)
	Pi_95 = DivData.delta_pi_hist.quantile(0.05)
	#"TajD_hist_salt"
	print(Fst_95,Dxy_95,Pi_95, sep='\t')

	Fst_outlier = DivData[DivData["Fst_hist"]>Fst_95]
	Dxy_outlier = DivData[DivData["Dxy_hist"]>Dxy_95]
	Pi_outlier = DivData[DivData["delta_pi_hist"]<Pi_95]
	print(len(Fst_outlier),len(Dxy_outlier),len(Pi_outlier))
	
	#FstOut = "SAVS_%s_Fst_outliers_7Dec21.csv" % pop
	#DxyOut = "SAVS_%s_Dxy_outliers_7Dec21.csv" % pop

	#Fst_outlier.to_csv(FstOut, na_rep='NaN', index=False)
	#Dxy_outlier.to_csv(DxyOut, na_rep='NaN', index=False)
	
	OutBed = "%s_Fst_outlierWindows_7Dec.bed" % pop
	InBed1 = "%s_Fst_outlierWindows_7Dec_sorted.bed" % pop
	InBed2 = "%s_Fst_outlierWindows_7Dec_AnnotatedOutliers.txt" % pop
	Fst_outlier.to_csv(OutBed, na_rep='NaN', sep='\t', index=False, header=False)

	#run bedtools slop to expand selection
	os.system("bedtools sort -i %s > %s" % (OutBed, InBed1))

	#then use resulting bed file for intersect
	os.system("bedtools intersect -wb -a FullData_WTSP_CaptureRegions_merged.bed -b %s > %s" % (InBed1, InBed2))
	
	DXYBed = "%s_Dxy_outlierWindows_7Dec.bed" % pop
	Dxy_inBed1 = "%s_Dxy_outlierWindows_7Dec_sorted.bed" % pop
	Dxy_inBed2 = "%s_Dxy_outlierWindows_7Dec_AnnotatedOutliers.txt" % pop
	Dxy_outlier.to_csv(DXYBed, na_rep='NaN', sep='\t', index=False, header=False)

	#run bedtools slop to expand selection
	os.system("bedtools sort -i %s > %s" % (DXYBed, Dxy_inBed1))

	#then use resulting bed file for intersect
	os.system("bedtools intersect -wb -a FullData_WTSP_CaptureRegions_merged.bed -b %s > %s" % (Dxy_inBed1, Dxy_inBed2))
	
	

	#add column describing whether snp is an outlier or not
	DivData["SNPtype"] = np.where((DivData["Fst_hist"]>Fst_95) & (DivData["Dxy_hist"]>Dxy_95) & (DivData["delta_pi_hist"]<Pi_95), "hist_outlier", "not_outlier")
	
	#add column with population information
	pop_name = np.full((len(DivData)), pop)
	DivData["Population"] = pop_name 
	DivData["delta_Fst"] = DivData["Fst_mod"] - DivData["Fst_hist"]

	#perform permutation test on differences between means of outlier and non outlier
	SNPtype_mean = DivData.groupby(["SNPtype"])["delta_Fst"].mean()
	#gT = np.abs(np.average(feat_vir[:,0]) - np.average(feat_ver[:,0]))
	gT = np.abs(SNPtype_mean[0]-SNPtype_mean[1])
	deltaFst = DivData["delta_Fst"]
	print(gT)
	
	#Copy pooled distribution:
	pS = copy.copy(deltaFst)
	#Initialize permutation:
	pD = []
	#Define p (number of permutations):
	p=100
	# Permutation loop:
	for i in range(0,p):
		print(i)
		# Shuffle the data:
		random.shuffle(pS)
		# Compute permuted absolute difference of your two sampled distributions and store it in pD
		pD.append(np.abs(np.average(pS[0:int(len(pS)/2)]) - np.average(pS[int(len(pS)/2):])))
		
	p_val = len(np.where(pD>=gT)[0])/p
	
	print(pop, p_val)
	

	Fst_windows = Fst_outlier["window_name"]
	Dxy_windows = Dxy_outlier["window_name"]
	Pi_windows = Pi_outlier["window_name"]
	
	#subset Fst_outliers for Dxy_95>0.95 and Pi_95 <0.05
	Overlapping_Outliers = DivData[DivData["SNPtype"]=="hist_outlier"]
	print(pop, len(Overlapping_Outliers))
	

	df_Fst_append.append(Fst_outlier)
	df_Dxy_append.append(Dxy_outlier)
	

	fig, ax = plt.subplots(figsize=(12, 6))
	ax = venn3([set(Fst_windows),set(Dxy_windows),set(Pi_windows)], set_labels = ("Fst outliers", "Dxy outliers", "Pi outliers"))
	fig.tight_layout()
	FigName = "%s_venn_outliers_Fst_Dxy_pi.pdf" % pop
	fig.savefig(FigName, dpi=300)

	OutBed = "SAVS_%s_Fst_outliers_17Nov21.bed" % pop
	#InBed1 = "%s_nevadensis_Modern_sorted.bed" % pop
	#InBed2 = "%s_nevadensis_Modern_buffered.bed" % pop
	InBed3 = "%s_nevadensis_Modern_FstAnnotatedOutliers.txt" % pop
	Fst_outlier.to_csv(OutBed, na_rep='NaN', sep='\t', index=False, header=False)
	#then use resulting bed file for intersect
	os.system("bedtools intersect -wb -a FullData_WTSP_CaptureRegions_merged.bed -b %s > %s" % (OutBed, InBed3))

	df_pop_append.append(DivData)
	

df_Fst = pandas.concat(df_Fst_append)
df_Fst.reset_index(drop=True, inplace=True)

df_Dxy = pandas.concat(df_Dxy_append)
df_Dxy.reset_index(drop=True, inplace=True)

df_Fst_fin = df_Fst.drop_duplicates('window_name', keep='first')
df_Dxy_fin = df_Dxy.drop_duplicates('window_name', keep='first')

df_Fst_fin.reset_index(drop=True, inplace=True)
df_Dxy_fin.reset_index(drop=True, inplace=True)

df_Fst_fin = df_Fst_fin[["Chrom", "Start", "End", "Fst_hist", "Dxy_hist", "window_name"]]
df_Dxy_fin = df_Dxy_fin[["Chrom", "Start", "End", "Fst_hist", "Dxy_hist", "window_name"]]

print(df_Fst_fin.head())
print(df_Dxy_fin.head())

OutBed = "SAVS_Fst_outlierWindows_17Nov.bed"
InBed1 = "SAVS_Fst_outlierWindows_17Nov_sorted.bed"
InBed2 = "SAVS_Fst_outlierWindows_17Nov_Modern_AnnotatedOutliers.txt"
df_Fst_fin.to_csv(OutBed, na_rep='NaN', sep='\t', index=False, header=False)

#run bedtools slop to expand selection
os.system("bedtools sort -i %s > %s" % (OutBed, InBed1))

#then use resulting bed file for intersect
os.system("bedtools intersect -wb -a FullData_WTSP_CaptureRegions_merged.bed -b %s > %s" % (InBed1, InBed2))


OutBed_dxy = "SAVS_Dxy_outlierWindows_17Nov.bed"
InBed1_dxy = "SAVS_Dxy_outlierWindows_17Nov_sorted.bed"
InBed2_dxy = "SAVS_Dxy_outlierWindows_17Nov_Modern_AnnotatedOutliers.txt"
df_Dxy_fin.to_csv(OutBed_dxy, na_rep='NaN', sep='\t', index=False, header=False)

#run bedtools slop to expand selection
os.system("bedtools sort -i %s > %s" % (OutBed_dxy, InBed1_dxy))

#then use resulting bed file for intersect
os.system("bedtools intersect -wb -a FullData_WTSP_CaptureRegions_merged.bed -b %s > %s" % (InBed1_dxy, InBed2_dxy))


OutFileName = "AllPopulations_WindowedDivergenceStats_5kb_25October2021_exons.csv"
df_full.to_csv(OutFileName , na_rep='NaN')

fig, ax = plt.subplots(figsize=(12, 6))
sns.despine(ax=ax, offset=5)
ax = sns.violinplot(x='Population', y='delta_Fst', hue='SNPtype', data=df_full)#, errwidth=0.25, capsize=0.1)	
ax.legend(bbox_to_anchor=(1, 1), loc='upper left', fontsize='medium')
ax.axhline(0, color="black", lw=0.5)
#ax.set_ylabel(ylab)
fig.tight_layout()
fig.savefig("HistoricModern_changeDivergence-5kbwindows_25Oct2021.pdf", dpi=300)
