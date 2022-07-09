import os
import numpy as np
import scipy
import pandas
import matplotlib as mpl
import matplotlib.pyplot as plt
import seaborn as sns
sns.set_style('white')
sns.set_style('ticks')
sns.set_context('notebook')
import h5py
import allel; print('scikit-allel', allel.__version__)

#Import hd5 file
SAVScapture_callset = h5py.File('/Users/phbenham/Dropbox/CA_SAVS_temporalGenomics_analyses_ms/Bioinformatics_analyses/Data/SAV_CA_downsampled_CT-GA_removed_13May21.h5', mode='r')


#import list of samples and other data
samples_fn = '/Users/phbenham/Dropbox/CA_SAVS_temporalGenomics_analyses_ms/Bioinformatics_analyses/Data/CA_sample_data_forGenDiv2.csv'
samples = pandas.read_csv(samples_fn)
print(samples.head())

#import bed file with capture regions that overlap SNPs in VCFfile
bed_regions = pandas.read_csv('/Users/phbenham/Dropbox/CA_SAVS_temporalGenomics_analyses_ms/Bioinformatics_analyses/Data/CA_SAVS_CaptureRegions_27May2021.bed.csv')
print(len(bed_regions))

#create index for each contig and position
chrom = SAVScapture_callset['variants/CHROM']
pos = SAVScapture_callset['variants/POS']
idx = allel.SortedMultiIndex(chrom,pos)
print(len(idx))
#make chunked variant table
vtbl = allel.VariantChunkedTable(SAVScapture_callset['variants'], names=['CHROM', 'POS', 'REF', 'ALT', 'QUAL'])

#make genotype table
genotypes = allel.GenotypeChunkedArray(SAVScapture_callset['calldata/GT'])
#get array of populations/time period
populations = ["Orange", "Morro", "BayArea", "Humboldt"]
times = samples.SamplePeriod.unique()

#set up empty df that includes time period and population 
#population
#historic modern
#then underneath each column list stats
mi = pandas.MultiIndex.from_product((populations,times), names=["Population","Time"])
df = pandas.DataFrame(columns=mi, index=range(len(idx)))

pos = vtbl['POS']

for pop in populations:
	
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

	acs = genotypes.count_alleles_subpops(subpops)
	acu = allel.AlleleCountsArray(acs["Fresh_Hist_pop"][:] + acs["Salt_Hist_pop"][:]+ acs["Fresh_Mod_pop"][:]+ acs["Salt_Mod_pop"][:])
	flt = acu.is_segregating() & (acu.max_allele() == 1) & (acu[:, :2].min(axis=1) > 1)
	print('retaining', np.count_nonzero(flt), 'SNPs')
	
	#get allele counts for modern and historic populations

	ac_freshHist = allel.AlleleCountsArray(acs["Fresh_Hist_pop"].compress(flt, axis=0)[:, :2])
	ac_saltHist = allel.AlleleCountsArray(acs["Salt_Hist_pop"].compress(flt, axis=0)[:, :2])
	ac_freshMod = allel.AlleleCountsArray(acs["Fresh_Mod_pop"].compress(flt, axis=0)[:, :2])
	ac_saltMod = allel.AlleleCountsArray(acs["Salt_Mod_pop"].compress(flt, axis=0)[:, :2])
	
	num_hist, den_hist = allel.hudson_fst(ac_freshHist, ac_saltHist)
	hist_fst_hudson = num_hist / den_hist
	
	num_mod, den_mod = allel.hudson_fst(ac_freshMod, ac_saltMod)
	mod_fst_hudson = num_mod / den_mod
	index_range = range(np.count_nonzero(flt))

	df.loc[index_range,(pop,"Historic")] = hist_fst_hudson
	df.loc[index_range,(pop,"Modern")] = mod_fst_hudson


DivData1 = df.unstack().unstack(level=1).reset_index(level=1,drop=True).reset_index() #.reset_index(level=1,drop=True).rename_axis('Population').reset_index()

DivData1['delta_stat'] = DivData1['Modern'] - DivData1['Historic']

grouped = DivData1.groupby(['Population'])

Outfile = open("Fst_saltmarsh_divergenceStats_31May2021.txt", "w")

SNP_pos = pandas.DataFrame({'CHROM':chrom, 'START': pos})
SNP_pos['END'] = SNP_pos['START'] + 1

append_df = []
for pops in populations:
	marsh_pop = grouped.get_group(pops)
	marsh_pop.reset_index(drop=True, inplace=True)
	#marsh_pop = marsh_pop[marsh_pop.index < 33520] 

	Fst_hist = np.asarray(marsh_pop['Historic'],dtype=np.float64)
	Fst_mod =  np.asarray(marsh_pop['Modern'],dtype=np.float64)
	print(type(Fst_hist))
	Fst_hist[Fst_hist<0]=0
	Fst_max_hist = np.nanmax(Fst_hist)
	Fst_sd_hist = np.nanstd(Fst_hist)
	Fst_mean_hist = np.nanmean(Fst_hist)
	Five_x_mean_hist = Fst_mean_hist + 5*Fst_sd_hist
	print(pops,'Historic', sep="\t", file=Outfile)
	print("Historic Mean Fst:", Fst_mean_hist, "Max Fst:",Fst_max_hist,"Mean Fst + 5*sd:", Five_x_mean_hist, sep="\t", file=Outfile)
	
	Fst_mod[Fst_mod<0]=0
	Fst_max_mod = np.nanmax(Fst_mod)
	Fst_sd_mod = np.nanstd(Fst_mod)
	Fst_mean_mod = np.nanmean(Fst_mod)
	Five_x_mean_mod = Fst_mean_mod + 5*Fst_sd_mod
	print(pops,'Modern', sep="\t", file=Outfile)
	print("Modern Mean Fst:", Fst_mean_mod, "Max Fst:",Fst_max_mod,"Mean Fst + 5*sd:", Five_x_mean_mod, sep="\t", file=Outfile)
	
	#if Fst_pos is greater than threshold new col it is outlier_hist, else non_outlier
	marsh_pop["SNPtype"] = np.where(marsh_pop['Historic']>Five_x_mean_hist, "hist_outlier", "not_outlier")
	#print(marsh_pop.head())
	#concatenate all 4 dataframes together
	append_df.append(marsh_pop)
	
	

	'''
	df_out = pandas.concat([SNP_pos, marsh_pop], axis=1)
	InBed = "%s_perSNP_FstDivergence-CTGAremoved_downsampled_31May2021.bed" % pops
	OutBed = "%s_perSNP_FstDivergence-CTGAremoved_downsampled_31May2021_annot.txt" % pops
	df_out.to_csv(InBed, na_rep='NaN', sep='\t', index=False, header=False)
	os.system("bedtools intersect -wb -a FullData_WTSP_CaptureRegions_merged.bed -b %s > %s" % (InBed, OutBed))
	'''
		
Outfile.close()

SNP_df = pandas.concat(append_df)

FinalSNP_df = SNP_df.astype({'Historic': 'float64','Modern': 'float64','delta_stat': 'float64'})
print(FinalSNP_df.dtypes)

fig, ax = plt.subplots(figsize=(12, 6))
sns.despine(ax=ax, offset=5)
ax = sns.violinplot(x='Population', y='delta_stat', hue='SNPtype', data=FinalSNP_df)#, errwidth=0.25, capsize=0.1)	
ax.legend(bbox_to_anchor=(1, 1), loc='upper left', fontsize='medium')
ax.axhline(0, color="black", lw=0.5)
#ax.set_ylabel(ylab)
fig.tight_layout()
fig.savefig("HistoricModern_changeDivergence-CTGAremoved_downsampled.pdf", dpi=300)

'''


#this code creates a data frame that can then be exported as a csv file 
#you can either export the full Fst data 'Fst_data' or only the outliers 'Fst_data_fin'
Fst_data = pandas.DataFrame({'Chrom':chrom2, 'Pos': pos2, 'Fst':Fst})
End_Pos = Fst_data['Pos'] + 1
#Fst_data['Index'] = pandas.Series(range(0,len(Fst_data)))

Fst_data.insert(loc=2, column='End', value=End_Pos)

print(Fst_data.head())
Fst_data_fin = Fst_data[(Fst_data['Fst']>Five_x_mean)]
OutBed = "SAVS_Morro_PerSNP_HudsonFst_Mod_outliers.bed"
InBed1 = "Morro_nevadensis_Modern_sorted.bed"
InBed2 = "Morro_nevadensis_Modern_buffered.bed"
InBed3 = "Morro_nevadensis_Modern_AnnotatedOutliers.txt"
Fst_data_fin.to_csv(OutBed, na_rep='NaN', sep='\t', index=False, header=False)
Fst_data.to_csv("SAVS_morro_nevadensis_Modern_allSNPs_HudsonFst.csv", na_rep='NaN')

#run bedtools slop to expand selection
os.system("bedtools sort -i %s > %s" % (OutBed, InBed1))
os.system("bedtools slop -i %s -g WTSPcontigs.txt -b 100 >%s" % (InBed1, InBed2))

#then use resulting bed file for intersect
os.system("bedtools intersect -wb -a Exons_Capturefasta.bed -b %s > %s" % (InBed2, InBed3))
'''