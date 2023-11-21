import demes
import msprime
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import os
import scipy
import scipy.stats as stats
import pandas as pd
import random
import copy
import seaborn as sns
sns.set_style('white')
sns.set_style('ticks')
sns.set_context('notebook')
import h5py
import allel; print('scikit-allel', allel.__version__)

demes_model = "Historic_bestMod_demes.yml"
graph = demes.load(demes_model)
demography = msprime.Demography.from_demes(graph)
#load best-fitting model from gadma as demes file to msprime.


#at time T add migration rate change between bay area and east cal
#to generate no change in migration simulation remove these add_migration_rate_change lines.
demography.add_migration_rate_change(time=0, rate=0.15, source="BayArea", dest="EastCal")
demography.add_migration_rate_change(time=0, rate=0.0, source="EastCal", dest="BayArea")
demography.add_migration_rate_change(time=15, rate=0.0001332, source="BayArea", dest="EastCal")
demography.add_migration_rate_change(time=15, rate=0.0001559, source="EastCal", dest="BayArea")
demography.sort_events()

#dd = demography.debug()
#print(dd)

samples = [msprime.SampleSet(num_samples=25, population="Orange", time=0, ploidy=2)] + \
		  [msprime.SampleSet(num_samples=25, population="BayArea", time=0, ploidy=2)] + \
		  [msprime.SampleSet(num_samples=25, population="EastCal", time=0, ploidy=2)] + \
		  [msprime.SampleSet(num_samples=25, population="Orange", time=50, ploidy=2)] + \
		  [msprime.SampleSet(num_samples=25, population="BayArea", time=50, ploidy=2)] + \
		  [msprime.SampleSet(num_samples=25, population="EastCal", time=50, ploidy=2)]

ts = msprime.sim_ancestry(samples=samples, demography=demography, sequence_length=100000000) #, num_replicates=60000)
#rep_ts = msprime.sim_ancestry(samples=samples, demography=demography, sequence_length=1000, num_replicates=5) #, num_replicates=1000)
mts = msprime.sim_mutations(ts, rate=4.6e-9,random_seed=12345)
with open("Msprime_windows_Migchange.vcf", "w") as vcf_file:
	mts.write_vcf(vcf_file)

vcf_file = "Msprime_windows_Migchange.vcf"
Filtered_VCF = "Msprime_windows_MigChange_filtered"
vcf_cmd = "vcftools --vcf %s --maf 0.05 --min-alleles 2 --max-alleles 2 --recode --out %s" % (vcf_file,Filtered_VCF)	
os.system(vcf_cmd)


#import filtered vcf file, convert to h5 file and then import to to allel
allel.vcf_to_hdf5("Msprime_windows_MigChange_filtered.recode.vcf", 'Msprime_windows_migChange_filtered.h5', fields='*', overwrite=True)

Simulated_capture = h5py.File('Msprime_windows_migChange_filtered.h5', mode='r')

#import sample data
Samples = pd.read_csv("Sim_sample_file.csv")


#create index for each contig and position
chrom = Simulated_capture['variants/CHROM']
pos = Simulated_capture['variants/POS']
#idx = allel.SortedMultiIndex(chrom,pos)

#make chunked variant table
vtbl = allel.VariantChunkedTable(Simulated_capture['variants'], names=['CHROM', 'POS', 'REF', 'ALT', 'QUAL'])

#make genotype table
genotypes = allel.GenotypeChunkedArray(Simulated_capture['calldata/GT'])
#get array of populations/time period
populations = ['Newport', 'BayArea']
times = Samples.Time.unique()
pop_time = Samples.Pop_Time.unique()

pos = vtbl['POS']

Fst_df = pd.DataFrame(columns=['Newport_fst_hist', 'Newport_fst_mod', 'BayArea_fst_hist', 'BayArea_fst_mod'])
counter = 0
for pop in populations:
	fresh_hist = "EastCal_Hist"
	salt_hist = "%s_Hist" % pop

	fresh_mod = "EastCal_Mod"
	salt_mod = "%s_Mod" % pop

	subpops = {
		"Fresh_Hist_pop": Samples[Samples.Pop_Time == fresh_hist].index,
		"Salt_Hist_pop": Samples[Samples.Pop_Time == salt_hist].index,
		"Fresh_Mod_pop": Samples[Samples.Pop_Time == fresh_mod].index,
		"Salt_Mod_pop": Samples[Samples.Pop_Time == salt_mod].index,
	}

	acs = genotypes.count_alleles_subpops(subpops)

	ac_freshHist = allel.AlleleCountsArray(acs["Fresh_Hist_pop"])
	ac_saltHist = allel.AlleleCountsArray(acs["Salt_Hist_pop"])
	ac_freshMod = allel.AlleleCountsArray(acs["Fresh_Mod_pop"])
	ac_saltMod = allel.AlleleCountsArray(acs["Salt_Mod_pop"])
	
	#estimate Fst and Dxy use bed file coordinates to set start and end position of window.	
	Fst_hist,windows,_ = allel.windowed_hudson_fst(pos, ac_freshHist, ac_saltHist, size=5000, step=2500)
	Fst_mod,_,_ = allel.windowed_hudson_fst(pos, ac_freshMod, ac_saltMod, size=5000, step=2500)
	#Dxy_hist,_,_,_ = allel.windowed_divergence(vtbl_chrom_pos, ac_freshHist, ac_saltHist, size=5000, step=2500)
	#Dxy_mod,_,_,_ = allel.windowed_divergence(vtbl_chrom_pos, ac_freshMod, ac_saltMod, size=5000, step=2500)
	if counter == 0:
		#Fst_df['FstWindows'] = windows
		Fst_df['Newport_fst_hist'] = Fst_hist
		Fst_df['Newport_fst_mod'] = Fst_mod
	else:
		Fst_df['BayArea_fst_hist'] = Fst_hist
		Fst_df['BayArea_fst_mod'] = Fst_mod
	counter += 1

Fst_df['Newport_deltaFst'] = Fst_df['Newport_fst_mod'] - Fst_df['Newport_fst_hist'] 
Fst_df['BayArea_deltaFst'] = Fst_df['BayArea_fst_mod'] - Fst_df['BayArea_fst_hist']
		
Fst_df.to_csv('SimulatedFst_MigChange.csv')		

FST_Data = pd.read_csv("Windowed_Fst_observed_data.csv",header=[0,1])
Af_fst_fin = FST_Data.stack(level=[0])
print(Af_fst_fin.head())

Af_fst_fin['Pop1'] = Af_fst_fin.index


temp_df= pd.DataFrame(Af_fst_fin['Pop1'].to_list(), columns=['num','Population'])

Af_fst_fin.insert(0,'Population', np.array(temp_df['Population']))
Af_fst_fin = Af_fst_fin.drop(columns='Pop1')
Af_fst_fin = Af_fst_fin.sort_values(by='Population')
Af_fst_fin.reset_index(drop=True, inplace=True)
print(Af_fst_fin.head())
populations = ['Newport', 'BayArea']
new_df = []
for pop in populations:
	pop_df = Af_fst_fin[Af_fst_fin['Population'] == pop]
	Fst_95 = pop_df.Fst_hist.quantile(0.95)
	deltaFst_med = pop_df.deltaFst.quantile(0.5)
	#label SNPs as outliers or not outliers
	pop_df["SNPtype"] = np.where((pop_df["Fst_hist"]>Fst_95), "hist_outlier_obs", "not_outlier_obs")
	print(pop,Fst_95,sep="\t")
	new_df.append(pop_df)

#remerge dataframes and create violin plot.
AF_FST_fin_df = pd.concat(new_df)
AF_FST_fin_df.reset_index(drop=True, inplace=True)

AF_FST_fin_df.to_csv("Windows_simulated_obs_Fst_data.csv")

fig, ax = plt.subplots(figsize=(12, 6))
sns.despine(ax=ax, offset=5)
ax = sns.violinplot(x='Population', y='deltaFst', hue='SNPtype', data=AF_FST_fin_df)
ax.legend(bbox_to_anchor=(1, 1), loc='upper left', fontsize='medium')
ax.axhline(0, color="black", lw=0.5)
#ax.set_ylabel(ylab)
fig.tight_layout()
fig.savefig("Fstchange_windows_095_1Jun_obs.pdf", dpi=300)
