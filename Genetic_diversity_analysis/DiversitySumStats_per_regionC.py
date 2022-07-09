#! /usr/bin/env python3

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


#VCF_path = "SAVS_CA_downsampled_CT-GAremoved_13May21.vcf"
#allel.vcf_to_hdf5(VCF_path, 'SAVS_CA_downsampled_CT-GAremoved_13May21.h5', fields='*', overwrite=True)


#import list of samples and other data
samples_fn = 'CA_sample_data_forGenDiv4.csv'
samples = pandas.read_csv(samples_fn)
print(samples.head())

#import bed file with capture regions
bed_regions = pandas.read_csv('SAVS_CA_CTGA_merged.bed', sep="\t")
print(len(bed_regions))
#Import hd5 file
SAVScapture_callset = h5py.File('SAVS_CA_downsampled_CT-GAremoved_13May21.h5', mode='r')

#create index for each contig and position
chrom = SAVScapture_callset['variants/CHROM']
pos = SAVScapture_callset['variants/POS']
idx = allel.SortedMultiIndex(chrom,pos)

#make chunked variant table
vtbl = allel.VariantChunkedTable(SAVScapture_callset['variants'], names=['CHROM', 'POS', 'REF', 'ALT', 'QUAL'])


#make genotype table
genotypes = allel.GenotypeChunkedArray(SAVScapture_callset['calldata/GT'])

'''
ac = gt.count_alleles()[:]
mult_allele = np.count_nonzero(ac.max_allele() > 1)
singleton = np.count_nonzero((ac.max_allele() == 1) & ac.is_singleton(1))
print("No. multiallelic snps", mult_allele, "No. singleton snps", singleton, sep='\t')
flt = (ac.max_allele() == 1) & (ac[:, :2].min(axis=1) > 1)
gf = gt.compress(flt, axis=0)
print(len(gf))
genotypes = gf.to_n_alt()
'''

#get array of populations/time period
populations = samples.Time_Pop.unique()
times = samples.SamplePeriod.unique()
#SampleDecades = samples.SampleDecade.unique()
stats = ["pi", "tajimaD", "wattersonTheta"]

#set up empty df that includes time period and population 
#population
#historic modern
#then underneath each column list stats
mi = pandas.MultiIndex.from_product((populations,stats), names=["Population","Statistic"])
df = pandas.DataFrame(columns=mi, index=range(len(bed_regions)))

print(df)

print(bed_regions.head())

pos = vtbl['POS']

ac = genotypes.count_alleles()

chrom = np.unique(bed_regions['CHROM'])
print(len(chrom))

pi_values = np.array([])
num_pi = 0
counter = 0
for item in chrom:
	print(counter)
	beds = bed_regions[bed_regions['CHROM'].isin([item])]
	if len(beds) > 1:
		windows = np.asarray((beds['START_POS'],beds['END_POS'])).T

	else:
		START = int(beds['START_POS'])
		END = int(beds['END_POS'])
		windows = np.array([[START,END]])
		
	#subset vtbl and genotypes based on chromosome
	loc_chrom = idx.locate_range(item)
	vtbl_chrom = vtbl[loc_chrom]
	vtbl_chrom_pos = vtbl_chrom['POS']
	
	gt_chrom = genotypes[loc_chrom]

	for pop in populations:
		#for time in times:
		#subset genotypes and do allele counts 
		#Population_time = '%s_%s' % (str(pop),str(time))
		sample_selection = samples.Time_Pop.isin({pop}).values
		samples_subset = samples[sample_selection]
		samples_subset.reset_index(drop=True, inplace=True)
		variant_selection = vtbl_chrom.eval('QUAL > 1')[:]
		#ac_pop = genotypes.subset(variant_selection, sample_selection)
		genotypes_subset = genotypes.subset(variant_selection, sample_selection)
		ac_pop = genotypes_subset.count_alleles()
		
		'''
		seg_sites = []
		for window in windows:
			loc_region = idx.locate_range(item,window[0],window[1])
			vtbl_region = vtbl[loc_region]
			variant_region_selection = vtbl_region.eval('QUAL > 1')[:]
			genotypes_window_subset = genotypes.subset(variant_selection, sample_selection)
			ac_window = genotypes_window_subset.count_alleles()
			seg_site = ac_window.count_segregating()
			seg_sites.append(seg_site)
	
		seg_sites2 = np.asarray(seg_sites)
		print(seg_sites2)
		'''
		
		#calculate summary stats
		pi, windows, n_bases, counts = allel.windowed_diversity(vtbl_chrom_pos,ac_pop,windows=windows)
		theta,_,_,_ = allel.windowed_watterson_theta(vtbl_chrom_pos,ac_pop,windows=windows)
		tajD,_,_ = allel.windowed_tajima_d(vtbl_chrom_pos,ac_pop,windows=windows)

		#append to pandas dataframe
		index_range = range(num_pi,num_pi+len(pi))
		df.loc[index_range,(pop,'pi')] = pi
		#df.loc[index_range,(pop,time,'S')] = seg_sites2
		df.loc[index_range,(pop,'wattersonTheta')] = theta
		df.loc[index_range,(pop,'tajimaD')] = tajD

		
	num_pi += len(pi)

	counter += 1
print(df)
print(len(df))
print(num_pi)

df.to_csv("CaptureRegions_SumStats_CTGAremoved_downsampled_TimeSeries_1Sept2021.csv", na_rep='NaN')

'''
#calculate some summary stats of Fst, I use weir-cockerham Fst here.
Fst = np.asarray(snp_fst_hudson)

#set all negative Fst values to 0
#and use the numpy package to calculate basic stats on Fst
Fst[Fst<0]=0
Fst_max = np.nanmax(Fst)
Fst_sd = np.nanstd(Fst)
Fst_mean = np.nanmean(Fst)
Five_x_mean = Fst_mean + 5*Fst_sd
print("Mean Fst:", Fst_mean, "Max Fst:",Fst_max,"Mean Fst + 5*sd:", Five_x_mean, sep="\t")

#plot Fst distribution with threshold in red
fig, ax = plt.subplots(figsize=(9, 7))
sns.despine(ax=ax,offset=10)
ax.hist(Fst,bins=25,edgecolor='black')
plt.axvline(x=Five_x_mean, color='red', linestyle='dashed')
ax.set_xlabel('Fst')
ax.set_ylabel('counts')
print(len(Fst))
Fig_name = "PerSNP_HudsonFst_distribution.png"
fig.savefig(Fig_name, dpi=300)	

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