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
SAVScapture_callset = h5py.File('/Users/phbenham/Dropbox/CA_SAVS_temporalGenomics_analyses_ms/Bioinformatics_analyses/Data/Passerculus_California_FinalCapData_Downsampled_Filtered_31March2021_over4.5x.sorted.h5', mode='r')


#import list of samples and other data
samples_fn = '/Users/phbenham/Dropbox/CA_SAVS_temporalGenomics_analyses_ms/Bioinformatics_analyses/Data/CA_sample_data_forGenDiv2.csv'
samples = pandas.read_csv(samples_fn)
print(samples.head())

#create index for each contig and position
chrom = SAVScapture_callset['variants/CHROM']
pos = SAVScapture_callset['variants/POS']
idx = allel.SortedMultiIndex(chrom,pos)
chrom2 = np.asarray(chrom)
pos2 = np.asarray(pos)

#make chunked variant table
vtbl = allel.VariantChunkedTable(SAVScapture_callset['variants'], names=['CHROM', 'POS', 'REF', 'ALT', 'QUAL'])

#make genotype table
genotypes = allel.GenotypeChunkedArray(SAVScapture_callset['calldata/GT'])
populations = ["Orange", "Morro", "BayArea", "Humboldt"]


Outfile = open("Fst_saltmarsh_divergenceStats_28July2021b.txt", "w")

append_df = []
for pop in populations:
	print(pop, file=Outfile)
	Mod1 = 'Fresh_Modern'
	Mod2 = '%s_Modern' % pop
	Hist1 = 'Fresh_Historic'
	Hist2 = '%s_Historic' % pop

	n_mod1 = np.count_nonzero(samples.Time_Pop == Mod1)
	n_mod2 = np.count_nonzero(samples.Time_Pop == Mod2)
	n_hist1 = np.count_nonzero(samples.Time_Pop == Hist1)
	n_hist2 = np.count_nonzero(samples.Time_Pop == Hist2)
	print(Mod1,n_mod1, sep='\t')
	print(Mod2,n_mod2, sep='\t')
	print(Hist1,n_hist1, sep='\t')
	print(Hist2,n_hist2, sep='\t')

	Total = n_mod1 + n_mod2 + n_hist1 + n_hist2
	print(Total)

	# dictionary mapping population names to sample indices
	subpops = {
		Mod1: samples[samples.Time_Pop == Mod1].index,
		Mod2: samples[samples.Time_Pop == Mod2].index,
		Hist1: samples[samples.Time_Pop == Hist1].index,
		Hist2: samples[samples.Time_Pop == Hist2].index,
	}
	# allele counts
	acs = genotypes.count_alleles_subpops(subpops)

	#singletons = np.count_nonzero(acs.is_singleton(1))
	#print(singletons)
	acu = allel.AlleleCountsArray(acs[Mod1][:] + acs[Mod2][:] + acs[Hist1][:] + acs[Hist2][:])
	flt = acu.is_segregating() & (acu[:, :2].min(axis=1) > 5)
	chrom3 = chrom2.compress(flt)
	pos3 = pos2.compress(flt)

	ac_mod1 = allel.AlleleCountsArray(acs[Mod1].compress(flt, axis=0)[:, :2])
	ac_mod2 = allel.AlleleCountsArray(acs[Mod2].compress(flt, axis=0)[:, :2])
	ac_hist1 = allel.AlleleCountsArray(acs[Hist1].compress(flt, axis=0)[:, :2])
	ac_hist2 = allel.AlleleCountsArray(acs[Hist2].compress(flt, axis=0)[:, :2])
	genotype = genotypes.compress(flt, axis=0)

	loc_asc = ac_mod1.is_segregating() & ac_mod2.is_segregating() & ac_hist1.is_segregating() & ac_hist2.is_segregating()

	n_snps = np.count_nonzero(loc_asc)

	chrom4 = chrom3.compress(loc_asc)
	pos4 = pos3.compress(loc_asc)

	acs_mod1 = ac_mod1.compress(loc_asc, axis=0)
	acs_mod2 = ac_mod2.compress(loc_asc, axis=0)
	acs_hist1 = ac_hist1.compress(loc_asc, axis=0)
	acs_hist2 = ac_hist2.compress(loc_asc, axis=0)

	#print fst pop to outfile
	fst_hudson_hist, se_hudson_hist, vb_hudson_hist, _ = allel.blockwise_hudson_fst(acs_hist1, acs_hist2, blen=1000)
	print('%.04f +/- %.04f (Hudson_hist)' % (fst_hudson_hist, se_hudson_hist), file=Outfile)

	fst_hudson_mod, se_hudson_mod, vb_hudson_mod, _ = allel.blockwise_hudson_fst(acs_mod1, acs_mod2, blen=1000)
	print('%.04f +/- %.04f (Hudson_mod)' % (fst_hudson_mod, se_hudson_mod), file=Outfile)

	fst_hudson_BA, se_hudson_BA, vb_hudson_BA, _ = allel.blockwise_hudson_fst(acs_mod2, acs_hist2, blen=1000)
	print('%.04f +/- %.04f (Hudson_BA)' % (fst_hudson_BA, se_hudson_BA), file=Outfile)


	num_hist, den_hist = allel.hudson_fst(acs_hist1, acs_hist2)
	Fst_hist = num_hist / den_hist

	num_mod, den_mod = allel.hudson_fst(acs_mod1, acs_mod2)
	Fst_mod = num_mod / den_mod

	Fst_hist[Fst_hist<0]=0
	Fst_max_hist = np.nanmax(Fst_hist)
	Fst_sd_hist = np.nanstd(Fst_hist)
	Fst_mean_hist = np.nanmean(Fst_hist)
	Five_x_mean_hist = Fst_mean_hist + 5*Fst_sd_hist

	#print these data independently for each outfile to outfile
	print("Historic Mean Fst:", Fst_mean_hist, "Max Fst:",Fst_max_hist,"Mean Fst + 5*sd:", Five_x_mean_hist, sep="\t", file=Outfile)
	
	#make histogram of Fst distribution add line for outlier threshold and ESR1 
	fig, ax = plt.subplots(figsize=(9, 7))
	sns.despine(ax=ax,offset=10)
	ax.hist(Fst_hist,bins=25,edgecolor='black')
	plt.axvline(x=Five_x_mean_hist, color='red', linestyle='dashed')
	plt.axvline(x=0.45, color='gray', linestyle='dashed')
	plt.axvline(x=0.29, color='gray', linestyle='dashed')
	ax.set_xlabel('Fst')
	ax.set_ylabel('counts')

	Fig_name = "%s_PerSNP_HudsonFst_distribution_July2021.pdf" % pop
	fig.savefig(Fig_name, dpi=300)	
	
	
	Fst_mod[Fst_mod<0]=0
	Fst_max_mod = np.nanmax(Fst_mod)
	Fst_sd_mod = np.nanstd(Fst_mod)
	Fst_mean_mod = np.nanmean(Fst_mod)
	Five_x_mean_mod = Fst_mean_mod + 5*Fst_sd_mod

	print("Modern Mean Fst:", Fst_mean_mod, "Max Fst:",Fst_max_mod,"Mean Fst + 5*sd:", Five_x_mean_mod, sep="\t", file=Outfile)
	OutlierSNPs = Fst_hist[Fst_hist>Five_x_mean_hist]
	ModOutlierSNPs = Fst_mod[Fst_mod>Five_x_mean_mod]
	NumberOutlierSNPs = len(OutlierSNPs)
	NumberModSNPs = len(ModOutlierSNPs)
	print("Total SNPs:", n_snps, "Number historic outlier SNPs:",NumberOutlierSNPs, file=Outfile)
	print("Number modern outlier SNPs:",NumberModSNPs, "\n", file=Outfile)
	
	delta_fst = Fst_mod-Fst_hist

	#fst_data.pop for each independent dataframe/population yada yada
	Fst_data = pandas.DataFrame({'Chrom':chrom4, 'Start': pos4, 'Fst_hist':Fst_hist, 'Fst_mod':Fst_mod, 'deltaFst':delta_fst})
	End_pos = Fst_data['Start'] + 1
	Fst_data.insert(loc=2, column='End', value=End_pos)
	Fst_data["SNPtype"] = np.where(Fst_data['Fst_hist']>Five_x_mean_hist, "hist_outlier", "not_outlier")
	Fst_data["Population"] = "%s" % pop
	append_df.append(Fst_data)
	
	InBed = "%s_perSNP_FstDivergence_downsampled_27July2021_0.05maf.bed" % pop
	OutBed = "%s_perSNP_FstDivergence_downsampled_27July2021_annot_0.05maf.txt"  % pop
	Fst_data.to_csv(InBed, na_rep='NaN', sep='\t', index=False, header=False)
	os.system("bedtools intersect -wb -a FullData_WTSP_CaptureRegions_merged.bed -b %s > %s" % (InBed, OutBed))
	
Outfile.close()

SNP_df = pandas.concat(append_df)	

print(SNP_df.head())


FinalSNP_df = SNP_df.astype({'Fst_hist': 'float64','Fst_mod': 'float64','deltaFst': 'float64'})
print(FinalSNP_df.dtypes)

fig, ax = plt.subplots(figsize=(12, 6))
sns.despine(ax=ax, offset=5)
ax = sns.violinplot(x='Population', y='deltaFst', hue='SNPtype', data=FinalSNP_df)#, errwidth=0.25, capsize=0.1)	
ax.legend(bbox_to_anchor=(1, 1), loc='upper left', fontsize='medium')
ax.axhline(0, color="black", lw=0.5)
#ax.set_ylabel(ylab)
fig.tight_layout()
fig.savefig("deltaFst_changeDivergence_downsampled_28July2021b.pdf", dpi=300)


'''

retain_snps = []
loc_asc_snps = []
Fst_mod_all = []
Fst_hist_all = []
Fst_BA_all = []
Fst_mod_seg = []
Fst_hist_seg = []
Fst_BA_seg = []

n = [1,2,3,4,5,6,7,8,9,10]

Fst_results = pandas.DataFrame(list(zip(n,retain_snps,loc_asc_snps,Fst_mod_all,Fst_hist_all,Fst_mod_seg,Fst_hist_seg,Fst_BA_all,Fst_BA_seg)),
columns =['n','Post-maf_SNPs', 'Ascertained_SNPs','Fst_mod_all','Fst_hist_all','Fst_mod_seg','Fst_hist_seg','Fst_BA_all','Fst_BA_seg'])

print(Fst_results)

Fst_results.to_csv("MAF_impact_Fst.csv")

#print('not restricted to snps segregating across all populations')
fst_hudson_hist, se_hudson_hist, vb_hudson_hist, _ = allel.blockwise_hudson_fst(ac_hist1, ac_hist2, blen=1000)
#print('%.04f +/- %.04f (Hudson_hist)' % (fst_hudson_hist, se_hudson_hist))
Fst_hist_all.append(fst_hudson_hist)

fst_hudson_mod, se_hudson_mod, vb_hudson_mod, _ = allel.blockwise_hudson_fst(ac_mod1, ac_mod2, blen=1000)
#print('%.04f +/- %.04f (Hudson_mod)' % (fst_hudson_mod, se_hudson_mod))
Fst_mod_all.append(fst_hudson_mod)

fst_hudson_BA, se_hudson_BA, vb_hudson_BA, _ = allel.blockwise_hudson_fst(ac_mod2, ac_hist2, blen=1000)
#print('%.04f +/- %.04f (Hudson_BA)' % (fst_hudson_BA, se_hudson_BA))
Fst_BA_all.append(fst_hudson_BA)
'''

