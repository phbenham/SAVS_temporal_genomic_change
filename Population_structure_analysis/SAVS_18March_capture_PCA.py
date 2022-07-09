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


########################################################################################
def fig_pca(samples,model,sub_name):
	fig,ax = plt.subplots(figsize=(7, 5))
	sns.scatterplot(x='PC1',y='PC2',hue='Population', style='Epoch', size='Epoch', edgecolor='black',markers=['o','^'],sizes=[50,125], alpha=1, data=samples)
	lgd = ax.legend(bbox_to_anchor=(1, 1), loc='upper left', fontsize='medium')
	ax.set_xlabel('PC1 (%.1f%%)' % (model.explained_variance_ratio_[0]*100))
	ax.set_ylabel('PC2 (%.1f%%)' % (model.explained_variance_ratio_[1]*100))
	OutputFileName = "SAVScapture_%s_subset_pca.png" % sub_name
	fig.savefig(OutputFileName, dpi=300, bbox_extra_artists=(lgd,), bbox_inches='tight')

########################################################################################
########################################################################################
def ld_prune(gn, size, step, threshold=.1, n_iter=5):
    for i in range(n_iter):
        loc_unlinked = allel.locate_unlinked(gn, size=size, step=step, threshold=threshold)
        n = np.count_nonzero(loc_unlinked)
        n_remove = gn.shape[0] - n
        print('iteration', i+1, 'retaining', n, 'removing', n_remove, 'variants')
        gn = gn.compress(loc_unlinked, axis=0)
    return gn            
########################################################################################
def gt_filtering(gt):
	#gt = genotype.take(samples, axis=1)
	ac = gt.count_alleles()[:]
	mult_allele = np.count_nonzero(ac.max_allele() > 1)
	singleton = np.count_nonzero((ac.max_allele() == 1) & ac.is_singleton(1))
	print("No. multiallelic snps", mult_allele, "No. singleton snps", singleton, sep='\t')
	flt = (ac.max_allele() == 1) & (ac[:, :2].min(axis=1) > 1)
	gf = gt.compress(flt, axis=0)
	print(len(gf))
	gn = gf.to_n_alt()
	gnu = ld_prune(gn, size=2500, step=200, threshold=0.25, n_iter=5)
	return gnu
########################################################################################
########################################################################################



VCF_path = "/Users/phbenham/Dropbox/CA_SAVS_temporalGenomics_analyses_ms/Bioinformatics_analyses/Data/SAVS_CA_Exons_CTandGA_removed.recode.vcf"
allel.vcf_to_hdf5(VCF_path, '/Users/phbenham/Dropbox/CA_SAVS_temporalGenomics_analyses_ms/Bioinformatics_analyses/Data/SAVS_CA_Exons_CTandGA_removed.h5', fields='*', overwrite=True)


#load hdf5 file with genotype calls and sample data
Capture_callset = h5py.File('/Users/phbenham/Dropbox/CA_SAVS_temporalGenomics_analyses_ms/Bioinformatics_analyses/Data/Passerculus_California_FinalCapData_Filtered_31March2021_over4.5x.h5', mode='r')

samples_fn = '/Users/phbenham/Dropbox/CA_SAVS_temporalGenomics_analyses_ms/Bioinformatics_analyses/Data/CA_sample_data.csv'
samples = pandas.read_csv(samples_fn)

'''
samples['State_Clade'] = samples['State'] + '_' + samples['Clade']
print(samples.head())

samples_select = samples.State.isin({'California'}).values
samples_Selected = samples[samples_select]
samples_Selected.reset_index(drop=True, inplace=True)
'''

samples_list = list(Capture_callset['samples'])
print(samples.Subspecies.value_counts())

variants = allel.VariantChunkedTable(Capture_callset['variants'], names=['POS', 'REF', 'ALT', 'DP', 'MQ'])

idx = np.arange(0,len(variants),1)
print(len(idx))
samples_callset_index = [samples_list.index(s) for s in samples['SampleID']]
samples['callset_index'] = samples_callset_index

#create a genotype array from callset. 2nd create allele counts file from genotype array
gt = allel.GenotypeChunkedArray(Capture_callset['calldata/GT'])


#definition that outputs final genotypes that have singletons/non-seg sites removed and has been ld-pruned
gn1 = gt_filtering(gt)
coords1, model1 = allel.pca(gn1, n_components=10)
samples['PC1'] = coords1[:,0]
samples['PC2'] = coords1[:,1]
samples['PC3'] = coords1[:,2]
samples['PC4'] = coords1[:,3]

#samples_Selected.to_csv('CA_Downsample_PCdata.csv')
colors = {'alaudinus':(0.067,0.467,0.2), 'beldingi':(0.667,0.267,0.6), 'brooksi':(0.2,0.133,0.533), 'nevadensis':(0.533,0.8,0.933)}

fig,(ax1,ax2) = plt.subplots(1,2,figsize=(12, 5))
#fig,ax = plt.subplots(figsize=(12, 5))

sns.scatterplot(ax=ax1, x='PC1',y='PC2',hue='Subspecies',style='SamplePeriod',palette = colors, markers=['^','o'], edgecolor='black',s=100, alpha=1, legend=False, data=samples)
#ax1.legend(bbox_to_anchor=(1, 1), loc='upper left', fontsize='medium')
ax1.set_xlabel('PC1 (%.1f%%)' % (model1.explained_variance_ratio_[0]*100))
ax1.set_ylabel('PC2 (%.1f%%)' % (model1.explained_variance_ratio_[1]*100))

sns.scatterplot(ax=ax2, x='PC3',y='PC4',hue='Subspecies', style='SamplePeriod',palette = colors, markers=['^','o'], edgecolor='black', s=100, alpha=1, data=samples)
ax2.legend(bbox_to_anchor=(1, 1), loc='upper left', fontsize='medium')
ax2.set_xlabel('PC3 (%.1f%%)' % (model1.explained_variance_ratio_[2]*100))
ax2.set_ylabel('PC4 (%.1f%%)' % (model1.explained_variance_ratio_[3]*100))


#plt.show()

fig.savefig("Passerculus_California_CT-GAremoved_nongenic_PCA_25Aug2021.png", dpi=300, bbox_inches='tight')
