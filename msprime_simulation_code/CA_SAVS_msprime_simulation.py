import demes
import msprime
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np

time_mig_change = [5,10,15,20,25]
mig_rate_change = [0.00015, 0.0015, 0.015, 0.15, 0.5]
demes_model = "Historic_bestMod_demes.yml"
graph = demes.load(demes_model)
demography = msprime.Demography.from_demes(graph)
outfile = open("ModelSNPs.txt", "w")
for time in time_mig_change:
	for mig in mig_rate_change:
		print("Now running model- gen:", time, "mig:", mig, sep="\t")
		#load best-fitting model from gadma as demes file to msprime.


		#at time T add migration rate change between bay area and 
		demography.add_migration_rate_change(time=0, rate=mig, source="BayArea", dest="EastCal")
		demography.add_migration_rate_change(time=0, rate=0.0, source="EastCal", dest="BayArea")
		demography.add_migration_rate_change(time=time, rate=0.0001332, source="BayArea", dest="EastCal")
		demography.add_migration_rate_change(time=time, rate=0.0001559, source="EastCal", dest="BayArea")
		demography.sort_events()

		#dd = demography.debug()
		#print(dd)

		samples = [msprime.SampleSet(num_samples=25, population="Orange", time=0, ploidy=2)] + \
				  [msprime.SampleSet(num_samples=25, population="BayArea", time=0, ploidy=2)] + \
				  [msprime.SampleSet(num_samples=25, population="EastCal", time=0, ploidy=2)] + \
				  [msprime.SampleSet(num_samples=25, population="Orange", time=50, ploidy=2)] + \
				  [msprime.SampleSet(num_samples=25, population="BayArea", time=50, ploidy=2)] + \
				  [msprime.SampleSet(num_samples=25, population="EastCal", time=50, ploidy=2)]

		rep_ts = msprime.sim_ancestry(samples=samples, demography=demography, sequence_length=650, num_replicates=60000)
		#rep_ts = msprime.sim_ancestry(samples=samples, demography=demography, sequence_length=650, num_replicates=5)
		genotypes = None
		for rep in rep_ts:
			mts = msprime.sim_mutations(rep, rate=4.6e-9,random_seed=12345)
			var = mts.genotype_matrix()
			SNP = len(var)
			if SNP > 0:
				idx = np.random.randint(SNP, size=1)
				rand_SNP = var[idx,:]
				if np.sum(rand_SNP) > 1:
					if genotypes is None:
						genotypes = rand_SNP
					else:
						genotypes = np.vstack((genotypes, rand_SNP))

		genotypes_modern_or = genotypes[:,0:25] + genotypes[:,25:50]
		genotypes_modern_ba = genotypes[:,50:75] + genotypes[:,75:100]
		genotypes_modern_ec = genotypes[:,100:125] + genotypes[:,125:150]
		genotypes_hist_or = genotypes[:,150:175] + genotypes[:,175:200]
		genotypes_hist_ba = genotypes[:,200:225] + genotypes[:,225:250]
		genotypes_hist_ec = genotypes[:,250:275] + genotypes[:,275:300]


		genotypes = np.hstack((genotypes_modern_or,
							   genotypes_modern_ba,
							   genotypes_modern_ec,
							   genotypes_hist_or,
							   genotypes_hist_ba,
							   genotypes_hist_ec))

		#sample_times = [50 for i in range(75)] + [0 for i in range(75)]
	
		Filename = "Sim_dystruct_Gen%s_Mig%s.geno" % (str(time), str(mig))
		print(Filename, len(genotypes), sep="\t", file=outfile)
		np.savetxt(Filename, genotypes, fmt="%i", delimiter='')
		#np.savetxt("Sim_dystruct_sampleTimes.txt", sample_times, fmt="%i", delimiter="\n")
outfile.close()
