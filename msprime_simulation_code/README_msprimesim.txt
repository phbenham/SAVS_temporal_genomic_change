#Readme msprime simulation code.
this folder includes code for performing all the msprime simulations included in the manuscript:
Spatial variation in population genomic responses to over a century of anthropogenic change within a tidal marsh songbird 

(1) Historic_bestMod_demes.yml
	Historic best fit demographic model from GADMA2. Used as input to parameterize coalescent
	simulations in msprime.
	
(2) CA_SAVS_msprime_simulation.py
	Runs msprime simulations with different levels of migration rate change at different
	times in the past. 
	Generates genotype files to be input into DyStruct
	
(3) dystruct_sim_script.sh
	Runs dystruct on simulated datasets. Sim_dystruct_sampleTimes.txt specifies generation times
	for samples as part of input.

(4) Msprime_simulationsParser.py
	Processes output from dystruct runs on the simulated datasets to plot estimated ancestry proportions
	of eastern california in the bay area for different simulated parameter combinations (Fig. 3c).
	
(5) Msprime_vcfout.py
	msprime simulations based on historic demographic model to produce vcffiles with either
	no change through time in demographic history (NoChange_thinned.recode.vcf)
	or a shift in migration rate (Msprime_mig_thinned.recode.vcf). Estimates delta Fst and hist fst
	for plotting against observed values. 