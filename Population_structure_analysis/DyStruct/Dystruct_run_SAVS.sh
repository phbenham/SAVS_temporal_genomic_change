#!/usr/bin/env bash
export OMP_NUM_THREADS=3

for K in $(seq 1 10); do 
	dystruct --input SAVS_CAplusoutgroups.geno \
				   --generation-times SAVS_CAplusOutgroups_times.txt \
				   --output ./dystruct_out_NE100k_run4/${K}_out \
				   --npops ${K} \
				   --nloci 65582 \
				   --pop-size 100000 \
				   --seed 21349 \
				   --hold-out-fraction 0.1 \
				   --hold-out-seed 9567
	done			   