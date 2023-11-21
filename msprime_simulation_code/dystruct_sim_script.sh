#!/usr/bin/env bash
#export OMP_NUM_THREADS=4

while getopts ":i:s:o:" options; do
	case "${options}" in
		i)
			INPUT=${OPTARG} ;;
		s)
			SAMPTIMES=${OPTARG} ;;
		o)
			OUT=${OPTARG} ;;	
	esac
done		

nSNPs=$(wc -l < "$INPUT")

#make folder for output of each 

#run k = 3
echo "./bin/dystruct --input $INPUT \
               --generation-times $SAMPTIMES \
               --output ${INPUT}_${OUT} \
               --npops 3 \
               --nloci $nSNPs \
               --pop-size 100000 \
               --seed 1145 \
               --hold-out-fraction 0.1 \
               --hold-out-seed 55307"
               
#run k = 4
echo "./bin/dystruct --input $INPUT \
               --generation-times $SAMPTIMES \
               --output ${INPUT}_${OUT} \
               --npops 4 \
               --nloci $nSNPs \
               --pop-size 100000 \
               --seed 1145 \
               --hold-out-fraction 0.1 \
               --hold-out-seed 55307"

#run k = 5               
echo "./bin/dystruct --input $INPUT \
               --generation-times $SAMPTIMES \
               --output ${INPUT}_${OUT} \
               --npops 5 \
               --nloci $nSNPs \
               --pop-size 100000 \
               --seed 1145 \
               --hold-out-fraction 0.1 \
               --hold-out-seed 55307"