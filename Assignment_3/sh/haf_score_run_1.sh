#!/bin/bash

#$ -pe smp 4
#$ -l short
#$ -N HAF_Scores_1
#$ -l h_vmem=4.0G
#$ -S /bin/bash
#$ -o /frazer01/home/joreyna/repos/CSE-280a/Assignment_3/err/HAF_Scores_1.out
#$ -e /frazer01/home/joreyna/repos/CSE-280a/Assignment_3/err/HAF_Scores_1.err

source activate hla

date 1>&2
cmd="python /frazer01/home/joreyna/repos/CSE-280a/Assignment_3/Generate_HAF_Scores.py 1 0 0.1"
echo Executing: $cmd 1>&2
eval $cmd
date 1>&2
