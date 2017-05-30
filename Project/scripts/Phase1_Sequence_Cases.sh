# Phase1_Sequence.py

build_sequence=/frazer01/home/joreyna/repos/CSE-280a/Project/scripts/Phase1_Sequence.py
# Currently case and numVNTR will be equal to one anther. The VNTR will
# be hard coded as well but will change when investigating sequence 
# diversity.

case=1
seqLen=10000
vntr='GCACGCTGCTGTGTAGTGGAGAAAGGGCAGGCAGCGAGCAAGCGTGTACAAGGTATATACGTGCC' 
numVNTR=1
numMuts=3

prefix=sequence_${seqLen}_case_${case}
case_dir="$(pwd)/../output/pipeline/sample/${prefix}/"
out_fn="${case_dir}/${prefix}.fa"
cmd="python $build_sequence --mutation_type group_random_mutations --rlen 150 -o $out_fn --gen_ref $seqLen $vntr $numVNTR $numMuts"
printf "Evaluating: %s\n" "$cmd"
eval $cmd

