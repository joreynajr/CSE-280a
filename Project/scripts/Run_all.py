import os
import subprocess 
# Run_All.py
proj_dir = '/frazer01/home/joreyna/repos/CSE-280a/Project/'
build_sequence = '/frazer01/home/joreyna/repos/CSE-280a/Project/scripts/Phase1_Sequence.py'
map_reads = '/frazer01/home/joreyna/repos/CSE-280a/Project/scripts/Phase2_Mapping.py'
cases = [1, 2,]#3, 4, 5, 6, 7, 8, 9, 10]
cases = [3, 4,]# 5, 6, 7, 8, 9, 10]
cases = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10]

coverages = [10, 20,] #30, 40, 50, 60, 70, 80, 90, 100]
coverages = [30, 40, 50, 60, 70, 80, 90, 100]
coverages = [10, 20, 30, 40, 50, 60, 70, 80, 90, 100]

#cases = [0:1]
#cases = [4]
seq_len = 10000
vntr='GCACGCTGCTGTGTAGTGGAGAAAGGGCAGGCAGCGAGCAAGCGTGTACAAGGTATATACGTGCC' 

for case in cases:

	num_vntr = case 
	num_muts = 3
	sample = 'sequence_{}_case_{}'.format(seq_len, case)
	case_dir = os.path.join(proj_dir, 'output/pipeline/sample/', sample)
	fq_fn = os.path.join(case_dir, sample + '.fa')

	#print('\n') 
	#print('sample: {}'.format(sample))
	#print('case_dir: {}'.format(case_dir))
	#print('fq_fn: {}'.format(fq_fn))

	#print('\n')
	cmd = "python {} --mutation_type group_random_mutations ".format(build_sequence)
	cmd += "--rlen 150 -o {} --gen_ref {} {} {} {}".format(fq_fn, seq_len, vntr, num_vntr, num_muts)
	#print("Evaluating: {}\n".format(cmd))
	subprocess.call(cmd, shell = True)

	# RUN sample at different coverages 
	for coverage in coverages:

			#print('\n')
			cmd = "python Phase2_Mapping.py "
			cmd += "--sample {} --coverage {} --copy-len {} --non-vntr-len {}".format(sample, coverage, len(vntr), seq_len)
			#print("Evaluating: {}\n".format(cmd))
			subprocess.call(cmd, shell=True)
