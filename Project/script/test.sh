# Deprecated 
#if [[ "$1" == "1" ]];
#	then
#		python Generating_VNTR_sequences.py --help	
#fi

seqLen=20
vntr=AAAAA 
numVNTR=3
numMuts=2

# Case 1

#printf "\n-------------------------------- Case 1 --------------------------------"
## Testing the outer padding, generating a sequence, mutating the VNTR and making copies.
## Using the default mutation type. 
#message="\nTesting a sequence of length $seqLen. The vntr is $vntr and will be copied"
#message="$message $numVNTR times with $numMuts mutations."
#cmd="python Generating_VNTR_sequences.py --outer_pad $seqLen $vntr $numVNTR $numMuts"
#eval $cmd

# Case 2

#printf "\n-------------------------------- Case 2 --------------------------------"
## Testing the outer padding, inner padding, generating a sequence, mutating the VNTR 
## and making copies. Using the default mutation type individual random mutations 
#message="Testing a sequence of length $seqLen. The vntr is $vntr and will be copied"
#message="$message $numVNTR times with $numMuts mutations. Adding a pad between copies"
#message="$message VNTR's"
#cmd="python Generating_VNTR_sequences.py --outer_pad --inner_pad $seqLen $vntr $numVNTR $numMuts"
#eval $cmd

# Case 3

#printf "\n-------------------------------- Case 3 --------------------------------"
## Testing the outer padding, inner padding, generating a sequence, mutating the VNTR 
## and making copies. Using the group_random_mutations mutation type. 
#message="Testing a sequence of length $seqLen. The vntr is $vntr and will be copied"
#message="$message $numVNTR times with $numMuts mutations. Adding a pad between copies"
#message="$message VNTR's"
#cmd="python Generating_VNTR_sequences.py --outer_pad --inner_pad  --mutation_type group_random_mutations $seqLen $vntr $numVNTR $numMuts"
#eval $cmd


# Case 4

seqLen=200
vntr=AAAAA 
numVNTR=3
numMuts=2

printf "\n-------------------------------- Case 4 --------------------------------"
cmd="python Generating_VNTR_sequences.py --outer_pad --mutation_type group_random_mutations $seqLen $vntr $numVNTR $numMuts"
eval $cmd
