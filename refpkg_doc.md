Each TIPP reference package contains the following files. This is an informal
description of what is in each one.

*CONTENTS.json*:
	Standard JSON file to define a reference package for TIPP. I just
		copied this one right from the 16S_bacteria package and changed
		the relevant values in the "files" section
		
*pasta.fasta*:
	Sequences identified by FetchMG were aligned as proteins in PASTA, then, codons from the corresponding nucleotide sequences were subbed in.

*pasta.hmm*:
	hmm trained on the whole alignment

*pasta.taxonomy*
	- This is the refined taxonomy used for placement by TIPP. This tree was fit by running RAxML with the unrefined taxonomy as a constraint tree. The tree was first fit without branch lengths, then branch lengths were found by a call to raxml to optimize model parameters.
	- The RAxML commands used to generate this tree was:
		```
		raxmlHPC-PTHREADS-AVX -m GTRCAT -F -T 20 -p 1111 -g unrefined.taxonomy.renamed -s pasta.fasta -n refined -w $work/raxml_output
		raxmlHPC-PTHREADS-AVX -m GTRGAMMA -f e -t raxml_output/RAxML_result.refined -T 20 -p 1111 -s pasta.fasta -n optimized -w $work/raxml_output
		```
*pasta.tree*:
	- Contains the unconstrained gene tree estimated on pasta.fasta (with duplicates). This tree was fit using the same raxml commands as above, but without the '-g' argument.

*pasta.taxonomy.RAxML_info*:
	- RAxML model information for the pasta.taxonomy tree. Note that this is NOT exactly the RAxML_info that was output by raxml. The current version of RAxML creates an _info file that is not compatible with pplacer, so it must be converted to the old form. We do this by using a dummy RAxML_info file in the old form and inputing the relevant data.

*pasta.tree.RAxML_info*
	- Likewise for pasta.tree. This is provided in case a user wants to use SEPP on the unconstrained Tree.

*species.mapping*
	- The comma-delimited mapping between sequence names in the tree/alignment and NCBI taxon IDs.

*taxonomy.table*
	- The NCBI taxonomy table, as TIPP expects.

	
*pasta.tree.dedup*:
	- Contains a sub-tree of 'pasta.tree' above but with duplicate sequences removed. RAxML is not perfect about having the duplicate sequences be exactly contained in a single clade of its output, but modulo branch length it effectively is. So this tree arbitrarily picks one representative sequence from each equivalence class (see below) and pulls
	out the subtree with those sequences. Name-mapping info is in the json file below.
	
*dedup_name_map.json*
	A data structure in JSON form representing the set of duplicates in 'pasta.fasta'. It defines a set of equivalence classes for sets of identical sequences, and provides some mappings to work with the dedup tree. See [here](https://github.com/MGNute/phylogeny_utilities/blob/0b05ad50f201585ae360d47d4cd91ea210ec8122/utilities.py#L1063) for a function creating this data structure.
	
	Importing this file in python using the json library will create a dictionary object with the following:
	- 'fasta_names_to_equiv_class': dictionary object mapping every sequence in the original pasta.fasta file to an equivalence-class.
		key:	orignal pasta.fasta sequence names (str)
		value:	equivalence-class-id (int)
		
	- 'fasta_equivalence_class_definitions': a dictionary mapping equivalence-class-id to relevant metadata about it. Each value is itself a dictionary:
		key:	equivalence-class-id (int)
		value:	dictionary object with the following items:
			'seq':		(str) sequence string.
			'members': 	(list) sequence names from original pasta.fasta sharing the same sequence.
			'copynum':	number of duplications (i.e. length of the 'members' list, provided for convenience).
	
	- 'deduped_name_to_equivalence_class': dictionary containing a key for every name in the tree 'pasta.tree.dedup' and a value with its equivalence-class-id.
	


