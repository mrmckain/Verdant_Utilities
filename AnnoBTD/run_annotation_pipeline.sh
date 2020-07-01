#!/bin/bash
for k in $1/*
do
	NAMESP=$(echo $k | cut -f2 -d "/")
	echo $NAMESP
	SHORT=$(echo $NAMESP | cut -f1 -d ".")
	echo $SHORT
	perl single_line.pl $k
	mv $k\_single.fa $NAMESP
	perl annoBTD_multiref/annotate_plastome.pl $NAMESP $SHORT
	COUNT=0;
	for file in $2/sequences/*
	do
    		COUNT=$(($COUNT + 1))
    		echo $file >> $SHORT\_references.txt
		REF=$(echo $file |cut -f3 -d "/")
    		perl annoBTD_multiref/get_annotated_regions_fromverdant.pl $file $2/annotations/$REF $SHORT $COUNT
	done

	annoBTD_multiref/run_multiblast.sh $SHORT\_annotated_regions_fromverdant_genes.fsa $SHORT\_annotated_regions_fromverdant_trnas.fsa $SHORT\_annotated_regions_fromverdant_rrnas.fsa $SHORT\_orffinder_seqs.fsa $NAMESP
	perl annoBTD_multiref/identify_best_ref_for_orf_withscoring_v0.7.pl $SHORT\_orffinder_seqs.fsa_trnas.blastn $SHORT\_annotated_regions_fromverdant_trnas.fsa $NAMESP tRNA $SHORT $SHORT\_orffinder_coordinates.txt

	perl annoBTD_multiref/identify_best_ref_for_orf_withscoring_v0.7.pl $SHORT\_orffinder_seqs.fsa_rrnas.blastn $SHORT\_annotated_regions_fromverdant_rrnas.fsa $NAMESP rRNA $SHORT $SHORT\_orffinder_coordinates.txt

	perl annoBTD_multiref/identify_best_ref_for_orf_withscoring_v0.7.pl $SHORT\_orffinder_seqs.fsa_genes.tblastx $SHORT\_annotated_regions_fromverdant_genes.fsa $SHORT\_orffinder_seqs.fsa protein $SHORT $SHORT\_orffinder_coordinates.txt

	perl annoBTD_multiref/match_orfs_to_blast_v2.6.pl $SHORT\_annotated_regions_fromverdant_genes.fsa $SHORT\_orffinder_seqs.fsa $SHORT\_best_orfs_for_refs_SCORE.txt $SHORT\_orffinder_coordinates.txt $NAMESP $SHORT\_annotated_regions_fromverdant_trnas.fsa $SHORT\_annotated_regions_fromverdant_rrnas.fsa $SHORT\_orffinder_seqs.fsa_trnas.blastn $SHORT\_orffinder_seqs.fsa_rrnas.blastn $SHORT $SHORT\_best_tRNA_ref_SCORE.txt $SHORT\_best_rRNA_ref_SCORE.txt
done
