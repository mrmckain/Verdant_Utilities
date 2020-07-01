#!/bin/bash
perl annoBTD/annotate_plastome.pl $1 $3
COUNT=0;
for file in $(echo $2 | tr "," " "); 
do
    COUNT=$(($COUNT + 1))
    perl annoBTD/get_annotated_regions_fromverdant.pl ./sequenceFiles/$file ./files/$file $3 $COUNT
done
./annoBTD/run_multiblast.sh $3_annotated_regions_fromverdant_genes.fsa $3_annotated_regions_fromverdant_trnas.fsa $3_annotated_regions_fromverdant_rrnas.fsa $3_test_seqs_orffinder.fsa $1
perl annoBTD/match_orfs_to_blast_NEW_MULTIPLE_ATTEMPT2.pl $3_annotated_regions_fromverdant_genes.fsa $3_test_seqs_orffinder.fsa $3_test_seqs_orffinder.fsa_genes.tblastx $3_test_orffinder.txt $1 $3_annotated_regions_fromverdant_trnas.fsa $3_annotated_regions_fromverdant_rrnas.fsa $3_test_seqs_orffinder.fsa_trnas.blastn $3_test_seqs_orffinder.fsa_rrnas.blastn $3


