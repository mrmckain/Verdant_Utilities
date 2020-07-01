#!/bin/bash

perl annotate_plastome.pl $1 $3

perl get_annotated_regions_fromverdant.pl Silene_vulgaris.fsa Silene_vulgaris $3 1

perl get_annotated_regions_fromverdant.pl Schizachyrium.fsa Schizachyrium $3 2

./run_multiblast.sh $3_annotated_regions_fromverdant_genes.fsa $3_annotated_regions_fromverdant_trnas.fsa $3_annotated_regions_fromverdant_rrnas.fsa $3_orffinder_seqs.fsa $1

perl identify_best_ref_for_orf_withscoring-3.pl $3_orffinder_seqs.fsa_trnas.blastn $3_annotated_regions_fromverdant_trnas.fsa Schizachyrium.fsa tRNA

perl identify_best_ref_for_orf_withscoring-3.pl $3_orffinder_seqs.fsa_rrnas.blastn $3_annotated_regions_fromverdant_rrnas.fsa Schizachyrium.fsa rRNA

perl identify_best_ref_for_orf_withscoring-3.pl $3_orffinder_seqs.fsa_genes.tblastx $3_annotated_regions_fromverdant_genes.fsa $3_orffinder_seqs.fsa protein

perl match_orfs_to_blast_v2.2.pl $3_annotated_regions_fromverdant_genes.fsa $3_orffinder_seqs.fsa best_orfs_for_refs_SCORE.txt $3_orffinder_coordinates.txt Schizachyrium.fsa $3_annotated_regions_fromverdant_trnas.fsa $3_annotated_regions_fromverdant_rrnas.fsa $3_orffinder_seqs.fsa_trnas.blastn $3_orffinder_seqs.fsa_rrnas.blastn $3 best_tRNA_ref_SCORE.txt best_rRNA_ref_SCORE.txt
