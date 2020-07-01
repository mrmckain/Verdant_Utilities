#!/bin/bash

for k in Fixed_Plastomes/*fsa
do
	NAMESP=$(echo $k | cut -f2 -d "/")
        echo $NAMESP
        SHORT=$(echo $NAMESP | cut -f1 -d ".")
        echo $SHORT	
	perl /mrm/Hosta_Phylogenomics/Annotations/annoBTD_multiref/verdant2genbank.pl $k VERDANT_files/$SHORT\_VERDANT_cleaned_annotation.txt $SHORT /mrm/Hosta_Phylogenomics/Annotations/annoBTD_multiref/Gene_products_plastomes.txt
done
