#!/bin/bash
use blastp2.3.0
makeblastdb -in $1 -dbtype nucl
makeblastdb -in $2 -dbtype nucl
makeblastdb -in $3 -dbtype nucl

blastn -query $5 -db $2 -outfmt 6 -evalue 1e-1 > $4_trnas.blastn
blastn -query $5 -db $3 -outfmt 6 -evalue 1e-1 > $4_rrnas.blastn
tblastx -query $4 -db $1 -outfmt 6 -evalue 1e-1 > $4_genes.tblastx

