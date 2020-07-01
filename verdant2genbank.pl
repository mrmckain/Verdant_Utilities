#!/usr/bin/perl -w

use strict;
#USAGE: perl verdant2genbank.pl sequencefile products_file annotation_file
my %products;
open my $file, "<", $ARGV[1];
while(<$file>){
	chomp;
	my @tarray = split/\t/;
	$products{$tarray[0]}=$tarray[1];
}

my %annotation;
my $genus;
my $species;
my $accession;

open my $seqfile, "<", $ARGV[0];
while(<$seqfile>){
	chomp;
	if(/>/){
		my @tarray = split /_/;
		$genus = substr($tarray[0],1);
		$species = $tarray[1];
		$accession = $tarray[6];
	}
}


open my $out_anno, ">", $ARGV[0] . ".tbl";
print $out_anno ">Features $genus\_$species\_$accession \"\[organism=$genus $species\]\"\n";
open my $anno_file, "<",   $ARGV[2];
while(<$anno_file>){
	chomp;
	if(/\~/){
		next;
	}
	if(/intron/){
		next;
	}
	if(/LSC/ || /IRA/ || /IRB/ || /SSC/ || /FULL/){
		next;
	}
	if(/rps12/){
		if(/rps12_5end/){
			my @tarray=split/\s+/;
			
			$annotation{rps12}{$tarray[3]}{"5end"}{Start}=$tarray[1];
			$annotation{rps12}{$tarray[3]}{"5end"}{End}=$tarray[2];	
		}
		if(/rps12_3end/){
			my @tarray=split/\s+/;
			
			$tarray[0] =~ /_(exon\d)/;
			$annotation{rps12}{$tarray[3]}{$1}{Start}=$tarray[1];
			$annotation{rps12}{$tarray[3]}{$1}{End}=$tarray[2];
		}
	}
	elsif(/exon/){
		  my @tarray=split/\s+/;
		  
			$tarray[0] =~ /(.*?)_(exon\d)/;
			$annotation{$1}{$tarray[3]}{$2}{Start}=$tarray[1];
			$annotation{$1}{$tarray[3]}{$2}{End}=$tarray[2];
	}
	else{
		my @tarray=split/\s+/;
			
			$tarray[0] =~ /_(exon\d)/;
			$annotation{$tarray[0]}{$tarray[3]}{Full}{Start}=$tarray[1];
			$annotation{$tarray[0]}{$tarray[3]}{Full}{End}=$tarray[2];
	}
}
my %final_anno;
my %gene_pieces;
my %gene_dir;
for my $geneid (keys %annotation){
	if($geneid eq "rps12"){
		my $neg="-";
		my $pos="+";
		my $five = "5end";
		print $out_anno "$annotation{rps12}{$neg}{$five}{End}\t$annotation{rps12}{$neg}{$five}{Start}\tgene\n";
		my %rps12_sets;
		$rps12_sets{$neg}{Min}=10000000000;
		$rps12_sets{$neg}{Max}=0;
		$rps12_sets{$pos}{Min}=10000000000;
		$rps12_sets{$pos}{Max}=0;
		for my $dir (keys %{$annotation{rps12}}){
			
			for my $exon (keys %{$annotation{rps12}{$dir}}){
				if($exon eq $five){
					next;
				}
				if($annotation{rps12}{$dir}{$exon}{Start} > $rps12_sets{$dir}{Max}){
					$rps12_sets{$dir}{Max}=$annotation{rps12}{$dir}{$exon}{Start};
				}
				if($annotation{rps12}{$dir}{$exon}{End} > $rps12_sets{$dir}{Max}){
					$rps12_sets{$dir}{Max}=$annotation{rps12}{$dir}{$exon}{End};
				}
				if($annotation{rps12}{$dir}{$exon}{Start} < $rps12_sets{$dir}{Min}){
					$rps12_sets{$dir}{Min}=$annotation{rps12}{$dir}{$exon}{Start};
				}
				if($annotation{rps12}{$dir}{$exon}{End} < $rps12_sets{$dir}{Min}){
					$rps12_sets{$dir}{Min}=$annotation{rps12}{$dir}{$exon}{End};
				}
				$rps12_sets{$dir}{$annotation{rps12}{$dir}{$exon}{Start}}=1;
				$rps12_sets{$dir}{$annotation{rps12}{$dir}{$exon}{End}}=1;
			}
		}
		 $neg="-";
		$pos="+";
		print $out_anno "$rps12_sets{$neg}{Max}\t$rps12_sets{$neg}{Min}\n";
		print $out_anno "\t\t\tgene\trps12\n";
		print $out_anno "\t\t\texception\ttrans-splicing\n";
		print $out_anno "$annotation{rps12}{$neg}{$five}{End}\t$annotation{rps12}{$neg}{$five}{Start}\tCDS\n";
		print $out_anno "$annotation{rps12}{$neg}{exon1}{End}\t$annotation{rps12}{$neg}{exon1}{Start}\n";
		print $out_anno "$annotation{rps12}{$neg}{exon2}{End}\t$annotation{rps12}{$neg}{exon2}{Start}\n";
		print $out_anno "\t\t\tproduct\t$products{rps12}\n";
		print $out_anno "$annotation{rps12}{$neg}{$five}{End}\t$annotation{rps12}{$neg}{$five}{Start}\texon\n";
		print $out_anno "\t\t\tnumber 1\n";
		print $out_anno "$annotation{rps12}{$neg}{exon1}{End}\t$annotation{rps12}{$neg}{exon1}{Start}\texon\n";
		print $out_anno "\t\t\tnumber 3\n";
		print $out_anno "$annotation{rps12}{$neg}{exon2}{End}\t$annotation{rps12}{$neg}{exon2}{Start}\texon\n";
		print $out_anno "\t\t\tnumber 2\n";
		print $out_anno "$annotation{rps12}{$neg}{$five}{End}\t$annotation{rps12}{$neg}{$five}{Start}\tgene\n";
		print $out_anno "$rps12_sets{$pos}{Min}\t$rps12_sets{$pos}{Max}\n";
		print $out_anno "\t\t\tgene\trps12\n";
		print $out_anno "$annotation{rps12}{$neg}{$five}{End}\t$annotation{rps12}{$neg}{$five}{Start}\tCDS\n";
		print $out_anno "$annotation{rps12}{$pos}{exon1}{Start}\t$annotation{rps12}{$pos}{exon1}{End}\n";
		print $out_anno "$annotation{rps12}{$pos}{exon2}{Start}\t$annotation{rps12}{$pos}{exon2}{End}\n";
		print $out_anno "\t\t\tproduct\t$products{rps12}\n";
		print $out_anno "$annotation{rps12}{$neg}{$five}{End}\t$annotation{rps12}{$neg}{$five}{Start}\texon\n";
		print $out_anno "\t\t\tnumber 1\n";
		print $out_anno "$annotation{rps12}{$pos}{exon1}{Start}\t$annotation{rps12}{$pos}{exon1}{End}\texon\n";
		print $out_anno "\t\t\tnumber 2\n";
		print $out_anno "$annotation{rps12}{$pos}{exon2}{Start}\t$annotation{rps12}{$pos}{exon2}{End}\texon\n";
		print $out_anno "\t\t\tnumber 3\n";

	}
	else{
		
		for my $dir (keys %{$annotation{$geneid}}){
			for my $exon (keys %{$annotation{$geneid}{$dir}}){
				if($dir eq "+"){
					my $exons = 1;			
					if($exon =~ /exon/){
						$exons = scalar keys %{$annotation{$geneid}{$dir}};
					}				
					$final_anno{$annotation{$geneid}{$dir}{$exon}{Start}}{$annotation{$geneid}{$dir}{$exon}{End}}{$geneid}=$exons;
					$gene_pieces{$geneid}=$exons;
					$gene_dir{$geneid}{$annotation{$geneid}{$dir}{$exon}{Start}}=$dir;
				}
				elsif($dir eq "-"){
					my $exons = 1;			
					if($exon =~ /exon/){
						$exons = scalar keys %{$annotation{$geneid}{$dir}};
					}				
					$final_anno{$annotation{$geneid}{$dir}{$exon}{End}}{$annotation{$geneid}{$dir}{$exon}{Start}}{$geneid}=$exons;
					$gene_pieces{$geneid}=$exons;
					$gene_dir{$geneid}{$annotation{$geneid}{$dir}{$exon}{End}}=$dir;
				}
			}
		}
	}
}
my %current_genes;
my $current_gene;
my $current_dir;
my $trnk_id;
for my $start (sort {$a<=>$b} keys %final_anno){
	for my $end (keys %{$final_anno{$start}}){
		for my $gene (keys %{$final_anno{$start}{$end}}){
			if ($gene =~ /trnK/){
				if($current_gene){
					$current_genes{$start}{$end}=1;
					my $max = 0;
					my $min =100000000;
					for my $tstart (keys %current_genes){
							for my $tend (keys %{$current_genes{$tstart}}){
								if($tstart < $min){
									$min = $tstart;
								}
								if($tend < $min){
									$min = $tend;
								}
								if($tstart > $max){
									$max = $tstart;
								}
								if($tend > $max){
									$max = $tend;
								}
						}
					}
					my $counter=1;
					print "$gene\n";
					
					$trnk_id .="$max\t$min\tgene\n\t\t\tgene\t$current_gene\n";
					for my $tstart (sort {$b<=>$a} keys %current_genes){
						for my $tend (keys %{$current_genes{$tstart}}){
							if($counter ==1 ){
								$trnk_id.= "$tstart\t$tend\ttRNA\n";
									$counter--;
								}
								else{
									$gene =~ /trn\w+-(\w+)/;
									my $anti=$1;
									$trnk_id .= "$tstart\t$tend\n";
									$trnk_id .= "\t\t\tnote\tanticodon\:$anti\n\t\t\tproduct\t$products{$current_gene}\n";
								}
							}
						}
					}
					else{
						$current_genes{$start}{$end}=1;
						$current_dir=$gene_dir{$gene}{$start};
						$current_gene=$gene;
					}
				}
			}
		}
	}	
			
$current_gene=();
$current_dir=();
%current_genes=();		
					
for my $start (sort {$a<=>$b} keys %final_anno){
	for my $end (keys %{$final_anno{$start}}){
		for my $gene (keys %{$final_anno{$start}{$end}}){
			if ($gene =~ /trnK/){
				if($trnk_id){
					print $out_anno "$trnk_id";
					$trnk_id=();
				}
				next;
			}
			if($current_gene){
				if($current_gene eq $gene){
					$current_genes{$start}{$end}=1;
					$current_dir=$gene_dir{$gene}{$start};
				}
				else{
					my $max = 0;
					my $min =100000000;
					for my $tstart (keys %current_genes){
							for my $tend (keys %{$current_genes{$tstart}}){
								if($tstart < $min){
									$min = $tstart;
								}
								if($tend < $min){
									$min = $tend;
								}
								if($tstart > $max){
									$max = $tstart;
								}
								if($tend > $max){
									$max = $tend;
								}
							}
						}
					my $counter=1;
					print "$gene\n";
					my $gene_type ="CDS";
					if($current_gene =~ /trn/){
							$gene_type="tRNA";
					}
					if($current_dir eq "-"){
						print $out_anno "$max\t$min\tgene\n\t\t\tgene\t$current_gene\n";
						for my $tstart (sort {$b<=>$a} keys %current_genes){
							for my $tend (keys %{$current_genes{$tstart}}){
								if($counter ==1 ){
									print $out_anno "$tstart\t$tend\t$gene_type\n";
									$counter--;
								}
								else{
									print $out_anno "$tstart\t$tend\n";
								}
							}
						}
						print $out_anno "\t\t\tproduct\t$products{$current_gene}\n";
						my $exon_counter=1;
						for my $tstart (sort {$b<=>$a} keys %current_genes){
							for my $tend (keys %{$current_genes{$tstart}}){
					
									print $out_anno "$tstart\t$tend\texon\n";
									print $out_anno "\t\t\tnumber\t$exon_counter\n";
									$exon_counter++;
								}
							}
							$current_gene=();
					$current_dir=();
					%current_genes=();
					}
					elsif($current_dir eq "+"){
						print $out_anno "$min\t$max\tgene\n\t\t\tgene\t$current_gene\n";
						for my $tstart (sort {$a<=>$b} keys %current_genes){
							for my $tend (keys %{$current_genes{$tstart}}){
								if($counter ==1 ){
									print $out_anno "$tstart\t$tend\t$gene_type\n";
									$counter--;
								}
								else{
									print $out_anno "$tstart\t$tend\n";
								}
							}
						}
						print $out_anno "\t\t\tproduct\t$products{$current_gene}\n";
						my $exon_counter=1;
						for my $tstart (sort {$a<=>$b} keys %current_genes){
							for my $tend (keys %{$current_genes{$tstart}}){
					
									print $out_anno "$tstart\t$tend\texon\n";
									print $out_anno "\t\t\tnumber\t$exon_counter\n";
									$exon_counter++;
								}
							}
							
					}
					$current_gene=();
					$current_dir=();
					%current_genes=();
				}

			}		
			print "$gene\n";
			if($gene_pieces{$gene} == 1){
				
				if($gene =~ /trn/){
					$gene =~ /trn\w+-(\w+)/;
					my $anti=$1;
					print $out_anno "$start\t$end\tgene\n\t\t\tgene\t$gene\n";
					print $out_anno "$start\t$end\ttRNA\n\t\t\tnote\tanticodon\:$anti\n\t\t\tproduct\t$products{$gene}\n";
				}
				elsif($gene =~ /rrn/){
					print $out_anno "$start\t$end\tgene\n$start\t$end\trRNA\n\t\t\tproduct\t$products{$gene}\n";
				}
				else{
					print $out_anno "$start\t$end\tgene\n\t\t\tgene\t$gene\n";
					print $out_anno "$start\t$end\tCDS\n\t\t\tproduct\t$products{$gene}\n";
				}
				
			}
			elsif($gene_pieces{$gene} > 1){
				$current_genes{$start}{$end}=1;
				$current_dir=$gene_dir{$gene}{$start};
				$current_gene=$gene;
			}
		
		}
	}
}
					
			
	
