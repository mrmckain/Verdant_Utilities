#!/usr/bin/perl -w
use strict;

#Plastome annotation 
#USAGE: 1: Full sequence 2: species

my $plastome;
open my $cp_file, "<", $ARGV[0];
while(<$cp_file>){
	chomp;
	s/\r//;
	if(/>/){
		next;
	}
	else{
		s/\s//g;
		$plastome .= $_;
	}
}
$plastome = uc $plastome;
#Identify ORFs in sequence in both directions
my %forward_orfs = &orf_finder($plastome);
my $rc_plastome = reverse($plastome);
$rc_plastome =~ tr/ATCGatcg/TAGCtagc/;
my %reverse_orfs = &orf_finder($rc_plastome);


my $pl_len = length($plastome);

#Combine ORFs into single forward-based coordinates
for my $orfid (keys %reverse_orfs){
	$reverse_orfs{$orfid} =~ /(\d+)-(\d+)/;
	my $rstart = $1;
	my $rend = $2;
	my $nstart = $pl_len - $rend-1;
	my $nend = $pl_len - $rstart-1;
	$orfid = $orfid . "-r";
	$forward_orfs{$orfid}= $nstart . "-" . $nend;
}
open my $testout, ">", "find_thebad.txt";
for my $testval (sort {$a cmp $b} keys %forward_orfs){
		print $testout "$testval\t$forward_orfs{$testval}\n";

}
#Filter ORFs to remove nested ORFs

my %position_orfs;
my %orf_size;
for my $orfid (keys %forward_orfs){
	my @tarray = split("-", $forward_orfs{$orfid});
	$orf_size{$orfid}=abs($tarray[1]-$tarray[0]);
	$position_orfs{$tarray[0]}{$tarray[1]}=$orfid;
}

my %orf_layout_genome;
for (my $j=0; $j<=length($plastome)-1; $j++){
	$orf_layout_genome{$j}=0;
}

my %good_orfs;
my @sorted_size_orfs = sort { $orf_size{$b} <=> $orf_size{$a} } keys %orf_size;
for my $sorted_orfs (@sorted_size_orfs){
	my @tarray = split("-", $forward_orfs{$sorted_orfs});
	my $covered=0;
	for (my $k=$tarray[0]; $k<=$tarray[1]; $k++){
		if($orf_layout_genome{$k} == 1){
			$covered++;
		}
	}
	if($covered/($tarray[1]-$tarray[0]) >= 2){
		next;
	}
	else{
		for (my $k=$tarray[0]; $k<=$tarray[1]; $k++){
			$orf_layout_genome{$k} = 1;
		}
		$good_orfs{$sorted_orfs}=$forward_orfs{$sorted_orfs};
	}
}


open my $outfile, ">", $ARGV[1] . "_orffinder_coordinates.txt";
open my $seqfile, ">", $ARGV[1] . "_orffinder_seqs.fsa";
for my $ssiv (sort keys %good_orfs){
		$good_orfs{$ssiv} =~ /(\d+)-(\d+)/;
		my $start = $1;
		my $end = $2;
		my $tseq = substr($plastome, $start, ($end-$start+1));
		print $outfile "$ssiv\t$forward_orfs{$ssiv}\n";
		print $seqfile ">$ssiv\n$tseq\n";
}


sub orf_finder
{
	my $cur_plastome = $_[0];

	my %codons=('TCA'=>'S','TCC'=>'S','TCG'=>'S','TCT'=>'S','TTC'=>'F','TTT'=>'F','TTA'=>'L','TTG'=>'L','TAC'=>'Y','TAT'=>'Y','TAA'=>'_','TAG'=>'_','TGC'=>'C','TGT'=>'C','TGA'=>'_','TGG'=>'W','CTA'=>'L','CTC'=>'L','CTG'=>'L','CTT'=>'L','CCA'=>'P','CCC'=>'P','CCG'=>'P','CCT'=>'P','CAC'=>'H','CAT'=>'H','CAA'=>'Q','CAG'=>'Q','CGA'=>'R','CGC'=>'R','CGG'=>'R','CGT'=>'R','ATA'=>'I','ATC'=>'I','ATT'=>'I','ATG'=>'M','ACA'=>'T','ACC'=>'T','ACG'=>'T','ACT'=>'T','AAC'=>'N','AAT'=>'N','AAA'=>'K','AAG'=>'K','AGC'=>'S','AGT'=>'S','AGA'=>'R','AGG'=>'R','GTA'=>'V','GTC'=>'V','GTG'=>'V','GTT'=>'V','GCA'=>'A','GCC'=>'A','GCG'=>'A','GCT'=>'A','GAC'=>'D','GAT'=>'D','GAA'=>'E','GAG'=>'E','GGA'=>'G','GGC'=>'G','GGG'=>'G','GGT'=>'G', 'GCN'=>'A', 'CGN'=>'R', 'GGN'=>'G', 'CCN'=>'P', 'TCN'=>'S', 'ACN'=>'T', 'GTN'=>'V');

	my %orfs;
	my $current_orf;
	my $cur_start;
	my $cur_end;
	my $orf_count=0;
	for (my $i = 0; $i<=2; $i++){
		if($current_orf){
			
			if(length($current_orf)> 18){
				$orf_count++;
				my $torf = "orf" . $orf_count;
				$orfs{$torf}=$cur_start . "-". $cur_end;
			
			}
			$current_orf=();
			$cur_start=();
			$cur_end=();
		}
		for (my $j=$i; $j<=(length($cur_plastome)-2);$j+=3){ 
			my $curcodon=substr($cur_plastome,$j,3);
			if(exists $codons{$curcodon}){
				if($current_orf){
					if($codons{$curcodon} eq "_"){
						$current_orf.=$curcodon;
						$cur_end = $j+2;
						if(length($current_orf)>18){
							$orf_count++;
							my $torf = "orf" . $orf_count;
							$orfs{$torf}=$cur_start . "-". $cur_end;
							
						}
						$current_orf=();
						$cur_start=();
						$cur_end=();

					}
					else{
						$current_orf.=$curcodon;
						$cur_end = $j+2;
					}
				}	
				else{
					if($codons{$curcodon} eq "_"){
						$current_orf=();
						$cur_start=();
						$cur_end=();
					}
					else{
						$current_orf = $curcodon;
						$cur_start=$j;
						$cur_end=$j+2;
					}
				}
			}
			else{
				if($current_orf){
					if(length($current_orf)> 18){
						$orf_count++;
						my $torf = "orf" . $orf_count;
						$orfs{$torf}=$cur_start . "-". $cur_end;
						
					}
					$current_orf=();
					$cur_start=();
					$cur_end=();
				}
			}

		}
	}
	return(%orfs)
}

