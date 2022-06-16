#!/usr/bin/perl -w
#use strict;
$ARGV[0] =~ /(.*?)\//;
my $prefix = $1;
my $cp_genome;
open my $file, "<", $ARGV[0]; #cp genome full
while(<$file>){
    chomp;
    if(/>/){
        next;
    }
    else{
        $cp_genome .= $_;
    }
}

my $lsc;
my $ssc;
my $irb;
my $ira;
my $full = length($cp_genome);

my $boundary1= substr(reverse($cp_genome), 0, 20);
$boundary1 =~ tr/ATCGatcg/TAGCtagc/;
my $start_irb = index($cp_genome, $boundary1);

$lsc = substr($cp_genome, 0, $start_irb);

my $i = $start_irb;

my $irb_end = substr($cp_genome, $i, 20);
$irb_end = reverse($irb_end);
$irb_end =~ tr/ATCGatcg/TAGCtagc/;

until($cp_genome !~ /$irb_end/){
    $i++;
    $irb_end = substr($cp_genome, $i, 20);
    $irb_end = reverse($irb_end);
    $irb_end =~ tr/ATCGatcg/TAGCtagc/;
}
$i--;
$irb_end = substr($cp_genome, $i, 20);
$irb = substr($cp_genome, $start_irb, ($i+20-$start_irb));
my $ira_seq = $irb_end;
$ira_seq = reverse($ira_seq);
$ira_seq =~ tr/ATCGatcg/TAGCtagc/;
my $ira_start = index($cp_genome, $ira_seq);
$ssc = substr($cp_genome, $i+20, $ira_start-($i+20));
$ira= substr($cp_genome, $ira_start);
my ($lsc_range, $ssc_range, $ira_range, $irb_range);
$lsc_range = (index($cp_genome, $lsc)+1) . "-" . (index($cp_genome, $lsc)+length($lsc));
$ssc_range = (index($cp_genome, $ssc)+1) . "-" . (index($cp_genome, $ssc)+length($ssc));
$ira_range = (index($cp_genome, $ira)+1) . "-" . (index($cp_genome, $ira)+length($ira));
$irb_range = (index($cp_genome, $irb)+1) . "-" . (index($cp_genome, $irb)+length($irb));

my $regions = 1;
if(!$irb || !$ira || !$ssc || !$lsc){
	$regions=();
}
#GFF3 to Verdant format parser
#Usage: 1-GFF3 file
open my $out, ">", $ARGV[1] . "_genbank_annotation.vt";
my $outid = $ARGV[1] . "_genbank_annotation.vt";
my %final_annotation;
my %temp_annotation;
my %direction;
my $plastome;
open $file, "<", $ARGV[1];
while(<$file>){
		chomp;
		my @tarray = split/\s+/;
		if(@tarray <2){
				next;
		}
		if(/^\#/){
			next;
		}
		if($tarray[2] =~ /region/i){
			my $id;
			if($regions && (/LSC/i || /SSC/i || /IR/i)){
                                        next;
                        }
			  if($tarray[8] =~ /gene/i){
                                        $tarray[8] =~ /gene\=(\w+)\;/i;
                                        $id = $1;
                                }
                                elsif($tarray[8] =~ /Name/i){
                                        $tarray[8] =~ /Name\=(\w+)\;/i;
                                        $id = $1;
                                }
                                if(/LSC/i || /SSC/i || /IR/i){
                                        if(/LSC/i){
                                                $lsc_range = $tarray[3] . "-" . $tarray[4];
                                        }
                                        elsif(/SSC/i){
                                                $ssc_range = $tarray[3] . "-" . $tarray[4];
                                        }
                                        elsif(/IRA/i){
                                                $ira_range = $tarray[3] . "-" . $tarray[4];
                                        }
                                        elsif(/IRB/i){
                                                $irb_range = $tarray[3] . "-" . $tarray[4];
                                        }
                                }
                                else{

                                        $temp_annotation{$id}{$tarray[3]}{$tarray[4]}=$tarray[6];
                                        $direction{$id}{$tarray[6]}{$tarray[3]}=$tarray[4];
                                }
		}
		
		if($tarray[2] =~ /^CDS$/i){
				my $id;
				if($regions && (/LSC/i || /SSC/i || /IR/i)){
					next;
				}
				if($tarray[8] =~ /gene/i){
					$tarray[8] =~ /gene\=(\w+)/i;
					$id = $1;
				}
				elsif($tarray[8] =~ /Name/i){
					$tarray[8] =~ /Name\=(\w+)/i;
					$id = $1;
				}
				if(/LSC/i || /SSC/i || /IR/i){
					if(/LSC/i){
						$lsc_range = $tarray[3] . "-" . $tarray[4];
					}
					elsif(/SSC/i){	
						$ssc_range = $tarray[3] . "-" . $tarray[4];
					}
					elsif(/IRA/i){
						$ira_range = $tarray[3] . "-" . $tarray[4];
					}
					elsif(/IRB/i){
						$irb_range = $tarray[3] . "-" . $tarray[4];
					}
				}
				else{
					
					$temp_annotation{$id}{$tarray[3]}{$tarray[4]}=$tarray[6];
					$direction{$id}{$tarray[6]}{$tarray[3]}=$tarray[4];
				}
		}
		if($tarray[2] =~ /^tRNA$/ ){
				my $id;
				if($tarray[8] =~ /gene/){
					$tarray[8] =~ /gene\=(trn\w+\-\w+)/;
					$id = $1;
				}
				elsif($tarray[8] =~ /Name/i){
					$tarray[8] =~ /Name\=(trn\w+-\w+)/i;
                                        $id = $1;
				}
				unless($id){
					close $out;
				}
				
				my $skip;
				if(exists $temp_annotation{$id}{$tarray[3]}){
					for my $tempend (keys %{$temp_annotation{$id}{$tarray[3]}}){
						if ($tempend > $tarray[4]){
							delete $temp_annotation{$id}{$tarray[3]}{$tempend};
						}
						else {
							$skip =1;
						}
					}
				}
				if($skip){
					next;
				}
				$temp_annotation{$id}{$tarray[3]}{$tarray[4]}=$tarray[6];
				$direction{$id}{$tarray[6]}{$tarray[3]}=$tarray[4];
				
		}
		if($tarray[2] =~ /^rRNA$/i){
				my $id;
				
                                if($tarray[8] =~ /gene/i){
                                        $tarray[8] =~ /gene\=(\w+)/;
                                        $id = $1;
                                }
                                elsif($tarray[8] =~ /Name/i){
                                        $tarray[8] =~ /Name\=(\w+)/i;
                                        $id = $1;
                                }
				$temp_annotation{$id}{$tarray[3]}{$tarray[4]}=$tarray[6];
				$direction{$id}{$tarray[6]}{$tarray[3]}=$tarray[4];
		}
}

for my $geneid (keys %direction){
   for my $dirtrue (keys %{$direction{$geneid}}){	
	my %temp_pos;
	my $posmin=10000000;
	my $posmax=-1;
	my $minmin=10000000;
	my $minmax=-1;
	
	for my $start (keys %{$temp_annotation{$geneid}}){
		for my $stop (keys %{$temp_annotation{$geneid}{$start}}){
			unless($temp_annotation{$geneid}{$start}{$stop} eq $dirtrue){
				next;
			}
			if(keys %{$direction{$geneid}{$temp_annotation{$geneid}{$start}{$stop}}} > 1){
				$temp_pos{$start}=$stop;
				if($temp_annotation{$geneid}{$start}{$stop} eq "+"){
					if($start>$posmax){
						$posmax=$start;
					}
					if($start<$posmin){
						$posmin=$start;
					}
					if($stop>$posmax){
						$posmax=$stop;
					}
					if($stop<$posmin){
						$posmin=$stop;
					}
				}
				if($temp_annotation{$geneid}{$start}{$stop} eq "-"){
					if($start>$minmax){
						$minmax=$start;
					}
					if($start<$minmin){
						$minmin=$start;
					}
					if($stop>$minmax){
						$minmax=$stop;
					}
					if($stop<$minmin){
						$minmin=$stop;
					}
				}			
			
		
		
				for my $tstart (keys %temp_pos){
					if($dirtrue eq "+"){
						if($tstart != $posmin || $temp_pos{$tstart} != $posmax){

						if($temp_annotation{$geneid}{$tstart}{$temp_pos{$tstart}} eq "+"){
							my $exonnumber = keys %{$direction{$geneid}{"+"}};
							if($tstart == $posmin){
								$final_annotation{$tstart}{$temp_pos{$tstart}}{$temp_annotation{$geneid}{$tstart}{$temp_pos{$tstart}}}=$geneid . "_exon1";
							}
							elsif($temp_pos{$tstart} == $posmax){
								$final_annotation{$tstart}{$temp_pos{$tstart}}{$temp_annotation{$geneid}{$tstart}{$temp_pos{$tstart}}}=$geneid . "_exon" . $exonnumber;
							}
							else{
								$final_annotation{$tstart}{$temp_pos{$tstart}}{$temp_annotation{$geneid}{$tstart}{$temp_pos{$tstart}}}=$geneid . "_exon2";
							}
						
						}
						}
					}
					elsif($tstart != $minmin || $temp_pos{$tstart} != $minmax){
							my $exonnumber = keys %{$direction{$geneid}{"-"}};
							if($tstart == $minmin){
								$final_annotation{$tstart}{$temp_pos{$tstart}}{$temp_annotation{$geneid}{$tstart}{$temp_pos{$tstart}}}=$geneid . "_exon" . $exonnumber;
							}
							elsif($temp_pos{$tstart} == $minmax){
								$final_annotation{$tstart}{$temp_pos{$tstart}}{$temp_annotation{$geneid}{$tstart}{$temp_pos{$tstart}}}=$geneid . "_exon1";
							}
							else{
								$final_annotation{$tstart}{$temp_pos{$tstart}}{$temp_annotation{$geneid}{$tstart}{$temp_pos{$tstart}}}=$geneid . "_exon2";
							}
						}
					
					}
				}
			
			else{
		
				$final_annotation{$start}{$stop}{$temp_annotation{$geneid}{$start}{$stop}}=$geneid;
			}
		}		
		
	}
}}
my $prev;
my $prev_end;
for my $start (sort {$a <=> $b} keys %final_annotation){
	
	for my $end (sort {$a <=> $b} keys %{$final_annotation{$start}}){
		
		for my $dir (keys %{$final_annotation{$start}{$end}}){
			my $tids = $final_annotation{$start}{$end}{$dir};
			if(!$prev){
					$final_annotation{1}{$start-1}{"+"}="Start\~" . $tids;
					$prev = $tids;
					$prev_end=$end;
			}
			else{
				if($start-$prev_end <=0){
					$prev = $tids;
					$prev_end=$end;
				}
				elsif ($prev =~ /exon/ && $tids =~ /exon/){
					$prev =~ /(.+)_exon/;
					my $pgene = $1;
					$tids =~ /(.+)_exon/;
					my $tgene = $1;
	
					if($tgene eq $pgene){
					  if($tgene eq "rps12"){
						$final_annotation{$prev_end+1}{$start-1}{"+"}=$tgene . "_intron1";
                                                $prev = $final_annotation{$start}{$end}{$dir};
                                                $prev_end = $end;
					  }
				       	  else{
						$prev =~ /_exon(\d)/;
						my $pnum = $1;
						$tids =~ /_exon(\d)/;
						my $tnum = $1;
						if ($pnum == 3 && $tnum == 2){
							$final_annotation{$prev_end+1}{$start-1}{"+"}=$tgene . "_intron2";
							$prev = $final_annotation{$start}{$end}{$dir};
							$prev_end = $end;
						}
						elsif ($pnum ==2 && $tnum ==1) { 
							$final_annotation{$prev_end+1}{$start-1}{"+"}=$tgene . "_intron1";
							$prev = $final_annotation{$start}{$end}{$dir};
							$prev_end = $end;
						}
						elsif ($pnum ==1 && $tnum ==2) { 
							$final_annotation{$prev_end+1}{$start-1}{"+"}=$tgene . "_intron1";
							$prev = $final_annotation{$start}{$end}{$dir};
							$prev_end = $end;
						}
						elsif ($pnum ==2 && $tnum ==3) { 
							$final_annotation{$prev_end+1}{$start-1}{"+"}=$tgene . "_intron2";
							$prev = $final_annotation{$start}{$end}{$dir};
							$prev_end = $end;
						}
					}}
					else{
						$final_annotation{$prev_end+1}{$start-1}{"+"}=$prev . "~" . $final_annotation{$start}{$end}{$dir};
						$prev = $final_annotation{$start}{$end}{$dir};
						$prev_end = $end;
					}
				}
				else{
					$final_annotation{$prev_end+1}{$start-1}{"+"}=$prev . "~" . $final_annotation{$start}{$end}{$dir};
					$prev = $final_annotation{$start}{$end}{$dir};
					$prev_end = $end;
				}
			}
		}
	}
}

$final_annotation{$prev_end+1}{$plastome}{"+"}=$prev . "~End";


for my $start (sort {$a <=> $b} keys %final_annotation){
	for my $end (keys %{$final_annotation{$start}}){
		for my $dir (keys %{$final_annotation{$start}{$end}}){
			print $out "$final_annotation{$start}{$end}{$dir}\t$start\t$end\t$dir\n";
		}
	}
}
$lsc_range =~ /(\d+)-(\d+)/;
print $out "LSC\t$1\t$2\t+\n";
$irb_range =~ /(\d+)-(\d+)/;
print $out "IRB\t$1\t$2\t+\n";
$ssc_range =~ /(\d+)-(\d+)/;
print $out "SSC\t$1\t$2\t+\n";
$ira_range =~ /(\d+)-(\d+)/;
print $out "IRA\t$1\t$2\t+\n";
print $out "FULL\t1\t$full\t+\n";
#`mv $outid $ARGV[1]`;
