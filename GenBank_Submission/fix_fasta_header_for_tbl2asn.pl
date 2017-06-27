#!/usr/bin/perl -w

my @tbl_files = <tbl_files/*.tbl>;

for my $tbl (@tbl_files){
		$tbl =~ /\/(.*?).tbl/;
		my $base_name = $1;

		my $line1 = `head -n 1 $tbl`;
		my @name = split(/\s+/, $line1);
		open my $out, ">", $base_name . ".fsa";
		open my $seqfile, "<", "seq_files/" . $base_name . ".fsa";
		while(<$seqfile>){
			if(/>/){
				print $out ">$name[1] $name[2] $name[3] [location=chloroplast]\n";
			}	
			else{
				print $out "$_";
			}
		}
		close $out;


}
