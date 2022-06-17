#!/usr/bin/perl -w
use strict;

my %seqs;
my $sid;
open my $file, "<", $ARGV[0];
while(<$file>){
	s/\r//g;
	chomp;
	if(/>/){
		$sid = $_;
	}
	else{
		$seqs{$sid}.=$_;
	}
}

open my $out, ">", $ARGV[0] . "_single.fa";
for my $si (keys %seqs){
	print $out "$si\n$seqs{$si}\n";
}
