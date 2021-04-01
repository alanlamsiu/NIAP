#!/bin/perl
#Author: ALin
#Purpose: To combine the intron gtf as indicated by a list. The resulted file will show whether an exon-exon junction is present in a sample or not. The list should contain the file name only.
#Usage: perl NIAP_merge_intron_v1.1.pl <input list> <output>
#Change log:
#	v1.1	2021-04	The script is a copy of combine_intron_gtf_v1.1.pl.

if((@ARGV != 2) && (@ARGV != 3)){
	print "Usage: perl NIAP_merge_intron_v1.1.pl <input list> <output> <path to bedtools>\n";
	exit;
}

use strict;

my $list = shift @ARGV;
my $out = shift @ARGV;
my $btp = shift @ARGV;

unless($btp){
	$btp = "bedtools";
}

unless(`$btp`){
	print "Incorrect path for bedtools!\n";
	print "Usage: perl combine_intron_gtf_v1.1.pl <input list> <output> <path to bedtools>\n";
	exit;
}

open(LIST, $list) or die "Cannot open $list!\n";
open(OUT, ">temp.gtf") or die "Cannot create $out!\n";

my %intron = ();

while(<LIST>){
	my $file = $_;
	chomp $file;
	open(IN, $file) or die "Cannot open $file!\n";
	while(<IN>){
		my $line = $_;
		chomp $line;
		my @line = split('\t', $line);
		$line[8] =~ /intron_id \"([^\;]+)\"/;
		my $id = $1;
		unless(exists $intron{$id}){
			@{$intron{$id}} = (0, $line[1]);
		}
		$line[8] =~ /depth \"([^\;]+)\"/;
		my $depth = $1;
		$intron{$id}[0] += $depth;
	}
	close IN;
}
close LIST;

my $new_line_flag = 0;

foreach my $id (keys %intron){
	if($new_line_flag == 0){
		$new_line_flag = 1;
	}
	else{
		print OUT "\n";
	}
	$id =~ /^(.+)_(\d+)-(\d+)_([+-]$)/;
	my ($chr, $left, $right, $strand) = ($1, $2, $3, $4);
	unless($chr && $left && $right && $strand){
		print "Incorrect intron ID: $id!\n";
		exit;
	}
	print OUT "$chr\t$intron{$id}[1]\tintron\t$left\t$right\t.\t$strand\t.\tintron_id \"$id\"\; depth \"$intron{$id}[0]\"\;";
}
close OUT;

system("bedtools sort -i temp.gtf > $out");
system("rm temp.gtf");












