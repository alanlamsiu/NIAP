#!/bin/perl
#Author: ALin
#Purpose: Convert an exon gtf to an intron gtf.
#Usage: perl NIAP_exon2intron_v1.1.pl [option]
#		-i <String> Input exon gtf
#		-o <String> Basename of output files
#		-btp <String> The path to bedtools
#Change log
#		v1.1	2021-04	It is the same as exon2intron_gtf_v1.2.pl.
#
#
use strict;
use Getopt::Long;

my $version = "v1.1";

my $in;
my $out;
my $btp;
my $help = 0;

GetOptions(
	'i=s'	=>	\$in,
	'o=s'	=>	\$out,
	'btp=s'	=>	\$btp,
        'h!'    =>      \$help,
);

if($help){
	print_usage();
	exit;
}

unless($btp){
	$btp = "bedtools";
}

unless(`$btp`){
	print "Incorrect path for bedtools!\n";
	print_usage();
	exit;
}


unless($in && $out){
	print_usage();
	exit;
}

open(IN, $in) or die "Cannot open $in!\n";
open(OUT, ">temp.gtf") or die "Cannot create temp.gtf!\n";
open(LEN, ">$out\_length.txt") or die "Cannot create $out\_length.txt!\n";

my %exon = ();
my %exon_count = ();
my %intron = ();
my %gene_info = ();

while(<IN>){
	chomp;
	my $line = $_;
	if($line =~ /^#/){
		next;
	}
	my @line = split('\t', $line);
	if($line[2] ne "exon"){
		next;
	}
	$line[8] =~ /gene_id \"([^\"]+)\"/;
	my $gene_id = $1;
	$line[8] =~ /gene_name \"([^\"]+)\"/;
	my $gene_name = $1;
	$line[8] =~ /gene_[^\"]*type \"([^\"]+)\"/;
	my $gene_biotype = $1;
	$line[8] =~ /transcript_id \"([^\"]+)\"/;
	my $transcript_id = $1;
	unless(exists $exon{$gene_id}){
		%{$exon{$gene_id}} = ();
		%{$exon_count{$gene_id}} = ();
		%{$intron{$gene_id}} = ();
		@{$gene_info{$gene_id}} = ($line[0], $line[1], $gene_name, $gene_biotype);
	}
	unless(exists $exon{$gene_id}{$line[6]}){
		%{$exon{$gene_id}{$line[6]}} = ();
		%{$exon_count{$gene_id}{$line[6]}} = ();
		%{$intron{$gene_id}{$line[6]}} = ();
	}
	unless(exists $exon{$gene_id}{$line[6]}{$transcript_id}){
		@{$exon{$gene_id}{$line[6]}{$transcript_id}} = ();
		$exon_count{$gene_id}{$line[6]}{$transcript_id} = 0;
	}
	$line[3] -= 1;
	$line[4] += 1;
	@{$exon{$gene_id}{$line[6]}{$transcript_id}[$exon_count{$gene_id}{$line[6]}{$transcript_id}]} = ($line[3], $line[4]);
	$exon_count{$gene_id}{$line[6]}{$transcript_id}++;
}
close IN;

foreach my $gene_id (keys %exon){
	STRAND:foreach my $strand (keys %{$exon{$gene_id}}){
		TRANSCRIPT:foreach my $transcript_id (keys %{$exon{$gene_id}{$strand}}){
			if($exon_count{$gene_id}{$strand}{$transcript_id} == 1){
				next TRANSCRIPT;
			}
			my $num_exon = $exon_count{$gene_id}{$strand}{$transcript_id};
			my @sorted_exon = sort {$a->[0] <=> $b->[0]} @{$exon{$gene_id}{$strand}{$transcript_id}};
			for(my $i = 0; $i <= $num_exon - 2; $i++){
				unless(exists $intron{$gene_id}{$strand}{$sorted_exon[$i][1]}){
					%{$intron{$gene_id}{$strand}{$sorted_exon[$i][1]}} = ();
				}
				unless(exists $intron{$gene_id}{$strand}{$sorted_exon[$i][1]}{$sorted_exon[$i+1][0]}){
					$intron{$gene_id}{$strand}{$sorted_exon[$i][1]}{$sorted_exon[$i+1][0]} = 1;
				}
			}
		}
	}
}

my $new_line_flag_gtf = 0;
my $new_line_flag_len = 0;

foreach my $gene_id (keys %intron){
	foreach my $strand (keys %{$intron{$gene_id}}){
		my $length = 0;
		foreach my $pos1 (keys %{$intron{$gene_id}{$strand}}){
			foreach my $pos2 (keys %{$intron{$gene_id}{$strand}{$pos1}}){
				my $temp_length = $pos2 - $pos1 + 1;
				$length += $temp_length;
				if($new_line_flag_gtf == 0){
					$new_line_flag_gtf = 1;
				}
				else{
					print OUT "\n";
				}
				print OUT "$gene_info{$gene_id}[0]\t$gene_info{$gene_id}[1]\tintron\t$pos1\t$pos2\t.\t$strand\t.\tintron_id \"" . $gene_info{$gene_id}[0] . "_" . $pos1 . "-" . $pos2 . "_" . $strand . "\"\; gene_id \"$gene_id\"\; gene_name \"$gene_info{$gene_id}[2]\"\; gene_type \"intron\"\; gene_type_original \"$gene_info{$gene_id}[3]\"\;";
			}
		}
		if($new_line_flag_len == 0){
			$new_line_flag_len = 1;
		}
		else{
			print LEN "\n";
		}
		print LEN "$gene_id\t$strand\t$length";
	}
}
close OUT;
system("bedtools sort -i temp.gtf > $out\_sorted.gtf");
system("rm temp.gtf");

sub print_usage{
	print "Usage: perl NIAP_exon2intron_${version}.pl [option]\n\t-i <String> Input exon gtf file\n\t-o <String> Basename of output files\n\t-btp <String> The path to bedtools\n";
}





