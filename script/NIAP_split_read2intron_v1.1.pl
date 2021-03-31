#!/bin/perl
#Author: ALin
#Purpose: To identify introns (exon-exon junctions) according to split read sites.
#Usage: perl NIAP_split_read2intron_v1.1.pl [option]
#		-bam <String> Input bam file
#		-o <String> Output
#		-s <String> Strand of reads (Default: unstranded; RF for dUTP method; D for direct RNA-seq long-read)
#		-q <Integer> Minimal mapping quality (Inclusive; default: 0)
#		-btp <String> The path to bedtools
#Change log:
#		v1.1	2021-03	This script is a copy of split_reads2intron_gtf_v1.4.pl, except that the source column in the output .gtf file was changed to "NIAP".

use strict;
use Getopt::Long;

my $in;
my $out;
my $lib_strand;
my $mapq = 0;
my $btp;
my $help = 0;

GetOptions(
	'bam=s'	=>	\$in,
	'o=s'	=>	\$out,
	's=s'	=>	\$lib_strand,
	'btp=s'	=>	\$btp,
        'q=i'   =>      \$mapq,
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

open(BAM, "samtools view $in |") or die "Cannot open $in!\n";
open(OUT, ">temp.gtf") or die "Cannot create temp.gtf!\n";

my %intron = ();

OUTTER:while(<BAM>){
	chomp;
	my $line = $_;
	my @line = split('\t', $line);
	unless(exists $intron{$line[2]}){
		%{$intron{$line[2]}} = ();
	}
	unless($line[5] =~ /\d+N/){
		next OUTTER;
	}
#	print "$line[0]\t$line[1]\n";
	if($line[4] < $mapq || $line[2] eq "*"){
		next OUTTER;
	}
#	unless($line[1] & 0x2){
		#Not properly paired
#		next OUTTER;
#	}
#	if($line[1] & 0x100){
		#Not primary alignment
#		next OUTTER;
#	}
	if($line[1] & 0x400){
		#Poor reads
		next OUTTER;
	}
	if($line[1] & 0x800){
		#Supplementary alignment
		next OUTTER;
	}
	my $strand;
	if($lib_strand eq "FR"){
#		print "FR\n";
		if($line[1] & 0x40){
			if($line[1] & 0x10){
				$strand = "-";
			}
			else{
				$strand = "+";
			}
		}
		elsif($line[1] & 0x80){
			if($line[1] & 0x10){
				$strand = "+";
			}
			else{
				$strand = "-";
			}
		}
	}
	elsif($lib_strand eq "RF"){
#		print "RF\n";
		if($line[1] & 0x40){
                        if($line[1] & 0x10){
				$strand = "+";
			}
			else{
				$strand = "-";
			}
		}
		elsif($line[1] & 0x80){
			if($line[1] & 0x10){
				$strand = "-";
			}
			else{
				$strand = "+";
			}
		}
	}
	elsif($lib_strand eq "D"){
		if($line[1] & 0x10){
			$strand = "-";
		}
		else{
			$strand = "+";
		}
	}
	else{
		$strand = ".";
	}
#	print "$strand\n";
	my @cigar = ();
	my $pos = $line[3] - 1;
	while($line[5] =~ /^(\d+)(\w)/){
		my $temp_length = $1;
		my $status = $2;
		push(@cigar, $status);
		push(@cigar, $temp_length);
		$line[5] =~ s/^(\d+)(\w)//;
	}
	my $num_cigar = @cigar;
#	print "$line[2]\n$line[0]\t$pos\n";
	for(my $i = 0; $i < $num_cigar; $i+=2){
#		print "$cigar[$i]$cigar[$i+1]\t";
		if($cigar[$i] eq "N"){
			#First deal with the left exon
			my $pos1 = $pos + 1;
			unless(exists $intron{$line[2]}{$pos1}){
				%{$intron{$line[2]}{$pos1}} =();
			}
			#Then deal with the right exon
			$pos += $cigar[$i+1];
			my $pos2 = $pos;
			unless(exists $intron{$line[2]}{$pos1}{$pos2}){
				%{$intron{$line[2]}{$pos1}{$pos2}} =();
			}
			unless(exists $intron{$line[2]}{$pos1}{$pos2}{$strand}){
				$intron{$line[2]}{$pos1}{$pos2}{$strand} = 0;
			}
			$intron{$line[2]}{$pos1}{$pos2}{$strand}++;
#			print "$line[2]\_$pos1\-$pos2\_$strand\t$intron{$line[2]}{$pos1}{$pos2}{$strand}\n";
		}
		elsif(($cigar[$i] eq "M") || ($cigar[$i] eq "D")){
			$pos += $cigar[$i+1];
		}
	}
}
close BAM;

my $new_line_flag = 0;

foreach my $ref (keys %intron){
	foreach my $pos1 (keys %{$intron{$ref}}){
		foreach my $pos2 (keys %{$intron{$ref}{$pos1}}){
			foreach my $strand (keys %{$intron{$ref}{$pos1}{$pos2}}){
				if($new_line_flag == 0){
					$new_line_flag = 1;
				}
				else{
					print OUT "\n";
				}
				if($strand eq "."){
					print OUT "$ref\tNIAP\tintron\t$pos1\t$pos2\t.\t+\t.\tintron_id \"" . $ref . "_" . $pos1 . "-" . $pos2 . "_+\"\; depth \"$intron{$ref}{$pos1}{$pos2}{$strand}\"\;";
					print OUT "\n$ref\tNIAP\tintron\t$pos1\t$pos2\t.\t-\t.\tintron_id \"" . $ref . "_" . $pos1 . "-" . $pos2 . "_-\"\; depth \"$intron{$ref}{$pos1}{$pos2}{$strand}\"\;";
				}
				else{
					print OUT "$ref\tNIAP\tintron\t$pos1\t$pos2\t.\t$strand\t.\tintron_id \"" . $ref . "_" . $pos1 . "-" . $pos2 . "_" . $strand . "\"\; depth \"$intron{$ref}{$pos1}{$pos2}{$strand}\"\;";
				}
			}
		}
	}
}
close OUT;
system("bedtools sort -i temp.gtf > $out");
system("rm temp.gtf");

sub print_usage{
	print "Usage: perl NIAP_split_read2intron_v1.1.pl [option]\n\t-bam <String> Input bam file\n\t-o <String> Output\n\t-s <String> Strand of reads (Default: unstranded; RF for dUTP method; D for direct RNA-seq long-read)\n\t-q <Integer> Minimal mapping quality (Inclusive; default: 0)\n\t-btp <String> The path to bedtools\n";
}





