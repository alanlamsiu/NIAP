#!/bin/perl
#Author: ALin
#Purpose: To collapse RNA-seq long reads, perform splice junction correction based on the nearest splice junction provided by NGS or reference if available, and generate a .gtf file for all transcripts with unique structures.
#Usage: perl NIAP_v1.1.pl [option]
#		-bam <String> Input bam file
#		-junc <String> The gtf file specifying splice junctions used for correction. An optional attribute, depth, can be included.
#		-injunc <Boolean> The intervals of splice junctions specified by -junc are introns. (Default: false for exons as input)
#		-minjuncdist <Int> The minimal distance for splice junction correction (Default: 100 bp, empirical)
#		-juncdist <Boolean> Output the distance of splice junction between the prediction and the reference (Default: false)
#		-enddist <Boolean> Output the distribution of the 5' and 3' end position for each transcript (Default: false)
#		-polya <String> The output from nanopolish polya
#		-prefix <String> The prefix of gene and transcirpt IDs (Default: NIAP)
#		-o <String> Output base
#		-q <Integer> Minimal mapping quality (Inclusive; default: 0)
#		-stp <String> The path to samtools
#		-t <Int> The number of threads to use (Default: 1)
#		-btp <String> The path to bedtools
#Change log:
#	v1.1		This code is adopted from GLRA_v1.1.pl .Comparison between the 5'/3' ends of the read data and those in a referenc annotation is no longer supported. Instead, it is not possible to investigate the distribution of 5' and 3' end position for each transcript. Poly-A information from Nanopolish can be incorporated. It will also consider the coverage of the reference exon-exon junction if applicable.

use strict;
use Getopt::Long;
use List::Util qw(max min);
use Benchmark qw(:hireswallclock);

my $total_starttime = Benchmark->new;
my $total_finishtime;
my $random = `tr -dc A-Za-z0-9 </dev/urandom | head -c 10`;
my $starttime;
my $finishtime;
my $timespent;

my $bam;
my $junc;
my $injunc;
my $minjuncdist = 100;
my $juncdist;
my $enddist;
my $polya;
my $prefix = "NIAP";
my $out;
my $mapq = 0;
my $stp;
my $thread = 1;
my $btp;
my $help;

GetOptions(
	'bam=s'	=>	\$bam,
	'junc=s'	=>	\$junc,
	'injunc!'	=>	\$injunc,
	'minjuncdist=i'	=>	\$minjuncdist,
	'juncdist!'	=>	\$juncdist,
	'enddist!'	=>	\$enddist,
	'polya=s'	=>	\$polya,
	'prefix=s'	=>	\$prefix,
	'o=s'	=>	\$out,
	'q=i'	=>	\$mapq,
	'stp=s'	=>	\$stp,
	't=i'	=>	\$thread,
	'btp=s'	=>	\$btp,
        'h!'    =>      \$help,
);

print "Start processing $bam...\n";
system("date");

if($help){
        print_usage();
        exit;
}

unless($stp){
	$stp = "samtools";
}

unless(`which $stp`){
	print "Incorrect path for samtools!\n";
	print_usage();
	exit;
}

unless($btp){
	$btp = "bedtools";
}

unless(`which $btp`){
	print "Incorrect path for bedtools!\n";
	print_usage();
	exit;
}

unless($bam && $out){
	print_usage();
	exit;
}

if($juncdist){
	unless($junc){
		print "Must specify -junc with -juncdist!\n";
		print_usage();
		exit;
	}
}

my $total_thread = `grep -c \"^processor\" /proc/cpuinfo`;

if($thread > $total_thread){
	print "Exceeding the maximal total number of threads! Setting the number of thread to $total_thread...";
}

my $out_gtf = $out . ".gtf";
my $temp_out_gtf = "temp_" . $random . ".gtf";
my $out_stat = $out . "_stat.txt";
my $out_junc_dist = $out . "_junc_dist.txt";
my $out_end_dist = $out . "_end_dist.txt";

my %ref_eej_array = ();
my %ref_eej_hash = ();
my %polya = ();

my $gtf_fh;
my $enddist_fh;

open($gtf_fh, ">$temp_out_gtf") or die "Cannot create $temp_out_gtf!\n";

print $gtf_fh "#perl NIAP_v1.1.pl -bam $bam -o $out -q $mapq -minjuncdist $minjuncdist -t $thread -stp $stp -btp $btp";

if($junc){
	$starttime = Benchmark->new;
	read_ref_eej($junc, \%ref_eej_array, \%ref_eej_hash);
	$finishtime = Benchmark->new;
	$timespent = timediff($finishtime, $starttime);
	print "Used ".timestr($timespent)." reading the reference junction.\n";
	print $gtf_fh " -junc $junc";
}
if($injunc){
	print $gtf_fh " -injunc";
}
if($juncdist){
	print $gtf_fh " -juncdist";
}
if($enddist){
	print $gtf_fh " -enddist";
	open($enddist_fh, ">$out_end_dist") or die "Cannot create $out_end_dist!\n";
        print $enddist_fh "Chr\tStrand\tGene_id\tTranscript_id\t5'\t3'";
}
if($polya){
	$starttime = Benchmark->new;
	read_polya($polya, \%polya);
	$finishtime = Benchmark->new;
	$timespent = timediff($finishtime, $starttime);
	print "Used ".timestr($timespent)." reading the polya information.\n";
	print $gtf_fh " -polya $polya";
}

my %loci;
my $current_chr;
my %juncdist;
my %read;
%{$read{'all'}} = ();
%{$read{'filtered'}} = ();
%{$read{'kept'}} = ();
%{$read{'low_quality'}} = ();
%{$read{'unmapped'}} = ();
%{$read{'failed'}} = ();
%{$read{'duplicated'}} = ();
%{$read{'secondary'}} = ();
%{$read{'supplementary'}} = ();
my $num_filtered_juncdist = 0;

open(BAM, "samtools sort -@ $thread $bam | samtools view |") or die "Cannot open $bam!\n";

while(<BAM>){
	chomp;
	my $line = $_;
	my @line = split('\t', $line);
	#print "$line\n";
	#Filter the read
	if(filter_bam_line(\%read, \@line)){
		next;
	}
	#Parse the read structure
	my ($chr, $strand, $leftmost, $rightmost, $eej, $polya_flag) = parse_bam_line(\@line);
	unless($current_chr){
		$starttime = Benchmark->new;
		$current_chr = $chr;
	}
	#Continue with the same chromosome
	if($current_chr eq $chr){
                add_loci(\%loci, $strand, $leftmost, $rightmost, $eej, $polya_flag);
	}
	#Print for the previous chromosme and intiate the record for the new chromosome
	else{
		merge_loci(\%loci);
		if($junc){
			$num_filtered_juncdist = correct_eej($current_chr, \%loci, \%ref_eej_array, \%juncdist, $num_filtered_juncdist);
		}
		print_gtf($current_chr, \%loci, $gtf_fh, $enddist_fh, $prefix);
		$finishtime = Benchmark->new;
		$timespent = timediff($finishtime, $starttime);
		print "Used " . timestr($timespent) . " processing $current_chr.\n";
		#Reset global variables
		$starttime = Benchmark->new;
		$current_chr = $chr;
		#Continue with the new chromosome
		add_loci(\%loci, $strand, $leftmost, $rightmost, $eej, $polya_flag);
	}
}
close BAM;

if(%loci){
	merge_loci(\%loci);
	if($junc){
		$num_filtered_juncdist = correct_eej($current_chr, \%loci, \%ref_eej_array, \%juncdist, $num_filtered_juncdist);
	}
	print_gtf($current_chr ,\%loci, $gtf_fh, $enddist_fh, $prefix);
	$finishtime = Benchmark->new;
	$timespent = timediff($finishtime, $starttime);
	print "Used " . timestr($timespent) . " processing $current_chr.\n";
}
close $gtf_fh;

if($enddist){
	close $enddist_fh;
}

print_stat(\%read, $out_stat, $bam, $num_filtered_juncdist);

if($juncdist){
	print_juncdist(\%juncdist, $out_junc_dist);
}

system("bedtools sort -i $temp_out_gtf > $out_gtf");
system("rm $temp_out_gtf");

$total_finishtime = Benchmark->new;
$timespent = timediff($total_finishtime, $total_starttime);
print "Used " . timestr($timespent) . " in total.\n";
system("date");

sub print_usage{
	print "Usage: perl NIAP_v1.1.pl [option]
	\t-bam <String> Input bam file
	\t-junc <String> The gtf file specifying splice junctions used for correction. An optional attribute, depth, can be included.
	\t-injunc <Bolean> The intervals of splice junctions specified by -junc are introns. (Default: False for exons as input)
	\t-minjuncdist <Int> The minimal distance for splice junction correction (Default: 100 bp, empirical)
	\t-juncdist <Boolean> Output the distance of splice junction between the prediction and the reference (Default: false)
	\t-enddist <Boolean> Output the distribution of the 5' and 3' end position for each transcript (Default: false)
	\t-polya <String> The output from nanopolish polya
	\t-prefix <String> The prefix of gene and transcirpt IDs (Default: NIAP)
	\t-o <String> Output base
	\t-q <Integer> Minimal mapping quality (Inclusive; default: 0)
	\t-stp <String> The path to samtools
	\t-t <Int> The number of threads to use (Default: 1)
	\t-btp <String> The path to bedtools\n";
}

sub print_stat{
	my $read_ref = $_[0];
	my $temp_out_stat = $_[1];
	my $temp_bam = $_[2];
	my $temp_num_filtered_juncdist = $_[3];
	my $temp_total = keys %{$read_ref->{'all'}};
	my $temp_filtered = keys %{$read_ref->{'filtered'}};
	$temp_filtered += $temp_num_filtered_juncdist;
	my $temp_kept = keys %{$read_ref->{'kept'}};
	$temp_kept -= $temp_num_filtered_juncdist;
	my $temp_low_quality = keys %{$read_ref->{'low_quality'}};
	my $temp_unmapped = keys %{$read_ref->{'unmapped'}};
	my $temp_failed = keys %{$read_ref->{'failed'}};
	my $temp_duplicated = keys %{$read_ref->{'duplicated'}};
	my $temp_secondary = keys %{$read_ref->{'secondary'}};
	my $temp_supplementary = keys %{$read_ref->{'supplementary'}};

	open(TEMP, ">$temp_out_stat") or die "Cannot create $temp_out_stat!\n";
	print TEMP "Summary statistics for $temp_bam
	\tTotal read: $temp_total
	\tRetained read: $temp_kept	
	\tFiltered read: $temp_filtered
	\tLow quality read: $temp_low_quality
	\tUnmapped read: $temp_unmapped
	\tFailed read: $temp_failed
	\tDuplicated read: $temp_duplicated
	\tSecondary alignment: $temp_secondary
	\tSupplementary alignment: $temp_supplementary
	\tIncorrect splice junction: $temp_num_filtered_juncdist\n";
	close TEMP;

	return 1;
}

#Read the splice junction gtf file if provided.
sub read_ref_eej{
	my $temp_junc = $_[0];
	my $ref_eej_array_ref = $_[1];
	my $ref_eej_hash_ref = $_[2];
	open(TEMP, $temp_junc) or die "Cannot open $temp_junc!\n";
	while(<TEMP>){
		chomp;
		my $temp_line = $_;
		$temp_line =~ s/\r//;
		my @temp_line = split('\t', $temp_line);
		unless(exists $ref_eej_array_ref->{$temp_line[0]}){
			%{$ref_eej_array_ref->{$temp_line[0]}} = ();
			%{$ref_eej_hash_ref->{$temp_line[0]}} = ();
		}
		unless(exists $ref_eej_array_ref->{$temp_line[0]}{$temp_line[6]}){
			@{$ref_eej_array_ref->{$temp_line[0]}{$temp_line[6]}} = ();
			%{$ref_eej_hash_ref->{$temp_line[0]}{$temp_line[6]}} = ();
		}
		my $temp_depth = 0;
		if($temp_line[8] =~ /depth \"(\d+)\"\;/){
			$temp_depth = $1;
		}
		if($injunc){
			$temp_line[3] -= 1;
			$temp_line[4] += 1;
		}
		my @temp_ref_eej = ($temp_line[3], $temp_line[4], $temp_depth);
		#print "$temp_line[0]\t$temp_line[6]\t$temp_line[3]\t$temp_line[4]\t$temp_depth\n";
		push(@{$ref_eej_array_ref->{$temp_line[0]}{$temp_line[6]}}, \@temp_ref_eej);
		unless(exists $ref_eej_hash_ref->{$temp_line[0]}{$temp_line[6]}{$temp_line[3]}){
			%{$ref_eej_hash_ref->{$temp_line[0]}{$temp_line[6]}{$temp_line[3]}} = ();
		}
		unless(exists $ref_eej_hash_ref->{$temp_line[0]}{$temp_line[6]}{$temp_line[3]}{$temp_line[4]}){
			$ref_eej_hash_ref->{$temp_line[0]}{$temp_line[6]}{$temp_line[3]}{$temp_line[4]} = $temp_depth;
		}
	}

	return 1;
}

#Read the poly-A result from nanopolish polya
sub read_polya{
	my $temp_polya = $_[0];
	my $polya_ref =$_[1];
	open(TEMP, $temp_polya) or die "Cannot open $temp_polya!\n";
	my $header = 0;
	while(<TEMP>){
		if($header == 0){
			$header = 1;
			next;
		}
		chomp;
		my $temp_line = $_;
		$temp_line =~ s/\r//;
		my @temp_line = split('\t', $temp_line);
		unless(exists $polya_ref->{$temp_line[0]}){
			$polya_ref->{$temp_line[0]} = $temp_line[9];
		}
	}
	close TEMP;
}

#Check whether a single line of the bam file need to be filtered or not.
sub filter_bam_line{
	my $read_ref = $_[0];
	my $temp_line_ref = $_[1];
	my $temp_flag = 0;
	unless(exists $read_ref->{'all'}{$temp_line_ref->[0]}){
		$read_ref->{'all'}{$temp_line_ref->[0]} = 1;	
	}
	if($temp_line_ref->[4] < $mapq){
		#Score filtering
		unless(exists $read_ref->{'low_quality'}{$temp_line_ref->[0]}){
			$read_ref->{'low_quality'}{$temp_line_ref->[0]} = 1;
		}
		$temp_flag = 1;
	}
	if($temp_line_ref->[1] & 0x4){
		#Unmapped read
		unless(exists $read_ref->{'unmapped'}{$temp_line_ref->[0]}){
			$read_ref->{'unmapped'}{$temp_line_ref->[0]} = 1;
		}
		$temp_flag = 1;
	}
	if($temp_line_ref->[1] & 0x200){
		#Failed read
		unless(exists $read_ref->{'failed'}{$temp_line_ref->[0]}){
			$read_ref->{'failed'}{$temp_line_ref->[0]} = 1;
		}
		$temp_flag = 1;
	}
	if($temp_line_ref->[1] & 0x100){
		#Secondary alignment
		unless(exists $read_ref->{'secondary'}{$temp_line_ref->[0]}){
			$read_ref->{'secondary'}{$temp_line_ref->[0]} = 1;
		}
		$temp_flag = 1;
        }
	if($temp_line_ref->[1] & 0x400){
		#PCR or optical duplicate
		unless(exists $read_ref->{'duplicated'}{$temp_line_ref->[0]}){
			$read_ref->{'duplicated'}{$temp_line_ref->[0]} = 1;
		}
		$temp_flag = 1;
        }
	if($temp_line_ref->[1] & 0x800){
		#Supplementary alignment
		unless(exists $read_ref->{'supplementary'}{$temp_line_ref->[0]}){
			$read_ref->{'supplementary'}{$temp_line_ref->[0]} = 1;
		}
		$temp_flag = 1;
	}
	if($temp_flag == 0){
		unless(exists $read_ref->{'kept'}{$temp_line_ref->[0]}){
			$read_ref->{'kept'}{$temp_line_ref->[0]} = 1;
		}
		if(exists $read_ref->{'filtered'}{$temp_line_ref->[0]}){
			delete $read_ref->{'filtered'}{$temp_line_ref->[0]};
		}
	}
	elsif($temp_flag == 1){
		unless(exists $read_ref->{'kept'}{$temp_line_ref->[0]}){
			unless(exists $read_ref->{'filtered'}{$temp_line_ref->[0]}){
				$read_ref->{'filtered'}{$temp_line_ref->[0]} = 1;
			}
		}
	}

	return $temp_flag;
}

#Parse a single line of the bam file and return the chromosome, strand, leftmost position, rightmost position and exon-exon junction string (for multi-exon transcript only).
sub parse_bam_line{
	my $temp_line_ref = $_[0];
	my $temp_chr = $temp_line_ref->[2];
	my $temp_cigar = $temp_line_ref->[5];
	my $temp_strand;
	if($temp_line_ref->[1] & 0x10){
		$temp_strand = "-";
	}
	else{
		$temp_strand = "+";
	}
	my $temp_eej = "";
        my $temp_leftmost = $temp_line_ref->[3];
	my $temp_rightmost = $temp_line_ref->[3] - 1;
	my $temp_polya_flag = "NA";

	if(exists $polya{$temp_line_ref->[0]}){
		$temp_polya_flag = $polya{$temp_line_ref->[0]};
	}
	
	while($temp_cigar =~ /^(\d+)(\w)/){
		my $temp_length = $1;
		my $temp_status = $2;
		if($temp_status eq "N"){
			my $temp_exon_left = $temp_rightmost + $temp_length + 1;
			$temp_eej .= "-${temp_rightmost},${temp_exon_left}";
			$temp_rightmost += $temp_length;
		}
		elsif($temp_status eq "D"){
			my $temp_exon_left = $temp_rightmost + $temp_length + 1;
			if(exists $ref_eej_hash{$temp_chr}{$temp_strand}{$temp_rightmost}{$temp_exon_left}){
				$temp_eej .= "-${temp_rightmost},${temp_exon_left}";
			}
			else{
				$temp_rightmost += $temp_length;
			}
		}
		elsif($temp_status eq "M"){
			$temp_rightmost += $temp_length;
		}
		else{
			if(($temp_status ne "S") && ($temp_status ne "H") && ($temp_status ne "I")){
				print "Unrecognized cigar flag: $temp_status, in the following line\n@{$temp_line_ref}\n";
				exit;
			}
		}
		#print "$temp_length\t$temp_status\t$temp_eej\n";
		$temp_cigar =~ s/^\d+\w//;
	}
	$temp_eej .= "-";
	#print "@{$temp_line_ref}\n$temp_line_ref->[0]\t$temp_chr\t$temp_strand\t$temp_leftmost\t$temp_rightmost\t$temp_eej\t$temp_polya_flag\n";

	return ($temp_chr, $temp_strand, $temp_leftmost, $temp_rightmost, $temp_eej, $temp_polya_flag);
}

#Add a loci to the loci hash as in the first argument
sub add_loci{
	my $loci_ref = $_[0];
	my $temp_strand = $_[1];
	my $temp_leftmost = $_[2];
	my $temp_rightmost = $_[3];
	my $temp_eej = $_[4];
	my $temp_polya_flag = $_[5];
	#Initiate the strand tab for the loci hash
	unless(exists $loci_ref->{$temp_strand}){
		@{$loci_ref->{$temp_strand}} = ();
	}
	my %temp_transcript = ();
	unless(exists $temp_transcript{$temp_eej}){
		@{$temp_transcript{$temp_eej}} = ();
	}
	my %temp_polya_flag = ();
	unless(exists $temp_polya_flag{$temp_polya_flag}){
		$temp_polya_flag{$temp_polya_flag} = 1;
	}
	#print "$temp_polya_flag\t$temp_polya_flag{$temp_polya_flag}\n";
	my %temp_left_end = ();
	unless(exists $temp_left_end{$temp_leftmost}){
		$temp_left_end{$temp_leftmost} = 0;
	}
	$temp_left_end{$temp_leftmost}++;
	my %temp_right_end = ();
	unless(exists $temp_right_end{$temp_rightmost}){
		$temp_right_end{$temp_rightmost} = 0;
	}
	$temp_right_end{$temp_rightmost}++;
	my @temp_transcript = ($temp_leftmost, $temp_rightmost, 1, \%temp_polya_flag, \%temp_left_end, \%temp_right_end);
	#print "@temp_transcript\n";
	push(@{$temp_transcript{$temp_eej}}, \@temp_transcript);
	my @temp_loci = ($temp_leftmost, $temp_rightmost, \%temp_transcript);
	push(@{$loci_ref->{$temp_strand}}, \@temp_loci);
	#print "$current_chr\t$temp_strand\t@temp_loci\n";
	#foreach my $i (keys %{$temp_loci[2]}){
		#print "$i\t$temp_loci[2]{$i}\n";
	#}

	return 1;
}

#Merge locus in the loci hash as in the first argument
sub merge_loci{
	my $loci_ref = $_[0];
	foreach my $temp_strand (keys %{$loci_ref}){
		#print "$current_chr\t$temp_strand\n@{$loci_ref->{$temp_strand}}\n";
		my @temp_loci = sort_loci($loci_ref->{$temp_strand});
		@{$loci_ref->{$temp_strand}} = @temp_loci;
	}
}

#Sort loci
sub sort_loci{
	my ($temp_sref, $temp_start) = @_;
	return if(ref($temp_sref) ne 'ARRAY');

	if(!defined $temp_start){
		if(wantarray){
			my @temp_sets = map {[@{$_}]} @{$temp_sref};
			$temp_sref = \@temp_sets;
		}
		@{$temp_sref} = sort {$a->[0]<=>$b->[0] || $a->[1]<=>$b->[1]} @{$temp_sref};
		$temp_start = 0;
	}
	my $temp_last = $temp_sref->[$temp_start];
	$temp_start++;

	if(@{$temp_last}){
		for(my $i = $temp_start; $i < @{$temp_sref}; $i++){
		my $temp_cur = $temp_sref->[$i];
			if(!@{$temp_cur}){
				next;
			}

			if ($temp_cur->[0] >= $temp_last->[0] && $temp_cur->[0] <= $temp_last->[1] ){
				if ($temp_cur->[1] > $temp_last->[1]){
					$temp_last->[1] = $temp_cur->[1];
				}

				foreach my $temp_cur_eej (keys %{$temp_cur->[2]}){
					if(exists $temp_last->[2]{$temp_cur_eej}){
						push(@{$temp_last->[2]{$temp_cur_eej}}, @{$temp_cur->[2]{$temp_cur_eej}});
						my @temp_new_transcript = merge_transcript($temp_last->[2]{$temp_cur_eej});
						@{$temp_last->[2]{$temp_cur_eej}} = @temp_new_transcript;
					}
					else{
						@{$temp_last->[2]{$temp_cur_eej}} = merge_transcript($temp_cur->[2]{$temp_cur_eej});
					}
				}
				@{$temp_cur} = ();
			}
			else{
				last;
			}
		}
	}
	if($temp_start < @{$temp_sref}){
		sort_loci($temp_sref, $temp_start);
	}
	if(wantarray){
		return sort {$a->[0] <=> $b->[0]} map {@{$_} ? $_ : () } @{$temp_sref};
	}
}

#Merge transcripts with the same exon-exon junction in the same locus
sub merge_transcript{
	my ($temp_sref, $temp_start) = @_;
	if(!defined $temp_start){
		if(wantarray){
			my @temp_sets = map {[@{$_}]} @{$temp_sref};
			$temp_sref = \@temp_sets;
		}
		@{$temp_sref} = sort {$a->[0]<=>$b->[0] || $a->[1]<=>$b->[1]} @{$temp_sref};
		$temp_start = 0;
	}
	my $temp_last = $temp_sref->[$temp_start];
	$temp_start++;

	if(@{$temp_last}){
		#print "@{$temp_last}\n";
		for(my $i = $temp_start; $i < @{$temp_sref}; $i++){
			my $temp_cur = $temp_sref->[$i];
			if(!@{$temp_cur}){
				next;
			}

			if ($temp_cur->[0] >= $temp_last->[0] && $temp_cur->[0] <= $temp_last->[1] ){
				if ($temp_cur->[1] > $temp_last->[1]){
					$temp_last->[1] = $temp_cur->[1];
				}
				$temp_last->[2] += $temp_cur->[2];
				my %temp_last_polya_flag = %{$temp_last->[3]};
				my %temp_cur_polya_flag = %{$temp_cur->[3]};
				foreach my $temp_polya_flag (keys %temp_cur_polya_flag){
					unless(exists $temp_last_polya_flag{$temp_polya_flag}){
						$temp_last_polya_flag{$temp_polya_flag} = 1;
					}
				}
				$temp_last->[3] = \%temp_last_polya_flag;
				my %temp_last_left_end = %{$temp_last->[4]};
				my %temp_cur_left_end = %{$temp_cur->[4]};
				foreach my $temp_left_end (keys %temp_cur_left_end){
					unless(exists $temp_last_left_end{$temp_left_end}){
						$temp_last_left_end{$temp_left_end} = 0;
					}
					$temp_last_left_end{$temp_left_end} += $temp_cur_left_end{$temp_left_end};
				}
				$temp_last->[4] = \%temp_last_left_end;
				my %temp_last_right_end = %{$temp_last->[5]};
				my %temp_cur_right_end = %{$temp_cur->[5]};
				foreach my $temp_right_end (keys %temp_cur_right_end){
					unless(exists $temp_last_right_end{$temp_right_end}){
						$temp_last_right_end{$temp_right_end} = 0;
					}
					$temp_last_right_end{$temp_right_end} += $temp_cur_right_end{$temp_right_end};
				}
				$temp_last->[5] = \%temp_last_right_end;
				@{$temp_cur} = ();
			}
			else{
				last;
			}
		}
	}
	if($temp_start < @{$temp_sref}){
		merge_transcript($temp_sref, $temp_start);
	}
	if(wantarray){
		return sort {$a->[0] <=> $b->[0]} map {@{$_} ? $_ : () } @{$temp_sref};
	}
}

#Find the nearest intervals from a list of a specifed interval
sub find_nearest{
	my @temp_query = @{$_[0]};
	my $temp_db_ref = $_[1];
	my $temp_strand = $_[2];
	my $juncdist_ref = $_[3];
	my $temp_nearest_dist = "NA";
	my @temp_nearest_splice_site = (0, 0);
	my $temp_coverage = 0;
	SEARCH:foreach my $temp_db (@{$temp_db_ref}){
		my $temp_dist = abs($temp_query[0] - $temp_db->[0]) + abs($temp_query[1] - $temp_db->[1]);
		if(($temp_dist < $temp_nearest_dist) || (($temp_dist == $temp_nearest_dist) && ($temp_db->[2] > $temp_coverage)) || ($temp_nearest_dist eq "NA")){
			@temp_nearest_splice_site = @{$temp_db};
			$temp_nearest_dist = $temp_dist;
			$temp_coverage = $temp_db->[2];
		}
	}
	#print "\tNearest overlap @temp_nearest_splice_site of distance $temp_nearest_dist\n";
	unless(exists $juncdist_ref->{$current_chr}){
		%{$juncdist_ref->{$current_chr}} = ();
	}
	unless($juncdist_ref->{$current_chr}{$temp_strand}){
		%{$juncdist_ref->{$current_chr}{$temp_strand}} = ();
	}
	my $temp_query = join(",", @temp_query);
	my $temp_nearest_splice_site = join(",", @temp_nearest_splice_site);
	unless(exists $juncdist_ref->{$current_chr}{$temp_strand}{$temp_query}){
		%{$juncdist_ref->{$current_chr}{$temp_strand}{$temp_query}} = ();
	}
	unless(exists $juncdist_ref->{$current_chr}{$temp_strand}{$temp_query}{$temp_nearest_splice_site}){
		$juncdist_ref->{$current_chr}{$temp_strand}{$temp_query}{$temp_nearest_splice_site} = $temp_nearest_dist;
	}

	if($temp_nearest_dist <= $minjuncdist){
		return \@temp_nearest_splice_site;
	}
	else{
		return 0;
	}
}

#Perform exon-exon junction correction
sub correct_eej{
	my $temp_chr = $_[0];
	my $loci_ref = $_[1];
	my $ref_eej_array_ref = $_[2];
	my $juncdist_ref = $_[3];
	my $temp_num_filtered_juncdist = $_[4];
	foreach my $temp_strand (keys %{$loci_ref}){
		#print "Strand:$temp_strand\n";
		foreach my $temp_loci (@{$loci_ref->{$temp_strand}}){
			#print "Loci:@{$temp_loci}\n";
			EEJ:foreach my $temp_eej (keys %{$temp_loci->[2]}){
				#print "\tEEJ_before:$temp_eej\n\t@{$temp_loci->[2]{$temp_eej}}\n";
				if($temp_eej eq "-"){
					#Skip the correction step for mono-exonic transcript
					#print "\tNo correction for mono-exonic transcript.\n";
					next EEJ;
				}
				my @temp_eej = split('-', $temp_eej);
				#Shift to eliminate the empty value at the first position of the array
				shift @temp_eej;
				my @temp_nearest_splice_site;
				my $temp_corrected_eej = "-";
				my $temp_cur_leftmost = $temp_loci->[2]{$temp_eej}[0][0];
				my $temp_cur_rightmost = $temp_loci->[2]{$temp_eej}[0][1];
				foreach my $temp_splice_site (@temp_eej){
					my @temp_splice_site = split(',', $temp_splice_site);
					#print "\tComparing @temp_splice_site and @{$ref_eej_array_ref->{$current_chr}{$temp_strand}}\n";
					my $temp_nearest_splice_site_ref = find_nearest(\@temp_splice_site, \@{$ref_eej_array_ref->{$temp_chr}{$temp_strand}}, $temp_strand, $juncdist_ref);
					#Can correct one of the splice junction
					if($temp_nearest_splice_site_ref){
						#print "\tFound the nearest splice junction:@{$temp_nearest_splice_site_ref}\n";
						#Valid splice junction
						if(($temp_nearest_splice_site_ref->[0] >= $temp_cur_leftmost) && ($temp_nearest_splice_site_ref->[1] <= $temp_cur_rightmost)){
							#Keep the valid junction
							#print "\tThis is a valid splice junction.\n";
							$temp_corrected_eej .= $temp_nearest_splice_site_ref->[0] . "," . $temp_nearest_splice_site_ref->[1] . "-";
							$temp_cur_leftmost = $temp_nearest_splice_site_ref->[1];
						}
						#Invalid splice junction
						else{
							#print "\tThis is an invalid splice junction with depth $temp_loci->[2]{$temp_eej}[2]\n";
							$temp_num_filtered_juncdist += $temp_loci->[2]{$temp_eej}[0][2];
							delete $temp_loci->[2]{$temp_eej};
							next EEJ;
						}
					}
					#Cannot correct one of the splice junctions
					else{
						#print "\tNot able to correct the splice junction with depth $temp_loci->[2]{$temp_eej}[2].\n";
						$temp_num_filtered_juncdist += $temp_loci->[2]{$temp_eej}[0][2];
						delete $temp_loci->[2]{$temp_eej};
						next EEJ;
					}
				}
				#print "\tEEJ_after:$temp_corrected_eej\n";
				if($temp_eej eq $temp_corrected_eej){
					#print "\tNo correction needed.\n";
					next EEJ;
				}
				#The corrected eej already in the record
				if(exists $temp_loci->[2]{$temp_corrected_eej}){
					#print "\tCorrected EEJ already exists\n";
					#print "\tExisted correct transcript:@{$temp_loci->[2]{$temp_corrected_eej}}\n\tTranscript to be corrected:";
					push(@{$temp_loci->[2]{$temp_corrected_eej}}, @{$temp_loci->[2]{$temp_eej}});
					my @temp_new_transcript = merge_transcript($temp_loci->[2]{$temp_corrected_eej});
					@{$temp_loci->[2]{$temp_corrected_eej}} = @temp_new_transcript;
				}
				#The corrected eej not in the record
				else{
					#print "\tCorrected EEJ not exist\n";
					@{$temp_loci->[2]{$temp_corrected_eej}} = @{$temp_loci->[2]{$temp_eej}};
				}
				delete $temp_loci->[2]{$temp_eej};
			}
		}
	}

	return $temp_num_filtered_juncdist;
}

#Print the record in gtf format
sub print_gtf{
	my $temp_chr = $_[0];
	my $loci_ref = $_[1];
	my $temp_gtf_fh = $_[2];
	my $temp_enddist_fh = $_[3];
	my $temp_prefix = $_[4];
	my %temp_strand_flag = ('+', 0, '-', 1);
	foreach my $temp_strand (keys %{$loci_ref}){
		#print "Strand:$temp_strand\n";
		my $temp_loci_count = 1;
		LOCI:foreach my $temp_loci (@{$loci_ref->{$temp_strand}}){
			#print "Loci:@temp_loci\n";
			my %temp_eej = %{$temp_loci->[2]};
			unless(%temp_eej){
				next LOCI;
			}
			print $temp_gtf_fh "\n$temp_chr\tNIAP\tgene\t$temp_loci->[0]\t$temp_loci->[1]\t.\t$temp_strand\t.\tgene_id \"${temp_prefix}_${temp_chr}_$temp_strand_flag{$temp_strand}_${temp_loci_count}\"\;";
			my $temp_transcript_count = 1;
			foreach my $temp_eej (sort keys %temp_eej){
				#print "\tEEJ:$temp_eej\n\t@{$temp_eej{$temp_eej}}\n";
				TRANSCRIPT:foreach my $temp_transcript (@{$temp_eej{$temp_eej}}){
					#print "@temp_transcript\n";
					my $temp_eej_new = $temp_transcript->[0] . $temp_eej . $temp_transcript->[1];
					#print "$temp_eej_new\n";
					my @temp_eej = split(',', $temp_eej_new);
					my $temp_num_exon = @temp_eej;
					my @temp_polya_flag = keys %{$temp_transcript->[3]};
					my $temp_polya_flag = join(',', @temp_polya_flag);
					print $temp_gtf_fh "\n$temp_chr\tNIAP\ttranscript\t$temp_transcript->[0]\t$temp_transcript->[1]\t.\t$temp_strand\t.\tgene_id \"${temp_prefix}_${temp_chr}_$temp_strand_flag{$temp_strand}_${temp_loci_count}\"\; transcript_id \"${temp_prefix}_${temp_chr}_$temp_strand_flag{$temp_strand}_${temp_loci_count}.${temp_transcript_count}\"\; depth \"$temp_transcript->[2]\"\; polya_flag \"$temp_polya_flag\"\;";
					if($temp_enddist_fh){
						print $temp_enddist_fh "\n$temp_chr\t$temp_strand\t${temp_prefix}_${temp_chr}_$temp_strand_flag{$temp_strand}_${temp_loci_count}\t${temp_prefix}_${temp_chr}_$temp_strand_flag{$temp_strand}_${temp_loci_count}.${temp_transcript_count}";
						my @temp_end_col = (4, 5);
						if($temp_strand eq "-"){
							@temp_end_col = (5, 4);
						}
						foreach my $i (@temp_end_col){
							print $temp_enddist_fh "\t";
							foreach my $temp_pos (sort {$a <=> $b}keys %{$temp_transcript->[$i]}){
								print $temp_enddist_fh "$temp_pos:$temp_transcript->[$i]{$temp_pos},";
							}
						}
					}
					for(my $i = 0; $i < $temp_num_exon; $i++){
						my $temp_cur_exon_num;
						my @temp_exon = split('-', $temp_eej[$i]);
						#Handle exon number for record at the negative strand
						if($temp_strand eq "-"){
							$temp_cur_exon_num = $temp_num_exon - $i;
						}
						#Handle exon number for record at the positive and no strand
						else{
							$temp_cur_exon_num = $i + 1;
						}
						print $temp_gtf_fh "\n$temp_chr\tNIAP\texon\t$temp_exon[0]\t$temp_exon[1]\t.\t$temp_strand\t.\tgene_id \"${temp_prefix}_${temp_chr}_$temp_strand_flag{$temp_strand}_${temp_loci_count}\"\; transcript_id \"${temp_prefix}_${temp_chr}_$temp_strand_flag{$temp_strand}_${temp_loci_count}.${temp_transcript_count}\"\; exon_number \"$temp_cur_exon_num\"\; depth \"$temp_transcript->[2]\"\; polya_flag \"$temp_polya_flag\"\;";
					}
					$temp_transcript_count++;
				}
			}
			$temp_loci_count++;
		}
		#Clear the old record
		delete($loci_ref->{$temp_strand});
	}

	return 1;
}

#Print the nearest distance between all unique raw splice junctions and their nearest referece splice junctions.
sub print_juncdist{
	my $juncdist_ref = $_[0];
	my $temp_out = $_[1];
	open(TEMP, ">$temp_out") or die "Cannot create $temp_out!\n";
	print TEMP "Chr\tStrand\tQuery\tNearest\tDistance";
	foreach my $temp_chr (keys %{$juncdist_ref}){
		foreach my $temp_strand (keys %{$juncdist_ref->{$temp_chr}}){
			foreach my $temp_query (keys %{$juncdist_ref->{$temp_chr}{$temp_strand}}){
				foreach my $temp_nearest (keys %{$juncdist_ref->{$temp_chr}{$temp_strand}{$temp_query}}){
					print TEMP "\n$temp_chr\t$temp_strand\t$temp_query\t$temp_nearest\t$juncdist_ref->{$temp_chr}{$temp_strand}{$temp_query}{$temp_nearest}";
				}
			}
		}
	}
	close TEMP;

	return 1;
}


