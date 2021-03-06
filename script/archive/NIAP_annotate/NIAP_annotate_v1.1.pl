#!/bin/perl
#Author: ALin
#Purpose: To compare a query gtf with a refence gtf and annotate the transcripts in the query gtf.
#Usage: NIAP_annotate_v1.1.pl [option]
#		-query <String> Input query gtf file
#		-ref <String> Input reference gtf file
#		-o <String> Output file
#		-s <Boolean> Strand-specific query. (Default: Non-strand-specific)
#		-m <Boolean> Add those missing reference exons in the query gtf as in the refence gtf file
#		-btp <String> Path for bedtools
#		-h <Boolean> Help
#Change logs:
#	v1.1	This script is a succession of compare_Cufflinks_gtf_with_reference_gtf_v1.6.pl. It can now bettern handle cases of suspected isoform.

use strict;
use Getopt::Long;
use List::Util qw(min max sum);

my $ref = "";
my $query = "";
my $out = "";
my $stranded = 0;
my $missing = 0;
my $btp;
my $help = 0;



GetOptions(
	'query=s'	=>	\$query,
	'ref=s'	=>	\$ref,
        'o=s'   =>      \$out,
	's!'	=>	\$stranded,
	'm!'	=>	\$missing,
	'btp'	=>	\$btp,
	'h!'    =>      \$help,
);

unless($btp){
	$btp = "bedtools";
}

unless(`$btp`){
	print "Incorrect path for bedtools!\n";
	print_usage();
	exit;
}

if($help){
        print_usage();
        exit;
}

unless($query && $ref && $out){
	print_usage();
	exit;
}

open(QUERY, $query) or die "Cannot open $query!\n";
open(REF, $ref) or die "Cannot open $ref!\n";
open(OUT, ">temp_out.gtf") or die "Cannot create temp_out.gtf!\n";

my %ref_exon = ();
my %query_exon = ();
my %missing_exon = ();
my %annotated_exon = ();
my %ref_gene_id = ();

while(<REF>){
	chomp;
	my $line = $_;
	my @line = split('\t', $line);
	unless($line[2] eq "exon"){
		next;
	}
	$line[8] =~ /transcript_id \"([^\"]+)\"/;
	my $transcript_id = $1;
	$line[8] =~ /gene_id \"([^\"]+)\"/;
	my $gene_id = $1;
	unless(exists $ref_exon{$transcript_id}){
		$ref_exon{$transcript_id} = 0;
	}
	$ref_exon{$transcript_id}++;
	unless(exists $missing_exon{$transcript_id}){
		@{$missing_exon{$transcript_id}} = ();
	}	
	push(@{$missing_exon{$transcript_id}}, $line);
	unless(exists $ref_gene_id{$transcript_id}){
		$ref_gene_id{$transcript_id} = $gene_id;
	}
}
close REF;

while(<QUERY>){
	chomp;
	my $line = $_;
	my @line = split('\t', $line);
	unless($line[2] eq "exon"){
		next;
	}
	$line[8] =~ /transcript_id \"([^\"]+)\"/;
	my $transcript_id = $1;
	unless(exists $query_exon{$transcript_id}){
		$query_exon{$transcript_id} = 0;
	}
	$query_exon{$transcript_id}++;
}
close QUERY;

system("bedtools sort -i $query >  temp_query_sorted.gtf");
system("bedtools sort -i $ref >  temp_ref_sorted.gtf");

my %ref_strand = ();
my %query = ();
my %map = ();
my %query_strand = ();

open(IN, "bedtools intersect -wao -nonamecheck -a temp_query_sorted.gtf -b temp_ref_sorted.gtf |");

while(<IN>){
	chomp;
	my $line = $_;
	my @line = split('\t', $line);
	#Remove any "contained_in" transcript, which is initially generated by Cufflinks
        if($line[8] =~ / contained_in /){
                next;
        }
	unless($line[2] eq "exon"){
		next;
	}
	$line[8] =~ /transcript_id \"([^\;]+)\"\;/;
	my $query_transcript = $1;
	#Initiate the hash for the original query transcripts
	unless(exists $query{$query_transcript}){
		%{$query{$query_transcript}} = ();
	}
	#If non-strand-specific, force all query transcripts with single exon to be with no strand information.
	if(($stranded == 0) && ($query_exon{$query_transcript} == 1) && ($line[6] ne "\.")){
		$line[6] = "\.";
	}
	#Keep the original query gtf exon.
	my $temp_query = join("\t", @line[0, 1, 2, 3, 4, 5, 6, 7, 8]);
	$query{$query_transcript}{$temp_query} = 1;
	#Get the strand information of query transcript
	unless(exists $query_strand{$query_transcript}){
		$query_strand{$query_transcript} = $line[6];
	}
	unless(($line[11] eq "exon") || ($line[11] eq "\.")){
		next;
	}
	#Get reference transcript name
	#Amendment is needed for alternative reference annotation.
	my $ref_transcript = "";
	if($line[17] =~ /transcript_id \"([^\"]+)\"/){
		$ref_transcript = $1;
	}
	else{
		#No hit
		next;
	}
	#Check whether the query is really strand-specific
	if(($stranded == 1) && ($query_strand{$query_transcript} eq "\.")){
		print "$query_transcript is non-strand-specific!\n";
		#exit;
	}
	#Get the strand information of the reference transcript
	unless(exists $ref_strand{$ref_transcript}){
		$ref_strand{$ref_transcript} = $line[15];
	}
	#Check whether strands match
	if($stranded == 0){
		#Stranded disable
		if(($ref_strand{$ref_transcript} ne $query_strand{$query_transcript}) && ($query_strand{$query_transcript} ne "\.")){
			#Wrong strand
			next;
		}
	}
	else{
		#Stranded enable
		if($ref_strand{$ref_transcript} ne $query_strand{$query_transcript}){
			#Wrong strand
			next;
		}
	}
	#Mapping query transcript to reference
	unless(exists $map{$query_transcript}){
		%{$map{$query_transcript}} = ();
	}
	unless(exists $map{$query_transcript}{$ref_transcript}){
		@{$map{$query_transcript}{$ref_transcript}} = (0, 0);
	}
	#Get the exon number of the current reference exon
	$line[17] =~ /exon_number [\"]{0,1}(\d+)[\"]{0,1}/;
	my $current_ref_exon = $1;
	#First or last reference exon
	if($current_ref_exon == 1 || $current_ref_exon == $ref_exon{$ref_transcript}){
		if($ref_exon{$ref_transcript} == 1 && $query_exon{$query_transcript} == 1){
			#Both query and reference transcripts are single exon.
			$map{$query_transcript}{$ref_transcript}[0]++;
		}
		else{
			#Either query or reference transcript has more than one exon.
			if($line[3] == $line[12] || $line[4] == $line[13]){
				#Matching either one end of the exon
				$map{$query_transcript}{$ref_transcript}[0]++;
			}
		}
	}
	#Not first or last reference exon
	else{
		if(($line[3] == $line[12]) && ($line[4] == $line[13])){
			#Matching both exon ends
			$map{$query_transcript}{$ref_transcript}[0]++;
		}
	}
	$map{$query_transcript}{$ref_transcript}[1] += $line[18];
}
close IN;
system("rm temp_ref_sorted.gtf");
system("rm temp_query_sorted.gtf");

#Evaluate mapping

my $new_line_flag = 0;

foreach my $query_transcript (keys %map){
	my $anno_flag = 0;
	my %gene_hit = ();
	#Loop through the map and identify the best annotation
	foreach my $ref_transcript (keys %{$map{$query_transcript}}){	#There may be a bug in this loop such that query transcripts spanning more than one annoatated genes will be assigned to one of the the genes randomly, according to the order in the hash. #But this bug is kind of solved now. #A change of the way to get gene ID.
		if(($map{$query_transcript}{$ref_transcript}[0] == $query_exon{$query_transcript}) && ($map{$query_transcript}{$ref_transcript}[0] == $ref_exon{$ref_transcript})){
			$anno_flag = 1;
		}
		unless(exists $gene_hit{$ref_gene_id{$ref_transcript}}){
			@{$gene_hit{$ref_gene_id{$ref_transcript}}} = (0, 0);
		}
		if($gene_hit{$ref_gene_id{$ref_transcript}}[0] < $map{$query_transcript}{$ref_transcript}[0]){
			@{$gene_hit{$ref_gene_id{$ref_transcript}}} = @{$map{$query_transcript}{$ref_transcript}};
		}
		elsif($gene_hit{$ref_gene_id{$ref_transcript}}[0] == $map{$query_transcript}{$ref_transcript}[0]){
			if($gene_hit{$ref_gene_id{$ref_transcript}}[1] < $map{$query_transcript}{$ref_transcript}[1]){
				@{$gene_hit{$ref_gene_id{$ref_transcript}}} = @{$map{$query_transcript}{$ref_transcript}};
			}
		}
		#print "$query_transcript\t$ref_transcript\t@{$map{$query_transcript}{$ref_transcript}}\t$anno_flag\n";
	}
	my $ref_gene_anno = "";
	my $ref_transcript_anno = "";
	my $num_gene_hit = (keys %gene_hit);
	my $num_transcript_hit = (keys %{$map{$query_transcript}});
	my @sorted_gene = sort { ($gene_hit{$b}[0] <=> $gene_hit{$a}[0]) || ($gene_hit{$b}[1] <=> $gene_hit{$a}[1]) } keys %gene_hit;
	my @sorted_transcript = sort { ($map{$query_transcript}{$b}->[0] <=> $map{$query_transcript}{$a}->[0]) || ($map{$query_transcript}{$b}->[1] <=> $map{$query_transcript}{$a}->[1]) } keys %{$map{$query_transcript}};
	if($anno_flag == 0){
		if($num_gene_hit == 1){
			$ref_gene_anno = $sorted_gene[0];
		}
		else{
			if($query_exon{$query_transcript} == 1){
				$ref_gene_anno = $sorted_gene[0];
			}
			else{
				if($gene_hit{$sorted_gene[1]}[0] == 0){
					$ref_gene_anno = $sorted_gene[0];
				}
			}
		}
	}
	else{
		if($num_gene_hit == 1){
			$ref_gene_anno = $sorted_gene[0];
			$ref_transcript_anno = $sorted_transcript[0];
		}
		else{
			if($query_exon{$query_transcript} == 1){
				my $sum = 0;
				foreach my $gene_hit (%gene_hit){
					$sum += $gene_hit{$gene_hit}[0];
				}
				$sum /= $num_gene_hit;
				my $ref_sum = 1 / $num_gene_hit;
				if($sum == $ref_sum){
					$ref_gene_anno = $sorted_gene[0];
					$ref_transcript_anno = $sorted_transcript[0];
				}
				else{
					$anno_flag = 0;
				}
			}
			else{
				$ref_gene_anno = $sorted_gene[0];
				$ref_transcript_anno = $sorted_transcript[0];
			}
		}
	}
	#Loop through the original query entry and print.
	foreach my $line (keys %{$query{$query_transcript}}){
		if($ref_gene_anno ne ""){
			$line =~ /gene_id \"([^\;]+)\"\;/;
			$line =~ s/$1/$ref_gene_anno/;
		}
		if($anno_flag == 1){
			#Perfect match
			$line =~ /transcript_id \"([^\;]+)\"\;/;
			my $oId = $1;
			unless(exists $annotated_exon{$ref_transcript_anno}){
				$annotated_exon{$ref_transcript_anno} = $oId;
			}
			next;
		}
		elsif($ref_gene_anno eq ""){
			#Strange transcript, e.g. a transcript that connects two adjunction transcripts.
			$line .= "status \"suspected_isoform\"\;";
		}
		else{
			#New isoform
			$line .= " status \"novel_isoform\"\;";
		}
		if($new_line_flag == 0){
			$new_line_flag = 1;
		}
		else{
			print OUT "\n";
		}
		print OUT $line;
	}
	delete $query{$query_transcript};
}

%map = ();

#Annotated transcript
foreach my $anno_ref_id (keys %annotated_exon){
	foreach my $line (@{$missing_exon{$anno_ref_id}}){
		$line .= " oId \"$annotated_exon{$anno_ref_id}\"\; status \"annotated_isoform\"\;";
		print OUT "\n$line";
	}
	if($missing){
		delete $missing_exon{$anno_ref_id};
	}
}

%annotated_exon = ();

#Query transcript left
foreach my $query_transcript (keys %query){
	foreach my $line (keys %{$query{$query_transcript}}){
		$line .= " status \"unknown\"\;";
		print OUT "\n$line";
	}
}

#Print missing exons in the reference
if($missing){
	foreach my $transcript_id (keys %missing_exon){
		foreach my $exon (@{$missing_exon{$transcript_id}}){
			print OUT "\n$exon status \"missing_isoform\"\;";
		}
	}
}
close OUT;
system("bedtools sort -i temp_out.gtf > $out");
system("rm temp_out.gtf");


sub print_usage{
	print "Usage: perl NIAP_annotate_v1.1.pl [option]\n\t-query <String> Input query gtf file\n\t-ref <String> Input reference gtf file\n\t-o <String> Output file\n\t-s <Boolean> Strand-specific query. (Default: Non-strand-specific)\n\t-m <Boolean> Add those missing reference exons in the query as in the reference gtf file\n\t-h <Boolean> Help\n";
}





