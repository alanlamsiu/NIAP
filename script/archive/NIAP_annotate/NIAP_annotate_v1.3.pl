#!/bin/perl
#Author: ALin
#Purpose: To compare a query gtf with a refence gtf and annotate the transcripts in the query gtf.
#Usage: NIAP_annotate_v1.3.pl [option]
#		-query <String> Input query gtf file
#		-ref <String> Input reference gtf file
#		-o <String> Output file
#		-s <Boolean> Strand-specific query. (Default: Non-strand-specific)
#		-m <Boolean> Add those missing reference exons in the query gtf as in the refence gtf file
#		-btp <String> Path for bedtools
#		-h <Boolean> Help
#Change logs:
#	v1.1	2021-02	This script is a succession of compare_Cufflinks_gtf_with_reference_gtf_v1.6.pl. It can now bettern handle cases of suspected isoform.
#	v1.2	2021-03	This script is now able to make detailed annotation of the unkown transcripts, as antisense, intergenic, intronic and overlapping, partial transcripts of the annotated ones, which may be truncated, novel transcripts, including those with alternative ends and novel structures.
#	v1.3	2021-03	It now considers the direction of the cases with novel ends.

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

system("bedtools sort -i $query >  temp_query_sorted.gtf");
system("bedtools sort -i $ref >  temp_ref_sorted.gtf");
open(QUERY, "temp_query_sorted.gtf") or die "Cannot open $query!\n";
open(REF, "temp_ref_sorted.gtf") or die "Cannot open $ref!\n";
open(OUT, ">temp_out.gtf") or die "Cannot create temp_out.gtf!\n";

my %ref_transcript = ();
my %query_transcript = ();
my %ref_gene = ();
my %query_gene = ();

#Load reference annotation
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
	unless(exists $ref_transcript{$transcript_id}){
		%{$ref_transcript{$transcript_id}} = ();
		$ref_transcript{$transcript_id}{'gene'} = $gene_id;
		$ref_transcript{$transcript_id}{'chr'} = $line[0];
		$ref_transcript{$transcript_id}{'source'} = $line[1];
		$ref_transcript{$transcript_id}{'left'} = $line[3];
		$ref_transcript{$transcript_id}{'right'} = $line[4];
		$ref_transcript{$transcript_id}{'strand'} = $line[6];
		@{$ref_transcript{$transcript_id}{'exon'}} = ();
		@{$ref_transcript{$transcript_id}{'original'}} = ();
	}
	my @exon = ($line[3], $line[4]);
	push(@{$ref_transcript{$transcript_id}{'exon'}}, \@exon);
	push(@{$ref_transcript{$transcript_id}{'original'}}, $line);
	unless(exists $ref_gene{$gene_id}){
		%{$ref_gene{$gene_id}} = ();
		$ref_gene{$gene_id}{'chr'} = $line[0];
		$ref_gene{$gene_id}{'source'} = $line[1];
		$ref_gene{$gene_id}{'left'} = $line[3];
		$ref_gene{$gene_id}{'right'} = $line[4];
		$ref_gene{$gene_id}{'strand'} = $line[6];
	}
	if(($ref_gene{$gene_id}{'strand'} eq "\.") && ($line[6] ne "\.")){
		$ref_gene{$gene_id}{'strand'} = $line[6];
	}
	if($ref_gene{$gene_id}{'left'} > $line[3]){
		$ref_gene{$gene_id}{'left'} = $line[3];
	}
	if($ref_gene{$gene_id}{'right'} < $line[4]){
		$ref_gene{$gene_id}{'right'} = $line[4];
	}
}
close REF;

#Create gene annotation for the reference
my $new_line_flag = 0;

open(RGENE, ">temp_ref_gene.gtf") or die "Cannot create temp_ref_gene.gtf!\n";

foreach my $gene_id (keys %ref_gene){
	if($new_line_flag == 0){
		$new_line_flag = 1;
	}
	else{
		print RGENE "\n";
	}
	print RGENE "$ref_gene{$gene_id}{'chr'}\t$ref_gene{$gene_id}{'source'}\tgene\t$ref_gene{$gene_id}{'left'}\t$ref_gene{$gene_id}{'right'}\t.\t$ref_gene{$gene_id}{'strand'}\t.\tgene_id \"$gene_id\"\;";
}
close RGENE;

#Load query annotation
while(<QUERY>){
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
	unless(exists $query_transcript{$transcript_id}){
		%{$query_transcript{$transcript_id}} = ();
		$query_transcript{$transcript_id}{'gene'} = $gene_id;
		$query_transcript{$transcript_id}{'chr'} = $line[0];
		$query_transcript{$transcript_id}{'source'} = $line[1];
		$query_transcript{$transcript_id}{'left'} = $line[3];
		$query_transcript{$transcript_id}{'right'} = $line[4];
		$query_transcript{$transcript_id}{'strand'} = $line[6];
		@{$query_transcript{$transcript_id}{'exon'}} = ();
		@{$query_transcript{$transcript_id}{'original'}} = ();
	}
	my @exon = ($line[3], $line[4]);
	push(@{$query_transcript{$transcript_id}{'exon'}}, \@exon);
        push(@{$query_transcript{$transcript_id}{'original'}}, $line);
	unless(exists $query_gene{$gene_id}){
		%{$query_gene{$gene_id}} = ();
		$query_gene{$gene_id}{'chr'} = $line[0];
		$query_gene{$gene_id}{'source'} = $line[1];
		$query_gene{$gene_id}{'left'} = $line[3];
		$query_gene{$gene_id}{'right'} = $line[4];
		$query_gene{$gene_id}{'strand'} = $line[6];
	}
	if(($query_gene{$gene_id}{'strand'} eq "\.") && ($line[6] ne "\.")){
		$query_gene{$gene_id}{'strand'} = $line[6];
	}
	if($query_gene{$gene_id}{'left'} > $line[3]){
		$query_gene{$gene_id}{'left'} = $line[3];
	}
	if($query_gene{$gene_id}{'right'} < $line[4]){
		$query_gene{$gene_id}{'right'} = $line[4];
	}
}
close QUERY;

#Create gene annotation for the query
$new_line_flag = 0;

open(QGENE, ">temp_query_gene.gtf") or die "Cannot create temp_query_gene.gtf!\n";

foreach my $gene_id (keys %query_gene){
	if($new_line_flag == 0){
		$new_line_flag = 1;
	}
	else{
		print QGENE "\n";
	}
	print QGENE "$query_gene{$gene_id}{'chr'}\t$query_gene{$gene_id}{'source'}\tgene\t$query_gene{$gene_id}{'left'}\t$query_gene{$gene_id}{'right'}\t.\t$query_gene{$gene_id}{'strand'}\t.\tgene_id \"$gene_id\"\;";
}
close QGENE;

system("bedtools sort -i temp_ref_gene.gtf > temp_ref_gene_sorted.gtf");
system("bedtools sort -i temp_query_gene.gtf > temp_query_gene_sorted.gtf");

#Classify unannotated genes in the query
my %class = ();

open(TEMP, "bedtools intersect -wa -S -nonamecheck -a temp_query_gene_sorted.gtf -b temp_ref_gene_sorted.gtf |");

while(<TEMP>){
	chomp;
	my $line = $_;
	#print "$line\n";
	my @line = split('\t', $line);
	$line[8] =~ /gene_id \"([^\;]+)\"\;/;
	my $gene_id = $1;
	#print "$gene_id\n";
	unless(exists $class{$gene_id}){
		$class{$gene_id} = "antisense";
	}
}
close TEMP;

open(TEMP, "bedtools intersect -v -wa -nonamecheck -a temp_query_gene_sorted.gtf -b temp_ref_gene_sorted.gtf |");

while(<TEMP>){
	chomp;
	my $line = $_;
	my @line = split('\t', $line);
	$line[8] =~ /gene_id \"([^\;]+)\"\;/;
	my $gene_id = $1;
	unless(exists $class{$gene_id}){
		$class{$gene_id} = "intergenic";
	}
}
close TEMP;

open(TEMP, "bedtools intersect -v -wa -nonamecheck -a temp_query_gene_sorted.gtf -b temp_ref_sorted.gtf |");

while(<TEMP>){
	chomp;
	my $line = $_;
	my @line = split('\t', $line);
	$line[8] =~ /gene_id \"([^\;]+)\"\;/;
	my $gene_id = $1;
	unless(exists $class{$gene_id}){
		$class{$gene_id} = "intronic";
	}
}
close TEMP;

system("rm temp_ref_gene.gtf");
system("rm temp_ref_gene_sorted.gtf");
system("rm temp_query_gene.gtf");
system("rm temp_query_gene_sorted.gtf");

#Map query transcirpts to reference
my %map = ();

open(IN, "bedtools intersect -wao -nonamecheck -a temp_query_sorted.gtf -b temp_ref_sorted.gtf |");

while(<IN>){
	chomp;
	my $line = $_;
	my @line = split('\t', $line);
	#Remove any "contained_in" transcript, which is initially generated by Cufflinks
        if($line[8] =~ / contained_in /){
                next;
        }
	if(($line[2] ne "exon") || ($line[11] ne "exon")){
		next;
	}
	
	$line[8] =~ /transcript_id \"([^\;]+)\"\;/;
	my $query_transcript_id = $1;
	my $query_exon_num = @{$query_transcript{$query_transcript_id}{'exon'}};
	#If non-strand-specific, force all query transcripts with single exon to be with no strand information.
	if(($stranded == 0) && ($query_exon_num == 1) && ($line[6] ne "\.")){
		$line[6] = "\.";
		$query_transcript{$query_transcript_id}{'strand'} = "\.";
	}
	#Get reference transcript name
	#Amendment is needed for alternative reference annotation.
	$line[17] =~ /transcript_id \"([^\"]+)\"/;
	my $ref_transcript_id = $1;
	my $ref_exon_num = @{$ref_transcript{$ref_transcript_id}{'exon'}};
	#Check whether the query is really strand-specific
	if(($stranded == 1) && ($query_transcript{$query_transcript_id}{'strand'} eq "\.")){
		print "$query_transcript_id is non-strand-specific!\n";
		#exit;
	}
	#Check whether strands match
	if($stranded == 0){
		#Stranded disable
		if(($line[6] ne $line[15]) && ($line[6] ne "\.")){
			#Wrong strand
			next;
		}
	}
	else{
		#Stranded enable
		if($line[6] ne $line[15]){
			#Wrong strand
			next;
		}
	}
	#Mapping query transcript to reference
	unless(exists $map{$query_transcript_id}){
		%{$map{$query_transcript_id}} = ();
	}
	unless(exists $map{$query_transcript_id}{$ref_transcript_id}){
		%{$map{$query_transcript_id}{$ref_transcript_id}} = ();
		@{$map{$query_transcript_id}{$ref_transcript_id}{'query'}} = (0) x $query_exon_num;
		@{$map{$query_transcript_id}{$ref_transcript_id}{'ref'}} = (0) x $ref_exon_num;
		$map{$query_transcript_id}{$ref_transcript_id}{'match'} = 0;
		$map{$query_transcript_id}{$ref_transcript_id}{'partial'} = 0;
		$map{$query_transcript_id}{$ref_transcript_id}{'query_mis'} = 0;
		$map{$query_transcript_id}{$ref_transcript_id}{'ref_mis'} = 0;
		$map{$query_transcript_id}{$ref_transcript_id}{'overlap'} = 0;
		$map{$query_transcript_id}{$ref_transcript_id}{'annotate'} = 0;
	}
	#Get the exon number of the current query exon
	my $cur_query_exon_num;
	SEARCH:for(my $i = 0; $i < $query_exon_num; $i++){
		my $hit_flag = 0;
		if($query_exon_num == 1){
			$hit_flag = 1;
		}
		elsif($i == 0){
			if($query_transcript{$query_transcript_id}{'exon'}[$i][1] == $line[4]){
				$hit_flag = 1;
			}
		}
		elsif($i == ($query_exon_num - 1)){
			if($query_transcript{$query_transcript_id}{'exon'}[$i][0] == $line[3]){
				$hit_flag = 1;
			}
		}
		elsif(($query_transcript{$query_transcript_id}{'exon'}[$i][0] == $line[3]) && ($query_transcript{$query_transcript_id}{'exon'}[$i][1] == $line[4])){
			$hit_flag = 1;
		}
		if($hit_flag == 1){
			if($line[6] eq "-"){
				$cur_query_exon_num = $query_exon_num - $i;
			}
			else{
				$cur_query_exon_num = $i + 1;
			}
			last SEARCH;
		}
	}
	#Get the exon number of the current reference exon
	my $cur_ref_exon_num;
	SEARCH:for(my $i = 0; $i < $ref_exon_num; $i++){
		my $hit_flag = 0;
		if($ref_exon_num == 1){
			$hit_flag = 1;
		}
		elsif($i == 0){
			if($ref_transcript{$ref_transcript_id}{'exon'}[$i][1] == $line[13]){
				$hit_flag = 1;
			}
		}
		elsif($i == ($ref_exon_num - 1)){
			if($ref_transcript{$ref_transcript_id}{'exon'}[$i][0] == $line[12]){
				$hit_flag = 1;
			}
		}
		elsif(($ref_transcript{$ref_transcript_id}{'exon'}[$i][0] == $line[12]) && ($ref_transcript{$ref_transcript_id}{'exon'}[$i][1] == $line[13])){
			$hit_flag = 1;
		}
		if($hit_flag == 1){
			if($line[15] eq "-"){
				$cur_ref_exon_num = $ref_exon_num - $i;
			}
			else{
				$cur_ref_exon_num = $i + 1;
			}
			last SEARCH;
		}
	}
	my $match = 0;
	if(($ref_exon_num == 1) && ($query_exon_num == 1)){
		#Both query and reference transcripts are monoexonic.
		$match = 1;
	}
	elsif(($ref_exon_num > 1) && ($query_exon_num > 1)){
		#Both query and reference transcripts are multiexonic.
		if(($cur_query_exon_num == 1) || ($cur_ref_exon_num == 1)){
			if($line[6] eq "\-"){
				#The rightmost case
				if($line[3] == $line[12]){
					$match = 1;
				}
			}
			else{
				#The leftmost case
				if($line[4] == $line[13]){
					$match = 1;
				}
			}
		}
		elsif(($cur_query_exon_num == $query_exon_num) || ($cur_ref_exon_num == $ref_exon_num)){
			if($line[6] eq "\-"){
				#The leftmost case
				if($line[4] == $line[13]){
					$match = 1;
				}
			}
			else{
				#The rightmost case
				if($line[3] == $line[12]){
					$match = 1;
				}
			}
		}
		else{
			#Internal exon
			if(($line[3] == $line[12]) && ($line[4] == $line[13])){
				#Matching both exon ends
				$match = 1;
			}
		}
	}
	if($match == 1){
		if($line[6] eq "\-"){
			$map{$query_transcript_id}{$ref_transcript_id}{'query'}[($query_exon_num - $cur_query_exon_num)] = 1;
			$map{$query_transcript_id}{$ref_transcript_id}{'ref'}[($ref_exon_num - $cur_ref_exon_num)] = 1;
		}
		else{
			$map{$query_transcript_id}{$ref_transcript_id}{'query'}[($cur_query_exon_num - 1)] = 1;
			$map{$query_transcript_id}{$ref_transcript_id}{'ref'}[($cur_ref_exon_num - 1)] = 1;
		}
		$map{$query_transcript_id}{$ref_transcript_id}{'match'}++;
	}
	$map{$query_transcript_id}{$ref_transcript_id}{'overlap'} += $line[18];
}
close IN;
system("rm temp_ref_sorted.gtf");
system("rm temp_query_sorted.gtf");

#Evaluate mapping
$new_line_flag = 0;

foreach my $query_transcript_id (keys %map){
	my $anno_flag = 0;
	my %gene_hit = ();
	my $query_exon_num = @{$query_transcript{$query_transcript_id}{'exon'}};
	#Loop through the map and identify the best annotation
	foreach my $ref_transcript_id (keys %{$map{$query_transcript_id}}){	#There may be a bug in this loop such that query transcripts spanning more than one annoatated genes will be assigned to one of the the genes randomly, according to the order in the hash. #But this bug is kind of solved now. #A change of the way to get gene ID.
		my $ref_gene_id = $ref_transcript{$ref_transcript_id}{'gene'};
		my $query_code = join('', @{$map{$query_transcript_id}{$ref_transcript_id}{'query'}});
		my $ref_code = join('', @{$map{$query_transcript_id}{$ref_transcript_id}{'ref'}});
		$map{$query_transcript_id}{$ref_transcript_id}{'query_mis'} = () = $query_code =~ /0/;
		$map{$query_transcript_id}{$ref_transcript_id}{'ref_mis'} = () = $ref_code =~ /0/;
		if(($query_code eq $ref_code) && !($query_code =~ /0/) && !($ref_code =~ /0/)){
			$map{$query_transcript_id}{$ref_transcript_id}{'annotate'} = 1;
			$anno_flag = 1;
		}
		elsif($ref_code =~ /$query_code/){
			$map{$query_transcript_id}{$ref_transcript_id}{'partial'} = 1;
		}
		#print "$query_transcript_id\t$query_code\t$ref_transcript_id\t$ref_code\t$map{$query_transcript_id}{$ref_transcript_id}{'match'}\t$map{$query_transcript_id}{$ref_transcript_id}{'partial'}\t$map{$query_transcript_id}{$ref_transcript_id}{'query_mis'}\t$map{$query_transcript_id}{$ref_transcript_id}{'ref_mis'}\t$map{$query_transcript_id}{$ref_transcript_id}{'overlap'}\t$map{$query_transcript_id}{$ref_transcript_id}{'annotate'}\n";
		unless(exists $gene_hit{$ref_gene_id}){
			%{$gene_hit{$ref_gene_id}} = ();
			$gene_hit{$ref_gene_id}{'match'} = 0;
			$gene_hit{$ref_gene_id}{'overlap'} = 0;
		}
		if($gene_hit{$ref_gene_id}{'match'} < $map{$query_transcript_id}{$ref_transcript_id}{'match'}){
			$gene_hit{$ref_gene_id}{'match'} = $map{$query_transcript_id}{$ref_transcript_id}{'match'};
			$gene_hit{$ref_gene_id}{'overlap'} = $map{$query_transcript_id}{$ref_transcript_id}{'overlap'};
		}
		elsif($gene_hit{$ref_gene_id}{'match'} == $map{$query_transcript_id}{$ref_transcript_id}{'match'}){
			if($gene_hit{$ref_gene_id}{'overlap'} < $map{$query_transcript_id}{$ref_transcript_id}{'overlap'}){
				$gene_hit{$ref_gene_id}{'match'} = $map{$query_transcript_id}{$ref_transcript_id}{'match'};
				$gene_hit{$ref_gene_id}{'overlap'} = $map{$query_transcript_id}{$ref_transcript_id}{'overlap'};
			}
		}
	}
	my $ref_gene_anno = "";
	my $ref_transcript_anno = "";
	my $num_gene_hit = (keys %gene_hit);
	#print "Number of genes hit: $num_gene_hit\n";
	my $num_transcript_hit = (keys %{$map{$query_transcript_id}});
	my @sorted_gene = sort { ($gene_hit{$b}{'match'} <=> $gene_hit{$a}{'match'}) || ($gene_hit{$b}{'overlap'} <=> $gene_hit{$a}{'overlap'}) } keys %gene_hit;
	my @sorted_transcript = sort { ($map{$query_transcript_id}{$b}{'annotate'} <=> $map{$query_transcript_id}{$a}{'annotate'}) || ($map{$query_transcript_id}{$b}{'match'} <=> $map{$query_transcript_id}{$a}{'match'}) || ($map{$query_transcript_id}{$b}{'partial'} <=> $map{$query_transcript_id}{$a}{'partial'}) || ($map{$query_transcript_id}{$a}{'query_mis'} <=> $map{$query_transcript_id}{$b}{'query_mis'}) || ($map{$query_transcript_id}{$a}{'ref_mis'} <=> $map{$query_transcript_id}{$b}{'ref_mis'}) || ($map{$query_transcript_id}{$b}{'overlap'} <=> $map{$query_transcript_id}{$a}{'overlap'}) || ($a cmp $b) } keys %{$map{$query_transcript_id}};
	#print "@sorted_transcript\n";
	if($anno_flag == 0){
		if($num_gene_hit == 1){
			$ref_gene_anno = $sorted_gene[0];
		}
		else{
			if($query_exon_num == 1){
				$ref_gene_anno = $sorted_gene[0];
			}
			else{
				if($gene_hit{$sorted_gene[1]}{'match'} == 0){
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
			if($query_exon_num == 1){
				my $sum = 0;
				foreach my $gene_hit (%gene_hit){
					$sum += $gene_hit{$gene_hit}{'match'};
				}
				$sum /= $num_gene_hit;
				if($sum == 1){
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
	#print "Annotation:$anno_flag\t$ref_gene_anno\t$ref_transcript_anno\n";
	if($anno_flag == 1){
		#Perfect match
		foreach my $line (@{$query_transcript{$query_transcript_id}{'original'}}){
			$line =~ /gene_id \"([^\;]+)\"\;/;
			$line =~ s/$1/$ref_gene_anno/;
			$line =~ /transcript_id \"([^\;]+)\"\;/;
			my $asm_id = $1;
			$line =~ s/$asm_id/$ref_transcript_anno/;
			if($new_line_flag == 0){
				$new_line_flag = 1;
			}
			else{
				print OUT "\n";
			}
			print OUT "$line asm_id \"$asm_id\"\; status \"annotated\"\;";
        }
		delete $query_transcript{$query_transcript_id};
		next;
	}
	my $status;
	if($ref_gene_anno){
		my $query_code = join('', @{$map{$query_transcript_id}{$sorted_transcript[0]}{'query'}});
		my $ref_code = join('', @{$map{$query_transcript_id}{$sorted_transcript[0]}{'ref'}});
		if($query_code =~ /0/){
			#Novel isoform with novel transcript structure
			$status = "novel_structure";
		}
		else{
			if($ref_code =~ /^0*1+0*$/){
				#print "Could be partial...\n";
				my $novel_5 = 0;
				my $novel_3 = 0;
				my $first_seg_start = 0;
				if($ref_code =~ /(^0+)1+0*$/){
					my $first_seg = $1;
					$first_seg_start = length($first_seg);
				}
				#print "First seg start: $first_seg_start\n";
				if($query_transcript{$query_transcript_id}{'exon'}[0][0] < $ref_transcript{$sorted_transcript[0]}{'exon'}[$first_seg_start][0]){
					if($query_transcript{$query_transcript_id}{'strand'} eq "-"){
						$novel_3 = 1;
					}
					else{
						$novel_5 = 1;
					}
				}
				my $last_seg_start = (length($ref_code) - 1);
				if($ref_code =~ /(^0*1+)0*$/){
					my $last_seg = $1;
					$last_seg_start = (length($last_seg) - 1);
				}
				#print "Last seg start: $last_seg_start\n";
				if($query_transcript{$query_transcript_id}{'exon'}[-1][1] > $ref_transcript{$sorted_transcript[0]}{'exon'}[$last_seg_start][1]){
					if($query_transcript{$query_transcript_id}{'strand'} eq "-"){
						$novel_5 = 1;
					}
					else{
						$novel_3 = 1;
					}
				}
				if(($novel_5 == 1) && ($novel_3 == 1)){
					#Novel isoforms with varied transcript ends
					$status = "novel_ends";
				}
				elsif(($novel_5 == 1) && ($novel_3 == 0)){
					#Novel isoform with varied 5' end
					$status = "novel_5'";
				}
				elsif(($novel_5 == 0) && ($novel_3 == 1)){
					#Novel isoform with varied 3' end
					$status = "novel_3'";
				}
				else{
					#Partial transcript
					$status = "partial";
				}
			}
			else{
				#Novel isoform with novel transcript structure
				$status = "novel_structure";
			}
		}
	}
	else{
		#Strange transcript, e.g. a transcript that connects multiple adjunction transcripts.
		$status = "suspected";
	}
	#Loop through the original query entry and print.
	foreach my $line (@{$query_transcript{$query_transcript_id}{'original'}}){
		if($ref_gene_anno){
			$line =~ /gene_id \"([^\;]+)\"\;/;
			$line =~ s/$1/$ref_gene_anno/;
		}
		if($new_line_flag == 0){
			$new_line_flag = 1;
		}
		else{
			print OUT "\n";
		}
		print OUT "$line status \"$status\"\;";
	}
	delete $query_transcript{$query_transcript_id};
}

%map = ();

#Query transcript left
foreach my $query_transcript_id (keys %query_transcript){
	foreach my $line (@{$query_transcript{$query_transcript_id}{'original'}}){
		#print "$line\n";
		if(exists $class{$query_transcript{$query_transcript_id}{'gene'}}){
			$line .= " status \"$class{$query_transcript{$query_transcript_id}{'gene'}}\"\;";
		}
		else{
			$line .= " status \"overlapping\"\;";
		}
		if($new_line_flag == 0){
			$new_line_flag = 1;
		}
		else{
			print OUT "\n";
		}
		print OUT "$line";
	}
}

#Print missing exons in the reference
if($missing){
	foreach my $ref_transcript_id (keys %ref_transcript){
		foreach my $line (@{$ref_transcript{$ref_transcript_id}{'original'}}){
			if($new_line_flag == 0){
				$new_line_flag = 1;
			}
			else{
				print OUT "\n";
			}
			print OUT "$line status \"missing\"\;";
		}
	}
}
close OUT;
system("bedtools sort -i temp_out.gtf > $out");
system("rm temp_out.gtf");


sub print_usage{
	print "Usage: perl NIAP_annotate_v1.3.pl [option]\n\t-query <String> Input query gtf file\n\t-ref <String> Input reference gtf file\n\t-o <String> Output file\n\t-s <Boolean> Strand-specific query. (Default: Non-strand-specific)\n\t-m <Boolean> Add those missing reference exons in the query as in the reference gtf file\n\t-h <Boolean> Help\n";
}





