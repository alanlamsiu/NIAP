# NIAP
NIAP (Naive Isoform Analysis Pipeline) is a gernalized pipeline for analyzing short-read and long-read RNA-seq data. It is naive because it doesn't drop any transcript regarding the structure and abundance. Instead, users can perform any filtering step with the output from NIAP. There are mainly three tasks performed by NIAP, including (1) transcript merging, (2) transcript annotation, (3) transcript comparison and (4) transcript consolidation for long-read RNA-seq data.

## Installation

`git clone https://github.com/alanlamsiu/NIAP.git`

Scripts performing the four tasks mentioned above can be found in ./NIAP/script.

## Computational environment

All scripts have been tested in Perl (v5.22.1) and Ubuntu.

bedtools (https://bedtools.readthedocs.io/en/latest/) is requierd for transcript merging, transcript annotation and transcript comparison. The path to bedtools can either be put in the $PATH or speecified when running idividual scripts.

samtools (https://github.com/samtools/samtools) is requierd. The path to bedtools can either be put in the $PATH or speecified when running idividual scripts.

## 1. Transcript merging

In the transcript merging step, the script `NIAP_merge_v1.1.pl` is used to merge transcripts with the same structure, either in one or multiple .gtf file(s). Multiexonic transcirpts with the exact set of exon-exon junction(s) and being transcribed from the same direction will be merged into one representative transcript. On the other hand, overlapping monoexonic transcripts will be merged to generate a representative transcript. In both cases, the longest 5' and 3' ends, among the transcripts to be merged, will be taken for the representative transcript. A schematic illustration is shown below.

![NIAP_merge](https://user-images.githubusercontent.com/34832128/112115696-70761500-8bf4-11eb-8191-70532280b44b.jpg)

### Usage

```
perl NIAP_merge_v1.1.pl [option]
  -in <String> Input list of .gtf file(s)
  -out <String> Output base
  -exp <String> The atribute for the getting expression value (cov, FPKM, TPM or etc.)
  -prefix <String> The prefix of gene and transcript IDs (Default: NIAP)
  -btp <String> The path to bedtools
  -h <Boolean> Help
```

The input is a plain text file specifying the .gtf file(s) to be merged, one for each row. The full path(s) to the file(s) is(are) needed. There can be a single file or multiple files.

A .gtf file of the representative transcirpts after merging and a plain text file, with a sufix of `"_id.txt"`, indicating the source of each representative transcripts, will be generated. The `"_id.txt"` file is consisted of three columns, where the first two columns indicating the gene ID and transcript ID of the representative transcripts, the third one indicating `","` delineated original IDs before merging.

If the `-exp` option is specified for extracting the expression levele, or any other quantitative properties, of the transcripts, a summary table with the suffix `_exp.txt` will be generated.

The IDs of a gene and a transcript after merging as `"<PREFIX>_<CHR>_[01]_<LOCI_NUM>"` and `"<PREFIX>_<CHR>_[01]_<LOCUS_NUM>.<TRANSCRIPT_NUM>"`, respectively. The `"<PREFIX>"` value is `"NIAP"` by default or specified by users and can be specified by the `-prefix` option. The `"<CHR>"` value indicates the chromosome name. The `"[01]"` value indicates the strand of the locus, or transcript, with `"0"` meaning the `"+"` strand and `"1"` meaning the `"-"` one. The `"<LOCI_NUM>"` indicates the positional order of a locus at a strand. A locus is defined as a region consisted of overlapping transcripts, with each overlapping region containing at least two transcripts. The `"<TRANSCRIPT_NUM>"` is a unique identifier of a transcript in the corresponding locus.

The `-btp` is not necessary if bedtools is already in $PATH.
  
## 2. Transcript annotation

In the transcript annotation step, the `NIAP_annotate_v1.3.pl` script is used to compare a query .gtf file with a reference .gtf file. The purpoes is similar to what can be done by GffCompare (https://ccb.jhu.edu/software/stringtie/gffcompare.shtml), but the classification is a bit different in NIAP as shown by the schematic illustration below.

![NIAP_annotate](https://user-images.githubusercontent.com/34832128/112925260-29cc7180-9144-11eb-9894-acca4d747196.jpg)

### Usage

```
perl NIAP_annotate_v1.3.pl [option]
  -query <String> Input query gtf file
  -ref <String> Input reference gtf file
  -o <String> Output file
  -s <Boolean> Strand-specific query. (Default: Non-strand-specific)
  -m <Boolean> Add those missing reference exons in the query as in the reference gtf file
  -btp <String> Path for bedtools
  -h <Boolean> Help
```

Both query and reference inputs should be in proper .gtf format.

The output is a .gtf file with the `"asm_id"` and `"status"` attributes. The `"asm_id"` can be only found in annotated transcripts indicating the transcript ID before the annotation step. The value of `"status"` can be referred to the classification shown above.

If the query transcripts are strand-specific, meaning there is no transcript with the strand column shown as a `"."`, the `-s` option should be specified, which will make a difference in handling monoexonic transcripts while leaving those multiexonic unaffected.

If users wish to include transcripts only present in the reference, whose `"status"` attribute will be shown as `"missing"`, the `-m` option can be specified.

The `-btp` is not necessary if bedtools is already in $PATH.

## 3. Transcript comparison

In the transcript comparison step, similar to transcript annotation, a query .gtf file is compared with a reference file, but in .gff format. The script `NIAP_compare_v1.1.pl` is used for this step. Unlike in transcript annotation, where detailed classification is performed, transcrpit comparison simply tells whether query transcripts are supported by reference transcripts. Transcript comparison is usually useful when there are transcripts (in .fa format) from external sources. It is recommended to aligned these transcripts to the reference genome first, using GMAP (http://research-pub.gene.com/gmap/), to generate the reference .gff file.

### USage

```
perl NIAP_compare_v1.1.pl [option]
  -query <String> Input .gtf file to be checked
  -ref <String> Input .gff file of external transcripts
  -o <String> Output file
  -s <Boolean> Strand-specific query. (Default: Non-strand-specific)
  -btp <String> Path for bedtools
  -h <Boolean> Help
```

The query file should be in .gtf format, while the reference file should be in .gff format generated by GMAP.

The output is a plain text file consisted of three columns, where the first two columns indicating the gene ID and transcript ID of the query transcripts, the third one indicating `","` delineated IDs of the external transcripts.

If the query transcripts are strand-specific, meaning there is no transcript with the strand column shown as a `"."`, the `-s` option should be specified, which will make a difference in handling monoexonic transcripts while leaving those multiexonic unaffected.

The `-btp` is not necessary if bedtools is already in $PATH.

## 4. Transcript consolidation for long-read RNA-seq data

The transcript consolidation step is mainly used for analyzing long-read RNA-seq data, including Nanopore direct RNA-seq (dRNA-seq) and Pac-Bio Iso-Seq (Iso-Seq) data. The `NIAP_v1.2.pl` script is used for this purpose. The overall workflow is shown in the schematic illustration below, which can be divided into Iso-Seq and dRNA-seq sessions.

![NIAP_consolidate](https://user-images.githubusercontent.com/34832128/113243647-f88ba700-92e5-11eb-8460-acd8991fcd94.jpg)

### Iso-Seq

Iso-Seq data are in general quite accurate. Therefore, after mapping the transcripts to the reference genome by the tool of choice (e.g. minimap2, https://github.com/lh3/minimap2), the mapped transcripts can be collapsed directly, based on their structures, to generate the representative transcripts in .gtf format.

### dRNA-seq

The same protocol can be used for dRNA-seq data. However, it is not recommended because dRNA-seq data are usually quite erroneous. It is encouraged to include the more accurate short-read RNA-seq data, or accurate exon-exon junctions from other sources, to correct those potentially erroneous exon-exon junctions in dRNA-seq data. Before correction, the short-read RNA-seq data need to be mapped to the reference genome. Then the exon-exon junctions can be detected by the `NIAP_split_read2intron_v1.1.pl` script, based on the read mapping result. Error correction will be enabled by the `-junc` option. In the error correction step, each exon-exon junction of the dRNA-seq transcripts will be compared with the more accurate exon-exon junctions detected from short-read RNA-seq. If the deviation between dRNA-seq and short-read RNA-seq exon-exon junctions does not exceed a user-specified threshold, the dRNA-seq junction will be placed by the short-read RNA-seq one. The deviation is defined as below.

`deviation = |left position of dRNA-seq junction - left position of short-read RNA-seq junciton| + |right position of dRNA-seq junction - right position of short-read RNA-seq junciton|`

If users would like to examine this deviation, NIAP is able to generate a table showing this information by specifying the `-juncdist` option.

Users should be cautious that this error correction step is not perfect, because an exon-exon junction from the short-read RNA-seq data is not guaranteed to correspond to the transcript from dRNA-seq, which uses it for correction. Such correspondence relationship is computationally intractable. As a result, it may lead to over-correction and generate false transcripts.

Apart from consoidating dRNA-seq data and generating representative transcripts, NIAP can incorporate the poly-A tailing information, which is generated by `polya` function of nanopolish (https://github.com/jts/nanopolish). The corresponding NIAP option is `-polya`.

### Distributions of transcripts 5'/3' end positions

While NIAP handles Iso-Seq and dRNA-seq data differently as mentioned above, it can summarize the distributions of transcripts 5'/3' end positions for both Iso-Seq and dRNA-seq data when the `-enddist` option is specified. This information may be useful if users would like to investigate alternative transcription start sites and polyadenylation sites of transcripts with a certain structure.

### Usage

```
perl NIAP_v1.2.pl [option]
  -bam <String> Input bam file
  -junc <String> The gtf file specifying splice junctions used for correction. An optional attribute, depth, can be included.
  -injunc <Bolean> The intervals of splice junctions specified by -junc are introns. (Default: False for exons as input)
  -minjuncdist <Int> The minimal distance for splice junction correction (Default: 100 bp, empirical)
  -juncdist <Boolean> Output the distance of splice junction between the prediction and the reference (Default: false)
  -enddist <Boolean> Output the distribution of the 5' and 3' end position for each transcript (Default: false)
  -polya <String> The output from nanopolish polya
  -prefix <String> The prefix of gene and transcirpt IDs (Default: NIAP)
  -o <String> Output base
  -q <Integer> Minimal mapping quality (Inclusive; default: 0)
  -stp <String> The path to samtools
  -btp <String> The path to bedtools
  -t <Int> The number of threads to use (Default: 1)
  -h <Boolean> Help  
```

It is recommended to use minimap2 to generate the .bam file input.

A .gtf file indicating the accurate exon-exon junctions for error correction can be specified, especially for dRNA-seq data. NIAP comes with two scripts, `NIAP_split_read2intron_v1.1.pl` and `NIAP_exon2intron_v1.1.pl`, to generate a such file using short-read RNA-seq data and reference genome annotation, respectively, which indicating the borders of introns.

The shematic illustration below describe when the `-injunc` option is needed, when speciying the `-junc` option for error correction. If the exon-exon junctions .gtf file is genreated by `NIAP_split_read2intron_v1.1.pl`, the `-injunc` option is required.

![NIAP_injunc](https://user-images.githubusercontent.com/34832128/113243661-ff1a1e80-92e5-11eb-9854-3bcbcdc9ecdf.jpg)

In the error correction step, an exon-exon junction from the dRNA-seq data will be corrected if the deviation does not exceed the threadshold, as specified by the `-minjuncdist` option. Otherwise, a transcript will be filtered if it fails for any of it exon-exon junctions.

When the `-juncdist` option is enabled, it will output a talbe with the suffix `"_junc_dist.txt"` summarizing the deviations between dRNA-seq and short-read RNA-seq exon-exon junctions.

With the `-enddist` option enabled, a table with the suffix `"_end_dist.txt"` will be generated to summarize the distributions of 5' and 3' end positions, in column five and six, respectively. The distribution of end positions is recored as follows.

`<position 1>:<read count 1>,<position 2>:<read count 2>,...`

Each position refers to the position in the corresponding strand of a chromosome.

The `-polya` option takes in a table generated by the `polya` function of nanopolish.

The IDs of a gene and a transcript after merging as `"<PREFIX>_<CHR>_[01]_<LOCI_NUM>"` and `"<PREFIX>_<CHR>_[01]_<LOCUS_NUM>.<TRANSCRIPT_NUM>"`, respectively. The `"<PREFIX>"` value is `"NIAP"` by default or specified by users and can be specified by the `-prefix` option. The `"<CHR>"` value indicates the chromosome name. The `"[01]"` value indicates the strand of the locus, or transcript, with `"0"` meaning the `"+"` strand and `"1"` meaning the `"-"` one. The `"<LOCI_NUM>"` indicates the positional order of a locus at a strand. A locus is defined as a region consisted of overlapping transcripts, with each overlapping region containing at least two transcripts. The `"<TRANSCRIPT_NUM>"` is a unique identifier of a transcript in the corresponding locus.

In the output .gtf file, it will contains all representative transcripts, with attributes `"depth"` indicating the read depth of the transcrpts and `"polya_flag"` listing the poly-A annotation generated by nanopolish. The "polya_flag" attribute may be used for filtering purpose. If users wish to merge the transcripts from multiple samples and generate the read counts, the `NIAP_merge_v1.1` can be used. The file with the suffix `"_stat.txt"` summaries the read statistics.

The long-reads can be filtered based on their mapping quality score by specifying the `-q` option.

The `-stt` is not necessary if samtools is already in $PATH.

The `-btp` is not necessary if bedtools is already in $PATH.

Multithreading will be supported in latter versions.

```
perl NIAP_split_read2intron_v1.1.pl [option]
  -bam <String> Input bam file
  -o <String> Output
  -s <String> Strand of reads (Default: unstranded; RF for dUTP method; D for Nanopore direct RNA-seq or Pac-Bio Iso-Seq long-read)
  -q <Integer> Minimal mapping quality (Inclusive; default: 0)
  -btp <String> The path to bedtools
  -h <Boolean> Help
```

If users would like to generate the .gtf file of accurate exon-exon junctions from short-read RNA-seq, the short-read data should first be mapped to the reference genome using a splice-awared mapper, e.g. HISAT2 (http://daehwankimlab.github.io/hisat2/manual/). The output file in .bam format will be used as the input specified by the `-bam` option.

The output file will be in .gtf format indicating the borders of each intron and read depth by the `"depth"` attribute.

The sohrt-read data can be unstranded or stranded. It is also acceptable to use long-read data though it is generally not recommended.

The short-reads can be filtered based on their mapping quality score by specifying the `-q` option.

The `-btp` is not necessary if bedtools is already in $PATH.

```
perl NIAP_exon2intron_v1.1.pl [option]
  -i <String> Input exon gtf file
  -o <String> Basename of output files
  -btp <String> The path to bedtools
```

Alternatively, if short-read RNA-seq data is not available, the splice sites from a reference genome annotation can also be used. However, it is not recommended because splice sites in the reference genome annotation may not completely correspond to those in long-read data.

The input file is simply a .gtf file of the reference genome annotation.

Two output file will be generated, a .gtf file of the introns and a .txt plain text file summarizing the length of each intron.

The `-btp` is not necessary if bedtools is already in $PATH.

```
perl NIAP_merge_intron_v1.1.pl <input list> <output> <path to bedtools>
```

When multiple short-read datasets for the same sample are available, it is possible to merge those data using the `NIAP_merge_intron_v1.1.pl` script before performing transcript consolidation.

The script takes in a plain text file specifying the intron .gtf file(s) to be merged, one for each row. The full path(s) to the file(s) is(are) needed.

The output will be a similar .gtf file indicating the intron borders.

The `<path to bedtools>` argument is not necessary if bedtools is already in $PATH.
