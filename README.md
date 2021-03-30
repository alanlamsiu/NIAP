# NIAP
NIAP (Naive Isoform Analysis Pipeline) is a gernalized pipeline for analyzing short-read and long-read RNA-seq data. It is naive because it doesn't drop any transcript regarding the structure and abundance. Instead, users can perform any filtering step with the output from NIAP. There are mainly three tasks performed by NIAP, including (1) transcript merging, (2) transcript annotation, (3) transcript comparison and (4) transcript consolidation for long-read RNA-seq data.

## Installation

`git clone https://github.com/alanlamsiu/NIAP.git`

Scripts performing the four tasks mentioned above can be found in ./NIAP/script.

## Computational environment

All scripts have been tested in Perl (v5.22.1) and Ubuntu.

bedtools (https://bedtools.readthedocs.io/en/latest/) is requierd. The path to bedtools can either be put in the $PATH or speecified when running idividual scripts.

## 1. Transcript merging

In the transcript merging task, the script NIAP_merge_v1.1.pl is used to merge transcripts with the same structure, either in one or multiple .gtf file(s). Multiexonic transcirpts with the exact set of exon-exon junction(s) and being transcribed from the same direction will be merged into one representative transcript. On the other hand, overlapping monoexonic transcripts will be merged to generate a representative transcript. In both cases, the longest 5' and 3' ends, among the transcripts to be merged, will be taken for the representative transcript. A schematic illustration is shown below.

![NIAP_merge](https://user-images.githubusercontent.com/34832128/112115696-70761500-8bf4-11eb-8191-70532280b44b.jpg)

### Usage

```
perl NIAP_merge_v1.1.pl [option]
  -in <String> Input list of .gtf file(s)
  -out <String> Output base
  -exp <String> The atribute for the getting expression value (cov, FPKM, TPM or etc.)
  -prefix <String> The prefix of gene and transcript IDs (Default: NIAP)
  -btp <String> The path to bedtools
```

The input is a plain text file specifying the .gtf file(s) to be merged, one for each row. The full path(s) to the file(s) is(are) needed. There can be a single file or multiple files. A .gtf file of the representative transcirpts after merging and a plain text file, with a sufix of `"_id.txt"`, indicating the source of each representative transcripts will be generated. The `"_id.txt"` file is consisted of three columns, where the first two columns indicating the gene ID and transcript ID of the representative transcripts, the third one indicating `","` delineated original IDs before merging. If the `-exp` option is specified for extracting the expression levele, or any other quantitative properties, of the transcripts, a summary table with the suffix `_exp.txt` will be generated. The IDs of a gene and a transcript after merging as `"<PREFIX>_<CHR>_[01]_<LOCI_NUM>"` and `"<PREFIX>_<CHR>_[01]_<LOCUS_NUM>.<TRANSCRIPT_NUM>"`, respectively. The `"<PREFIX>"` value is `"NIAP"` by default or specified by users and can be specified by the `-prefix` option. The `"<CHR>"` value indicates the chromosome name. The `"[01]"` value indicates the strand of the locus, or transcript, with `"0"` meaning the `"+"` strand and `"1"` meaning the `"-"` one. The `"<LOCI_NUM>"` indicates the positional order of a locus at a strand. A locus is defined as a region consisted of overlapping transcripts, with each overlapping region containing at least two transcripts. The `"<TRANSCRIPT_NUM>"` is a unique identifier of a transcript in the corresponding locus. The `-btp` is not necessary if bedtools is already in $PATH.
  
## 2. Transcript comparison




![NIAP_annotate](https://user-images.githubusercontent.com/34832128/112115718-766bf600-8bf4-11eb-9db9-9a92619ff62b.jpg)


