# NIAP
NIAP (Naive Isoform Analysis Pipeline) is a gernalized pipeline for analyzing short-read and long-read RNA-seq data. It is naive because it doesn't drop any transcript regarding the structure and abundance. Instead, users can perform any filtering step with the output from NIAP. There are mainly three tasks performed by NIAP, including (1) transcript merging, (2) transcript annotation, (3) transcript comparison and (4) transcript consolidation for long-read RNA-seq data.

## Installation

`git clone https://github.com/alanlamsiu/NIAP.git`

Scripts performing the four tasks mentioned above can be found in ./NIAP/script.

## Computational environment

All scripts have been tested in Perl (v5.22.1) and Ubuntu.

bedtools (https://bedtools.readthedocs.io/en/latest/) is requierd. The path to bedtools can either be put in the $PATH or speecified when running idividual scripts.

## 1. Transcript merging

In the transcript merging task, the script NIAP_merge_vx.x.pl is used to merge transcripts with the same structure, either in one or multiple .gtf file(s). Multiexonic transcirpts with the exact set of exon-exon junction(s) and being transcribed from the same direction will be merged into one representative transcript. On the other hand, overlapping monoexonic transcripts will be merged to generate a representative transcript. In both cases, the longest 5' and 3' ends, among the transcripts to be merged, will be taken for the representative transcript. A schematic illustration is shown below.

![NIAP_merge](https://user-images.githubusercontent.com/34832128/112115696-70761500-8bf4-11eb-8191-70532280b44b.jpg)

## 2. Transcript comparison




![NIAP_annotate](https://user-images.githubusercontent.com/34832128/112115718-766bf600-8bf4-11eb-9db9-9a92619ff62b.jpg)


