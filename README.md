# NIAP
NIAP (Naive Isoform Analysis Pipeline) is a gernalized pipeline for analyzing short-read and long-read RNA-seq data. It is naive because it doesn't drop any transcript regarding the structure and abundance. Instead, users can perform any filtering step with the output from NIAP. There are mainly three tasks performed by NIAP, including (1) transcript merging, (2) transcript comparison and (3) transcript consolidation for long-read RNA-seq data.

## 1. Transcript merging

In the transcript merging task, the script NIAP_merge_vx.x.pl is used to merge transcripts with the same structure, either in one or multiple .gtf file(s). For 



![NIAP_merge](https://user-images.githubusercontent.com/34832128/112115696-70761500-8bf4-11eb-8191-70532280b44b.jpg)


![NIAP_annotate](https://user-images.githubusercontent.com/34832128/112115718-766bf600-8bf4-11eb-9db9-9a92619ff62b.jpg)


