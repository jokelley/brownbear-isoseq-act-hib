# brownbear-isoseq-act-hib
Code for Long-read isoform sequencing reveals tissue-specific isoform expression between active and hibernating brown bears (Ursus arctos)

![workflow](https://github.com/Magdoll/images_public/blob/master/bear_figures/20210512_bear_figures_design.png?raw=true)


## 1. Iso-Seq Analysis

### (1a) Pooled Iso-Seq Analysis

We pooled all 18 SMRT Cells 1M into a single dataset and ran through Iso-Seq Analysis (v8.1) in [SMRTLink](https://www.pacb.com/products-and-services/analytical-software/). The input to the Iso-Seq Analysis was the pooled HiFi (CCS) reads and the output of the Iso-Seq Analysis was the high-quality, full-length transcript sequences (`hq_transcripts.fasta`). The [Iso-Seq SMRTLink Analysis Report PDF can be found here](https://github.com/jokelley/brownbear-isoseq-act-hib/blob/main/isoseq_figs/SL50279_18cell_bear_IsoSeqJob.pdf). 

### (1b) Mapping and Collapsing

We mapped the HQ transcripts to the bear genome using [minimap2](https://github.com/lh3/minimap2) and collapsed it using [Cupcake](https://github.com/Magdoll/cDNA_Cupcake), in particular the [post-Iso-Seq processing tutorial](https://github.com/Magdoll/cDNA_Cupcake/wiki/Cupcake:-supporting-scripts-for-Iso-Seq-after-clustering-step).

```
minimap2 -ax splice -t 30 -uf --secondary=no -C5 GCF_003584765.1_ASM358476v1_genomic.fna hq_transcripts.fasta > hq_transcripts.fasta.sam
sort -k 3,3 -k 4,4n hq_transcripts.fasta.sam > hq_transcripts.fasta.sorted.sam
collapse_isoforms_by_sam.py  --input hq_transcripts.fasta -s hq_transcripts.fasta.sorted.sam -c 0.99 -i 0.95 --dun-merge-5-shorter -o hq.no5merge
get_abundance_post_collapse.py hq.no5merge.collapsed cluster_report.csv
filter_away_subset.py hq.no5merge.collapsed
```

### (1c) Extracting Per-Sample Counts

We create a custom `classify_report.csv` where each full-length (FLNC) read has the proper sample labeling (ex: `CF1N`). Samples are named by the convention `[bear][tissue-F:fat,M:muscle,L:liver][1N:hibernation or 3N:active]`. 

```
python <path_to_cupcake>/post_isoseq_cluster/demux_isoseq_with_genome.py --mapped_fafq hq.no5merge.collapsed.filtered.rep.fa --read_stat hq.no5merge.collapsed.read_stat.txt --classify_csv classify_report.csv -o hq.no5merge.collapsed.filtered.mapped_fl_count.txt
```

### (1d) Classification and Filtering using SQANTI3

We used [SQANTI3](https://github.com/ConesaLab/SQANTI3/) to classify and filter the collapsed transcripts against the bear annotation.

```
python ~/GitHub/SQANTI3/sqanti3_qc.py \
                 hq.no5merge.collapsed.filtered_classification.filtered_lite.gtf \
                 GCF_003584765.1_ASM358476v1_genomic.gtf \
                 GCF_003584765.1_ASM358476v1_genomic.fna \
                 --fl_count hq.no5merge.collapsed.filtered.mapped_fl_count.txt \
                 -c splices_brownbear_shortread.tab \
                 --genename -n 20 --isoAnnotLite
     
             
python <path_to_sqanti3>/sqanti3_RulesFilter.py \
                 --faa hq.no5merge.collapsed.filtered_corrected.faa \
                 hq.no5merge.collapsed.filtered_classification.txt \
                 hq.no5merge.collapsed.filtered_corrected.fasta \
                 hq.no5merge.collapsed.filtered_corrected.gtf
```


