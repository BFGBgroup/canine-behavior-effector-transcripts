# Data inputs

This directory contains code, two sample manifests, and a six-dog DoGA preprocessing task list. It does not contain raw sequencing files, genotype panels, reference genomes, author-supplied GWAS summary statistics, intermediate results, or final figure source tables.

| Input group | Required information |
| --- | --- |
| RNA-seq | BarkBase paired FASTQ files and dog/tissue identifiers |
| WGS | Five BarkBase and six DoGA dog WGS inputs, plus reference genome |
| Reference annotation | CanFam3.1 reference FASTA/GTF and chromosome-name conventions used by the scripts |
| Transcript support | CAGE, chromatin/promoter context, poly(A) resources, and sequence inputs named in the endpoint scripts |
| Phasing | Required reference panel and the genotyping/phasing inputs named in stage 05 |
| GWAS | Complete 14-trait summary statistics supplied by the authors of the original study, trait metadata, and GRLS/PILOT LD-reference genotypes |
| Downstream analyses | The upstream tables, matrices, counts, phased genotypes, and credible-set files explicitly selected by each script |
| Figures | The specific tables and interval/reference resources referenced by each plot |

`config/dog11_wgs_manifest.tsv` retains the eleven-dog sample metadata and uses placeholder paths for files. `config/five_dogs_fastq_pair_stats.tsv` retains the nineteen dog-tissue entries and paired FASTQ naming. These are configured-example manifests based on the archived inputs, not new sample selections. Non-path fields are unchanged.

`config/STRT2_seq_fastq_task_list.tsv` is consumed by the four DoGA WGS preprocessing launchers (fastp, BWA, GATK, and normalization). Its single column, `animal_id`, contains six dogs in historical dispatch order: BC1642, JPK344, PI153, RW356, RW589, and SLK1048. It selects WGS tasks; sequencing files must be supplied separately. It does not add tissue-level replicates. Do not substitute `dog11_wgs_manifest.tsv`, which has a different schema.

## Data provenance in the accompanying SI Appendix

The SI Appendix provides the detailed sequencing accessions, reference versions, resource selection, and public access points. Use the following section and table index together with the script-specific input paths:

| Resource | SI Appendix location |
| --- | --- |
| BarkBase RNA-seq and matched five-dog WGS (PRJNA396033) | S1, Data sources and study design; Table S1, RNA-seq samples and matched WGS information |
| Six-dog DoGA WGS (PRJNA907518) and run-to-dog mapping | S1, Data sources and study design |
| CanFam3.1 (GCF_000002285.3) and RefSeq Annotation Release 105 | S2, Reference genome and candidate transcript construction; S11.1, Reference annotation and endpoint resources |
| CAGE, chromatin/promoter, and poly(A) resources | S3; S11.1; Table S2, Selection and coordinate processing of endpoint support resources; Table S4, Public software and resource access points |
| Phasing reference panel and genetic map | S4, Genotype construction and haplotype phasing; S11.3, Genotype phasing and historical RNA counts; Table S4 |
| Behavioral GWAS and GRLS/PILOT LD-reference data | S1; S7, Region definition and fine-mapping of external behavioral GWAS data; S11.5, Behavioral GWAS fine-mapping and sensitivity comparisons; Table S4; Supplementary References 5-6 |
| Software versions and implementation settings | S2-S8; S11; Table S3, Software versions and implementation details; Table S4 |
| Regional fine-mapping sensitivity comparisons | S11.5; the accompanying Dataset S1 workbook and its field descriptions |

The complete behavioral GWAS summary statistics were supplied by the authors of [Alex et al. (2025), GWAS for behavioral traits in golden retrievers identifies genes implicated in human temperament, mental health, and cognition](https://doi.org/10.1073/pnas.2421757122) and are not distributed with this source package. Contact the original-study authors for access under their sharing terms. The publicly available GRLS genotype resource is listed separately in SI Appendix Table S4 and Supplementary Reference 6.

As explained in SI Appendix S7 and S11.5, the GWAS and the separately constructed LD reference share a GRLS cohort source; the individual overlap is not established. The LD reference is not an independent validation cohort.

Preserve machine-readable field names and sample identifiers when configuring inputs.
