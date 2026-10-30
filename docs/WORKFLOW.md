# Workflow

The directory numbers describe analysis components; parallel branches and historical inputs mean that they are not a single serial command sequence. See `metadata/file_inventory.tsv` for the complete script list.

| Stage | Purpose | Main inputs | Main outputs |
| --- | --- | --- | --- |
| `01_rna_preprocessing` | RNA-seq preprocessing with fastp | BarkBase paired-end RNA-seq and sample manifest | Cleaned paired FASTQ files |
| `02_reference_and_discovery_alignment` | Reference attributes, STAR index, and discovery alignment | Reference FASTA/GTF and cleaned RNA-seq | Reference resources and non-WASP discovery BAMs |
| `03_candidate_catalog` | Assembly, reference comparison, candidate construction, and endpoint annotation | Discovery BAMs, StringTie models, GffCompare/FEELnc results, and external endpoint resources | Candidate loci and transcripts, candidate GTF, expanded annotation, and endpoint evidence |
| `04_wgs_preprocessing` | BarkBase and DoGA WGS preprocessing | Paired WGS FASTQ and reference genome | Aligned/deduplicated BAMs and historical per-dog VCFs |
| `05_joint_genotyping_phasing` | Joint genotyping and phasing | Eleven-dog WGS BAMs and reference phasing panel | Harmonized phased genotypes |
| `06_historical_ase_counts` | Historical RNA allele counting | RNA-seq and historical per-dog heterozygous SNPs | WASP-filtered alignments and ASEReadCounter output |
| `07_cnv` | CNV expression input, BAF input, and CaSpER | RNA expression, annotation, and allele counts | Tissue-specific CNV results |
| `08_expression_ase_integration` | Expression quantification, phase matching, and CNV-related filtering | Expanded annotation, RNA-seq, historical allele counts, phased genotypes, and CNV results | Expression matrices and GeneID-level haplotype counts |
| `09_allelome_link` | Allelic links to coding partners | Filtered haplotype counts and interval annotations | LINK pairs, evaluability information, and annotated results |
| `10_molecular_finemapping` | Exploratory PLASMA analysis | Total expression, allele counts, and phased genotypes | Molecular posterior probabilities, credible sets, and candidate catalogs |
| `11_gwas_finemapping` | Behavioral GWAS processing and SuSiE fine-mapping | Author-supplied GWAS summary statistics and GRLS/PILOT LD-reference genotypes | Regions, PIPs, credible sets, and main/strict sensitivity results |
| `12_evidence_integration` | Integrate molecular, GWAS, and LINK evidence | Molecular/GWAS credible sets, positional overlap, and LINK results | Candidate evidence tables and annotated-lncRNA context |
| `13_figures` | Plotting | Analysis outputs and source tables referenced by each script | Manuscript figures and source tables |
| `90_historical_support` | Download and naming provenance | Historical download addresses and file names | Download or renaming operations; not part of a new analysis by default |

## Relationships that must remain distinct

- The five BarkBase dogs contribute the tissue data. The six additional DoGA dogs support joint genotyping and phasing; they do not add tissue-level biological replicates.
- Candidate-discovery alignments and historical WASP/ASE-counting alignments are separate branches.
- Historical RNA allele counts are retained and subsequently matched to harmonized phasing. This workflow does not imply that RNA counting was repeated after joint genotyping.
- Candidate construction precedes downstream allelic and behavior-association integration. FEELnc and endpoint information annotate the candidate catalog according to the preserved script logic.
- In stage 10, `8.1` is the orientation-check branch and `8.3` runs the expression-plus-ASE joint analysis using an external PLASMA installation; the historical directory names are retained. See `PLASMA.md` for installation configuration and the local implementation differences.
- In stage 11, the retained LD-reference construction entry is `1.5`. Earlier `1.3`/`1.4` branches were not included as the final LD-reference workflow.
- The final integration uses both molecular/GWAS shared-variant evidence with LINK and positional GWAS credible-set overlap with LINK. Raw-GWAS comparison scripts are included as supporting/background analyses, not additional final selection paths.

## CNV alignment branch

The CNV launcher in `07_cnv/01_historical_star_no_xs` calls the bundled no-XS STAR worker in the same directory. It writes to `5.cnv/1.star`, which supplies the CNV expression-matrix step. The transcript-discovery worker in `02_reference_and_discovery_alignment/discovery_star` writes to `5.cnv/1.star-v5` and is a separate branch.

The CaSpER `rendered.R` file is an archived generated script; it is not a separate analysis entry point.

## Figures

Historical directories `figure1` through `figure5` correspond to manuscript Figures 2 through 6. The names have not been renumbered. `metadata/figure_script_map.tsv` gives the panel-level mapping. No computational script for the manuscript Figure 1 workflow diagram is asserted by this release.
