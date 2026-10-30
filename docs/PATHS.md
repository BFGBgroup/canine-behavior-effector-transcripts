# Placeholder paths

The root prefixes below are examples, not existing server locations or download addresses. Paths under different roots refer to distinct historical data/resource layouts. Configure them separately; do not replace all roots with a single data directory without checking the remaining suffixes.

| Placeholder root | Meaning |
| --- | --- |
| `/path/to/repository` | Absolute location of this source release; known project-script calls point to the bundled files |
| `/path/to/project/legacy` | Shared legacy root containing `five-dogs` and `reference` |
| `/path/to/project/analysis_20260609` | Historical analysis layout used by the CaSpER stage |
| `/path/to/project/analysis_20260622` | Earlier analysis layout referenced by a preparation step |
| `/path/to/project/analysis_20260709` | Historical candidate, genotype/phasing, and GWAS resources |
| `/path/to/project/analysis_20260817` | Later expression, molecular, integration, and figure resources |
| `/path/to/project/doga_six_dogs` | DoGA six-dog WGS resources |
| `/path/to/software/miniconda3` | Main Conda installation |
| `/path/to/software/legacy_miniconda3` | Separate legacy Conda installation |
| `/path/to/software/environments` | Named standalone software environments |
| `/path/to/software/legacy_envs` | Additional legacy environment layout |
| `/path/to/software/venvs` | Python virtual environments used by plotting |
| `/path/to/software/plasma` | External PLASMA installation containing `run_plasma.py` and its adjacent package; see `PLASMA.md` |
| `/path/to/software/Allelome.LINK` | External Allelome.LINK R modules |
| `/path/to/repository/support/figure_alignment` | Bundled figure-alignment validation module |
| `/path/to/software/R/figure3-library` | Additional R library location |
| `/path/to/fonts/Arial` | Required Arial font files for the plotting scripts that explicitly check them |

`config/path_placeholders.tsv` lists literal placeholder locations, retaining complete shell variable expressions such as `${ANIMAL}/${TERM}`. Additional paths are assembled from these roots and relative suffixes at runtime; inspect the relevant launcher and worker together.

Standard operating-system paths, such as `/usr/bin/env`, `/bin/bash`, `/dev/null`, and `/tmp`, remain where used.

The BarkBase WGS batch launchers use the bundled workers under `scripts/04_wgs_preprocessing/barkbase5`.

BarkBase WGS input and output paths use `/path/to/project/legacy/five-dogs/data/WGS`, including `WGS_fastp`, `WGS_bwa`, `WGS_gatk`, and other retained stage suffixes. These are data paths, not locations under `scripts/90_historical_support`. The WGS workers, manifest builder, historical SNP-count input, FASTQ-renaming helper, and sample manifest use this same data root.

The two sample manifests in `config/` show the expected input layout. Supply a configured copy at the path expected by each entry point, or explicitly change that entry point's manifest setting. Manifest paths still identify the corresponding data roots; they are not automatically redirected to the example files in this repository.

The four DoGA WGS preprocessing launchers directly use `/path/to/repository/config/STRT2_seq_fastq_task_list.tsv`. Replace the repository root as usual. The bundled file contains the `animal_id` column and preserves the six historical dispatch identifiers and their order; the worker data paths remain separately configured. See `DATA_INPUTS.md` for the task-list schema and its distinction from the eleven-dog manifest.

Old timestamps and output basenames are preserved because they distinguish actual inputs. A new run may produce different directory names. Configure downstream inputs explicitly instead of assuming that a historical timestamp selects your new output.
