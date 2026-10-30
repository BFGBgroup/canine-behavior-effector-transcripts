# Software and dependency sources

The workflow calls FastQC, MultiQC, fastp, STAR, StringTie, GffCompare, FEELnc, BWA, samtools, bcftools, bgzip, tabix, GATK, WhatsHap, SHAPEIT5, PLINK/PLINK2, CaSpER, Allelome.LINK, PLASMA, and SuSiE. See SI Appendix S2–S8, S11, and Tables S3–S4 for reported versions and settings, and `metadata/import_inventory.tsv` for explicit Python and R imports. Environment names in the scripts are not version-lock files; configure the interpreter and executable paths for each stage.

## Command-line dependencies and configuration

| Dependency | Workflow use | Acquisition and executable configuration |
| --- | --- | --- |
| FastQC | RNA-seq quality reports in stage 01 | Obtain from the [official FastQC site](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/); make `fastqc` available on the worker's `PATH`. |
| MultiQC | RNA-seq report aggregation in stage 01 | Follow the [official MultiQC documentation](https://docs.seqera.io/multiqc); make `multiqc` available on the worker's `PATH`. |
| bgzip and tabix | Compressed VCF output and indexing in stage 04 | Both utilities are distributed with [HTSlib](https://www.htslib.org/download/); make `bgzip` and `tabix` available on the worker's `PATH`. |
| WhatsHap | Read-backed phasing in stage 05 | Follow the [official installation documentation](https://whatshap.readthedocs.io/en/stable/installation.html); configure the executable path used by the stage-05 launcher. |
| SHAPEIT5 | Common- and rare-variant phasing in stage 05 | Configure the `phase_common` and `phase_rare` executable paths used by the stage-05 launchers. |

## Allelome.LINK

Official source: [AndergassenLab/Allelome.LINK](https://github.com/AndergassenLab/Allelome.LINK).

Obtain [00_src/Allelome.LINK.R at commit 18c67e7d829b38fa777f042557c2e29b0a8e893d](https://github.com/AndergassenLab/Allelome.LINK/blob/18c67e7d829b38fa777f042557c2e29b0a8e893d/00_src/Allelome.LINK.R) and configure its location as `/path/to/software/Allelome.LINK/Allelome.LINK.R`, or update the wrapper's `LINK` setting. See the upstream [MIT license](https://github.com/AndergassenLab/Allelome.LINK/blob/18c67e7d829b38fa777f042557c2e29b0a8e893d/LICENSE).

The stage-09 wrapper directly invokes this R file with `--window-size 500`, `--total-reads 10`, `--allelic-bias 0.6`, and `--duplicates TRUE`. These explicit settings take precedence over upstream defaults. This LINK script directly loads the R package `optparse`.

## PLASMA

PLASMA is an external dependency. The local implementation differs from upstream. This package supplies a patch and file checksums for reconstruction against the specified upstream commit, rather than a complete source tree. See [External PLASMA dependency](PLASMA.md) for the official source, method citation, source lineage, local changes, and configuration requirements.

The stage-10 launcher uses `PLASMA_SCRIPT=/path/to/software/plasma/run_plasma.py` and `PLASMA_PYTHON=/path/to/software/environments/plasma_py27/bin/python2.7`. Keep the `plasma/` package beside the entry point and provide NumPy, SciPy, and PyVCF in that interpreter's environment.

## Figure-layout helper

`support/figure_alignment/audit_panel_alignment.py` is bundled with this release. Three scripts import it from `/path/to/repository/support/figure_alignment` to produce manuscript Figures 2C, 3A, and 3C (historical Figures 1C, 2A, and 2C). It checks rendered panel geometry and may stop figure export when its strict checks fail.

## Fonts and distribution scope

Scripts that request Arial require the corresponding font files and configured font environment. Font files, the complete PLASMA source tree, and complete upstream repositories are not distributed in this directory.
