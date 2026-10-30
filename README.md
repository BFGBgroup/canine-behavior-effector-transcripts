# Canine behavior-associated loci: analysis code

Analysis code accompanying the manuscript **Integrating Brain Transcriptomic and Allele-Specific Evidence to Prioritize Candidate Effector Transcripts at Canine Behavior-Associated Loci**.

Source package version: **1.0.0** ([VERSION](VERSION)).

This package contains 180 project scripts covering RNA-seq and WGS preprocessing, candidate construction, allelic analysis, fine-mapping, evidence integration, and figure generation. Original stage and script names are retained.

Before running a script, replace its `/path/to/...` placeholders, install its dependencies, and supply the required inputs. Stages require separate configuration.

## Contents

- `scripts/`: 97 Shell, 60 Python, and 23 R scripts; see the [workflow](docs/WORKFLOW.md) and [file inventory](metadata/file_inventory.tsv).
- `config/`: two sample manifests and the six-dog DoGA preprocessing task list. See [data inputs and SI Appendix references](docs/DATA_INPUTS.md).
- `patches/`: PLASMA local-changes patch and reconstructed-file checksums; see [setup instructions](docs/PLASMA.md).
- `support/figure_alignment/`: the layout-check module used by three plotting scripts.
- [Path configuration](docs/PATHS.md) and [placeholder locations](config/path_placeholders.tsv).
- [Dependencies](docs/DEPENDENCIES.md), [source records](metadata/dependency_sources.tsv), and [external PLASMA setup](docs/PLASMA.md).
- [Figure map](metadata/figure_script_map.tsv): historical directory names mapped to manuscript panels.
- [Release notes](docs/RELEASE_NOTES.md): packaging changes.

## Configuration and execution

1. Select a stage and read its launcher and worker together.
2. Set `/path/to/repository` to this repository's absolute location, then configure data and software roots using [Paths](docs/PATHS.md).
3. Supply the inputs listed in [Data inputs](docs/DATA_INPUTS.md). Update downstream paths to the intended outputs; historical timestamps do not automatically select a new run.
4. Install the listed dependencies. PLASMA and Allelome.LINK are external. PLASMA requires the implementation differences documented in [PLASMA setup](docs/PLASMA.md); an unmodified upstream checkout is not equivalent to the local implementation.
5. Check arguments, output locations, concurrency, and worker bindings. `90_historical_support` includes download and file-renaming operations and is not part of a new analysis by default.

Sequencing data, genotype panels, author-supplied GWAS statistics, analysis outputs, and the complete PLASMA source tree are not included. Packaging checks did not rerun the analyses or establish numerical reproducibility.

## Citation

Cite the accompanying manuscript named above and identify source package version 1.0.0 and the repository release or commit used. Cite the third-party methods listed in [Dependencies](docs/DEPENDENCIES.md) and [PLASMA setup](docs/PLASMA.md) separately.
