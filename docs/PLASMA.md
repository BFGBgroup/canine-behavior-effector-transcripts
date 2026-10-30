# External PLASMA dependency

The complete PLASMA source tree is not bundled. The stage-10 project launcher invokes an independently installed PLASMA implementation.

Official source: [austintwang/plasma](https://github.com/austintwang/plasma).
Recorded upstream base: [c5fe317d2a59814c0a9275c2b3216e3c978f7cc9](https://github.com/austintwang/plasma/tree/c5fe317d2a59814c0a9275c2b3216e3c978f7cc9). This identifies source lineage, not a verified snapshot of every historical run.

Method reference: A. T. Wang et al., *Allele-Specific QTL Fine Mapping with PLASMA*. **American Journal of Human Genetics** 106, 170–187 (2020). [doi:10.1016/j.ajhg.2019.12.011](https://doi.org/10.1016/j.ajhg.2019.12.011).

## Installation and configuration

This package supplies [a local-changes patch](../patches/plasma-local-changes.patch) against the recorded upstream base and [SHA-256 checksums](../patches/plasma-local-changes.sha256) for the two reconstructed files. The patch was generated from the inspected local files and verified to reconstruct those files byte for byte. This is a file-level check, not an execution or numerical reproducibility test. Downloading upstream alone does not include these changes.

In a clean checkout, apply and verify the patch as follows. Replace `/path/to/repository` with this project package location. The commands below are setup instructions and do not run the analysis.

```sh
git clone https://github.com/austintwang/plasma.git plasma-local
cd plasma-local
git checkout --detach c5fe317d2a59814c0a9275c2b3216e3c978f7cc9
git apply --check /path/to/repository/patches/plasma-local-changes.patch
git apply /path/to/repository/patches/plasma-local-changes.patch
sha256sum -c /path/to/repository/patches/plasma-local-changes.sha256
```

Apply the patch only once to the specified base. Configure `PLASMA_SCRIPT` to the absolute path of `run_plasma.py` in this checkout. The patch retains the inspected local implementation, including its existing help text; it introduces no additional algorithm changes.

Configure these existing settings in `scripts/10_molecular_finemapping/8.3_no_gwas_run_expression_ase_plasma_joint/8.3.run_no_gwas_expression_ase_plasma_joint_batch.sh`:

| Setting | Placeholder | Meaning |
| --- | --- | --- |
| `PLASMA_SCRIPT` | `/path/to/software/plasma/run_plasma.py` | Absolute path to the external entry point |
| `PLASMA_PYTHON` | `/path/to/software/environments/plasma_py27/bin/python2.7` | Interpreter for that installation |

Keep the external `plasma/` package beside `run_plasma.py`. The local implementation requires `__init__.py`, `finemap.py`, `evaluator.py`, `parser.py`, `reader.py`, and `sim_locus.py`, together with NumPy, SciPy, and the `vcf` module from PyVCF.

The launcher's positional inputs are two haplotype matrices, two allele-count vectors, and an output directory. It also supplies `--total_exp_path`, `--search_mode`, `--max_causal`, `--confidence`, and `--overdispersion_global`. Its defaults are exhaustive search, a maximum of one causal variant per configuration, confidence 0.95, and global overdispersion 0. Preflight checks require the configured interpreter and external entry point to be present.

## Differences in the inspected local implementation

The local implementation differed from the recorded upstream base in two files. The differences below are captured in the supplied patch; the complete modified source tree is not bundled.

| File | Local change |
| --- | --- |
| `run_plasma.py` | Parse `overdispersion_global` and `confidence` as floating-point values and `max_causal` as an integer. |
| `run_plasma.py` | Read the overdispersion-file option from the parsed arguments and load the specified overdispersion file, rather than the haplotype file. |
| `run_plasma.py` | Represent the global overdispersion value as a floating-point vector with one entry per sample, so that the model can apply its sample mask. The retained value 0 produces a vector of zeros. |
| `run_plasma.py` | When no total-expression file is provided, use the sum of the two allele-count vectors. The retained joint-analysis launcher supplies a separate total-expression file and does not use this fallback. |
| `run_plasma.py` | Create the output directory if it is absent. The retained project launcher already creates each result directory before invoking PLASMA. |
| `plasma/finemap.py` | Correct construction of the default ASE sample mask so that it contains one Boolean `True` value per sample when no mask is supplied. |

Before running the project launcher, reconstruct and verify the implementation as above under the applicable terms, then configure its runtime environment. The unmodified upstream code is not asserted to be interchangeable with the local implementation.
