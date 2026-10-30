# Packaging changes

Source package version: **1.0.0**.

- Replaced personal server paths with configurable placeholders and redirected project-script calls to bundled workers. The CNV launcher uses the no-XS worker that writes to `5.cnv/1.star`; transcript-discovery alignment remains a separate branch.
- Translated diagnostics, help, and descriptive report text into English. Removed ordinary comments and documentation strings while retaining interpreter/encoding declarations, machine-readable fields, identifiers, formatting expressions, and evidence qualifications.
- Updated twelve preprocessing batch launchers to report worker, dispatch, and incomplete-task failures. Task return values are normalized to 0 or 1; original exit codes remain in logs. Task selection, concurrency, and worker arguments are preserved.
- Restored the FEELnc launcher's complete AWK filter from the collected local script: skip GTF comments and exclude chrM records from the coding-potential training input.
- Added the six-dog DoGA WGS task list in historical dispatch order and configured the four preprocessing launchers to use it.
- Removed two internal guidance notes from the stage-1.2 required inputs and launcher arguments; GWAS thresholds, sample sizes, and candidate calculations are unchanged.
- Bundled the required figure-alignment helper. PLASMA remains external; its source lineage, local implementation differences, and configuration are described in [PLASMA setup](PLASMA.md).
- Added a PLASMA patch against the recorded upstream commit, application instructions, and SHA-256 checksums. Verified reconstruction against the two inspected local files without running PLASMA.
- Added workflow, input, dependency, and path documentation, SI Appendix cross-references, and script/figure inventories.

These changes concern the publication copy. The original analysis files and results are preserved. Packaging checks did not rerun the analyses or establish numerical reproducibility.
