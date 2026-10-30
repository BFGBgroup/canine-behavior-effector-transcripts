#!/usr/bin/bash


set -euo pipefail
exec /usr/bin/python3.10 -B -u "/path/to/repository/scripts/12_evidence_integration/4.allelome-compare-raw-gwas/scripts/1.1_compare_dog5_allelome_raw_gwas/1.1.compare_dog5_allelome_raw_gwas.py" "$@"
