#!/usr/bin/env bash
set -euo pipefail

ROOT="/path/to/project/legacy/five-dogs/data/raw-rna-data"
OUT="/path/to/project/legacy/five-dogs/data/statistics/five_dogs_fastq_pair_stats.tsv"

mkdir -p "$(dirname "$OUT")"

printf "animal_id\tterm\traw_data_name\n" > "$OUT"

find "$ROOT" -mindepth 2 -maxdepth 2 -type f \( -name "*.fastq.gz" -o -name "*.fq.gz" \) | sort | \
awk -F'/' '
BEGIN{
  OFS="\t"
}
{
  file=$NF

  if (file !~ /_[12]\.(fastq|fq)\.gz$/) next

  name=file
  sub(/\.(fastq|fq)\.gz$/, "", name)

  n=split(name, a, "_")
  if (n < 4) next

  read=a[n]
  animal_id=a[1] "_" a[2]

  term=a[3]
  for(i=4; i<=n-1; i++){
    term=term "_" a[i]
  }

  key=animal_id OFS term

  if (read=="1") r1[key]=file
  else if (read=="2") r2[key]=file
}
END{
  for (k in r1) {
    if (k in r2) print k, r1[k] "|" r2[k]
    else print k, r1[k]
  }
  for (k in r2) {
    if (!(k in r1)) print k, r2[k]
  }
}
' | sort -k1,1 -k2,2 >> "$OUT"

echo "statistics file generated: $OUT"
echo "data row count(excluding the header): $(tail -n +2 "$OUT" | wc -l)"
