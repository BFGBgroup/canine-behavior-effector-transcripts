#!/bin/bash
set -euo pipefail


SRC_GTF="/path/to/project/legacy/reference/GCF_000002285.3_CanFam3.1_genomic_chr.gtf"
OUT_DIR="/path/to/project/legacy/five-dogs/5.cnv/0.gtf"

OUT_GTF="${OUT_DIR}/GCF_000002285.3_CanFam3.1_genomic_chr.gene_id_patched.gtf"
PATCHED_TSV="${OUT_DIR}/GCF_000002285.3_CanFam3.1_genomic_chr.gene_id_patched.features.tsv"
UNPATCHED_TSV="${OUT_DIR}/GCF_000002285.3_CanFam3.1_genomic_chr.gene_id_unpatched.features.tsv"
CHECK_LOG="${OUT_DIR}/patch_canFam_gtf_missing_geneid.check.log"

mkdir -p "${OUT_DIR}"

if [[ ! -s "${SRC_GTF}" ]]; then
  echo "[ERROR] Source GTF not found or empty: ${SRC_GTF}" >&2
  exit 1
fi

perl - "${SRC_GTF}" "${OUT_GTF}" "${PATCHED_TSV}" "${UNPATCHED_TSV}" "${CHECK_LOG}" <<'PERL'
use strict;
use warnings;

my ($src_gtf, $out_gtf, $patched_tsv, $unpatched_tsv, $check_log) = @ARGV;

sub get_attr {
    my ($attr, $key) = @_;
    return $1 if defined($attr) && $attr =~ /(?:^|;\s*)\Q$key\E\s+"([^"]+)"/;
    return "";
}

sub has_gene_id {
    my ($attr) = @_;
    return defined($attr) && $attr =~ /(?:^|;\s*)gene_id\s+"[^"]+"/;
}

sub insert_gene_id {
    my ($attr, $gene_id) = @_;
    if ($attr =~ /(transcript_id\s+"[^"]+"\s*;?)/) {
        $attr =~ s/(transcript_id\s+"[^"]+"\s*;?)/$1 gene_id "$gene_id";/;
        return $attr;
    }

    $attr =~ s/\s+$//;
    $attr .= ";" unless $attr =~ /;\s*$/;
    $attr .= qq( gene_id "$gene_id";);
    return $attr;
}

my %tid_to_gid;

open my $in1, "<", $src_gtf or die "Cannot read $src_gtf: $!\n";
while (my $line = <$in1>) {
    next if $line =~ /^#/;
    chomp $line;
    my @f = split /\t/, $line, 9;
    next unless @f == 9;

    my $tid = get_attr($f[8], "transcript_id");
    my $gid = get_attr($f[8], "gene_id");
    if ($tid ne "" && $gid ne "") {
        $tid_to_gid{$tid} = $gid;
    }
}
close $in1;

open my $in2, "<", $src_gtf or die "Cannot read $src_gtf: $!\n";
open my $out, ">", $out_gtf or die "Cannot write $out_gtf: $!\n";
open my $patched, ">", $patched_tsv or die "Cannot write $patched_tsv: $!\n";
open my $unpatched, ">", $unpatched_tsv or die "Cannot write $unpatched_tsv: $!\n";

print $patched join("\t", qw(line_number feature chr start end transcript_id patched_gene_id method original_attribute)), "\n";
print $unpatched join("\t", qw(line_number feature chr start end transcript_id reason original_attribute)), "\n";

my %feature_total;
my %feature_missing_before;
my %feature_patched;
my %feature_unpatched;
my %method_count;

my $line_no = 0;
my $total_features = 0;
my $missing_before = 0;
my $patched_count = 0;
my $unpatched_count = 0;

while (my $line = <$in2>) {
    $line_no++;

    if ($line =~ /^#/) {
        print {$out} $line;
        next;
    }

    chomp $line;
    my @f = split /\t/, $line, 9;
    if (@f != 9) {
        print {$out} $line, "\n";
        next;
    }

    $total_features++;
    my ($chr, $source, $feature, $start, $end, $score, $strand, $frame, $attr) = @f;
    $feature_total{$feature}++;

    if (!has_gene_id($attr)) {
        $missing_before++;
        $feature_missing_before{$feature}++;

        my $tid = get_attr($attr, "transcript_id");
        my $gid = "";
        my $method = "";

        if ($tid ne "" && exists $tid_to_gid{$tid}) {
            $gid = $tid_to_gid{$tid};
            $method = "transcript_id_lookup";
        } elsif ($tid =~ /^gene-/) {
            $gid = $tid;
            $method = "transcript_id_gene_prefix";
        }

        if ($gid ne "") {
            my $orig_attr = $attr;
            $attr = insert_gene_id($attr, $gid);
            $f[8] = $attr;
            $patched_count++;
            $feature_patched{$feature}++;
            $method_count{$method}++;
            print $patched join("\t", $line_no, $feature, $chr, $start, $end, $tid, $gid, $method, $orig_attr), "\n";
        } else {
            $unpatched_count++;
            $feature_unpatched{$feature}++;
            my $reason = ($tid eq "") ? "missing_transcript_id" : "no_gene_id_mapping_for_transcript_id";
            print $unpatched join("\t", $line_no, $feature, $chr, $start, $end, $tid, $reason, $attr), "\n";
        }
    }

    print {$out} join("\t", @f), "\n";
}

close $in2;
close $out;
close $patched;
close $unpatched;

open my $log, ">", $check_log or die "Cannot write $check_log: $!\n";
print $log "[patch_gtf] source_gtf=$src_gtf\n";
print $log "[patch_gtf] output_gtf=$out_gtf\n";
print $log "[patch_gtf] patched_features=$patched_tsv\n";
print $log "[patch_gtf] unpatched_features=$unpatched_tsv\n";
print $log "[patch_gtf] total_features=$total_features\n";
print $log "[patch_gtf] missing_gene_id_before=$missing_before\n";
print $log "[patch_gtf] patched_missing_gene_id=$patched_count\n";
print $log "[patch_gtf] missing_gene_id_after=$unpatched_count\n";
print $log "[patch_gtf] transcript_id_mappings_loaded=" . scalar(keys %tid_to_gid) . "\n";

print $log "[patch_gtf] patch_methods\n";
for my $method (sort keys %method_count) {
    print $log "[patch_gtf]   $method=$method_count{$method}\n";
}

print $log "[patch_gtf] feature_summary feature total missing_before patched unpatched\n";
for my $feature (sort keys %feature_total) {
    print $log join(
        " ",
        "[patch_gtf]",
        $feature,
        $feature_total{$feature} // 0,
        $feature_missing_before{$feature} // 0,
        $feature_patched{$feature} // 0,
        $feature_unpatched{$feature} // 0
    ), "\n";
}

if ($unpatched_count > 0) {
    print $log "[patch_gtf] status=WARN_UNPATCHED_FEATURES_REMAIN\n";
} else {
    print $log "[patch_gtf] status=OK\n";
}
close $log;

print "[patch_gtf] source_gtf=$src_gtf\n";
print "[patch_gtf] output_gtf=$out_gtf\n";
print "[patch_gtf] missing_gene_id_before=$missing_before\n";
print "[patch_gtf] patched_missing_gene_id=$patched_count\n";
print "[patch_gtf] missing_gene_id_after=$unpatched_count\n";
print "[patch_gtf] check_log=$check_log\n";

exit($unpatched_count > 0 ? 2 : 0);
PERL

echo
echo "[next]"
echo "1) Inspect check log:"
echo "   cat ${CHECK_LOG}"
echo "2) If missing_gene_id_after=0, rebuild STAR index with:"
echo "   STAR --runMode genomeGenerate --sjdbGTFfile ${OUT_GTF} ..."
echo "3) Then rerun STAR --quantMode GeneCounts and downstream CNV steps."
