#!/usr/bin/env bash
# Acceptance test for the oxo-flow-rnaseq port.
# Usage: ./test/run.sh            (uses ./main.oxoflow)
set -euo pipefail
cd "$(dirname "$0")/.."
OXO=${OXO:-oxo-flow}

echo "==> validate"
"$OXO" validate main.oxoflow

echo "==> lint (warnings are acceptable, errors are not)"
"$OXO" lint main.oxoflow

echo "==> dry-run with default config"
# oxo-flow v0.11.0 prints the plan to stderr; capture both streams
"$OXO" dry-run main.oxoflow --samples first:1 > /tmp/oxo-dryrun-$$.txt 2>&1
grep -q "would execute" /tmp/oxo-dryrun-$$.txt

echo "dry-run plan lines: $(wc -l < /tmp/oxo-dryrun-$$.txt)"

echo "==> DAG ordering: MultiQC rules come after their producers"
# The dry-run prints rules in DAG execution order ("  N. name  [run: ...]").
# multiqc has an input edge to every producer (featureCounts, RSeQC, Salmon,
# StringTie, DESeq2 QC) so it must be scheduled last; multiqc_custom_content
# must precede it. Regression guard for the {config.out_dir} expand_inputs
# wiring (literal-prefix patterns keep the DAG edges exact).
producers="bam_qc::featurecounts bam_qc::rseqc_bam_stat quantification::salmon_quant quantification::stringtie quantification::deseq2_qc"
for p in $producers; do
    pline=$(grep -nE "^  [0-9]+\. ${p}[^ ]*(  \[run.*)? *$" /tmp/oxo-dryrun-$$.txt | head -1 | cut -d: -f1 || true)
    [ -n "$pline" ] || { echo "producer rule '$p' missing from dry-run"; exit 1; }
    mline=$(grep -nE "^  [0-9]+\. multiqc(  \[run.*)? *$" /tmp/oxo-dryrun-$$.txt | head -1 | cut -d: -f1 || true)
    [ -n "$mline" ] || { echo "multiqc rule missing from dry-run"; exit 1; }
    [ "$pline" -lt "$mline" ] || { echo "multiqc scheduled before producer '$p' (line $pline vs $mline)"; exit 1; }
done
cc_line=$(grep -nE "^  [0-9]+\. multiqc_custom_content(  \[run.*)? *$" /tmp/oxo-dryrun-$$.txt | head -1 | cut -d: -f1 || true)
[ -n "$cc_line" ] || { echo "multiqc_custom_content rule missing from dry-run"; exit 1; }
[ "$cc_line" -lt "$mline" ] || { echo "multiqc_custom_content not scheduled before multiqc"; exit 1; }
echo "  multiqc_custom_content (line $cc_line) < multiqc (line $mline); all producers precede multiqc"

echo "==> debug: expanded commands contain no literal {wildcards}"
"$OXO" debug main.oxoflow 2>&1 | grep -q '{sample}' && { echo "unexpanded wildcards in debug output"; exit 1; } || true

echo "==> ribodetector branch: dry-run with remove_ribo_rna + ribo_removal_tool=ribodetector"
# Exclusive-gate sanity: flipping the ribo config to ribodetector must activate
# the seqkit_stats/ribodetector producers and the _ribodetector read-source
# twins instead of the default-path rules (no silent no-op like before the
# ribodetector wiring).
sed -e 's/^remove_ribo_rna = false$/remove_ribo_rna = true/' \
    -e 's/^ribo_removal_tool = "sortmerna"$/ribo_removal_tool = "ribodetector"/' \
    main.oxoflow > .ribodetector-test-tmp.oxoflow
grep -q '^remove_ribo_rna = true$' .ribodetector-test-tmp.oxoflow
grep -q '^ribo_removal_tool = "ribodetector"$' .ribodetector-test-tmp.oxoflow
trap 'rm -f .ribodetector-test-tmp.oxoflow' EXIT
"$OXO" dry-run .ribodetector-test-tmp.oxoflow --samples first:1 > /tmp/oxo-dryrun-ribo-$$.txt 2>&1
# Plan lines instantiate wildcards (e.g. fastq_qc::ribodetector_samples_S1),
# so match rule-name prefixes and require [run state for the positive checks.
grep -qE "^  [0-9]+\. fastq_qc::seqkit_stats[^ ]*  \[run" /tmp/oxo-dryrun-ribo-$$.txt \
    || { echo "ribodetector branch: seqkit_stats not scheduled"; cat /tmp/oxo-dryrun-ribo-$$.txt | tail -30; exit 1; }
grep -qE "^  [0-9]+\. fastq_qc::ribodetector[^ ]*  \[run" /tmp/oxo-dryrun-ribo-$$.txt \
    || { echo "ribodetector branch: ribodetector not scheduled"; exit 1; }
if grep -qE "^  [0-9]+\. fastq_qc::sortmerna[^ ]*  \[run" /tmp/oxo-dryrun-ribo-$$.txt; then
    echo "ribodetector branch: sortmerna unexpectedly scheduled"; exit 1
fi
seq_line=$(grep -nE "^  [0-9]+\. fastq_qc::seqkit_stats[^ ]*  \[run" /tmp/oxo-dryrun-ribo-$$.txt | head -1 | cut -d: -f1)
rd_line=$(grep -nE "^  [0-9]+\. fastq_qc::ribodetector[^ ]*  \[run" /tmp/oxo-dryrun-ribo-$$.txt | head -1 | cut -d: -f1)
[ "$seq_line" -lt "$rd_line" ] || { echo "ribodetector scheduled before seqkit_stats (line $rd_line vs $seq_line)"; exit 1; }
rm -f .ribodetector-test-tmp.oxoflow
trap - EXIT
echo "  seqkit_stats (line $seq_line) < ribodetector (line $rd_line); sortmerna off"

echo "PASS"
