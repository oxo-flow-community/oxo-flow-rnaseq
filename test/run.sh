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

echo "PASS"
