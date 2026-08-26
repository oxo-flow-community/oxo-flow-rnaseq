#!/usr/bin/env python

# Workflow-level MultiQC custom content for nf-core/rnaseq 3.26.0 (port of
# subworkflows/local/multiqc_rnaseq/main.nf).
#
# Generates, into --out-dir:
#   fail_trimmed_samples_mqc.tsv        samples below --min-trimmed-reads
#   fail_mapped_samples_mqc.tsv         samples below --min-mapped-reads
#   strand_check_summary_mqc.json       strandedness inference summary table
#   strand_check_composition_mqc.json   sense/antisense/unstranded bargraph
#   name_replacement.txt                fastq simpleName -> <id>_1 / <id>_2
#   multiqc_sample_merge.yml            table_sample_merge for PE samples
#   nf_core_rnaseq_software_mqc_versions.yml  pinned tool versions
#
# Behaviour mirrors the upstream Groovy functions:
#   getTrimGaloreReadsAfterFiltering (R2 report: total - length-cutoff reads)
#   getInferexperimentStrandedness + calculateStrandedness
#   classifyStrand (provided = config strandedness; salmon is null in the port)
#   strandSummaryCells / strandCheckSummaryYaml / strandCheckCompositionYaml
#   multiqcNameReplacements / multiqcSampleMergeYaml / loadMultiqcAsset
#
# Deviation from upstream: the merged-mode software versions file is static
# (tools are version-pinned in envs/*.yaml) instead of runtime-collated from
# per-process versions.yml files. workflow_summary_mqc.yaml and
# methods_description_mqc.yaml (Nextflow-param-rendered sections) are not
# generated; see the README fidelity table.

import argparse
import json
import math
import os
import re
import sys

import yaml

# ---------------------------------------------------------------------------
# Parsing helpers (ported from the upstream Groovy)
# ---------------------------------------------------------------------------

TOTAL_READS_RE = re.compile(r"([\d.]+)\ssequences processed in total")
FILTERED_READS_RE = re.compile(r"shorter than the length cutoff[^:]+:\s*([\d.]+)")
UNIQUELY_MAPPED_RE = re.compile(r"Uniquely mapped reads %\s*\|\s*([\d.]+)%")
BOWTIE2_OVERALL_ALIGNMENT_RE = re.compile(r"(\d+\.\d+)% overall alignment rate")

INFER_RE = [
    re.compile(r"Fraction of reads failed to determine:\s([\d.]+)"),
    re.compile(r'Fraction of reads explained by "\+\+,--":\s([\d.]+)'),
    re.compile(r'Fraction of reads explained by "\+-,-\+":\s([\d.]+)'),
    re.compile(r'Fraction of reads explained by "1\++,1--,2\+-,2-\+":\s([\d.]+)'),
    re.compile(r'Fraction of reads explained by "1\+-,1-\+,2\+\+,2--":\s([\d.]+)'),
]

REGEX_ESCAPE_RE = re.compile(r'[\\^$.|?*+()\[\]{}/]')


def get_trimgalore_reads_after_filtering(log_file):
    """Upstream getTrimGaloreReadsAfterFiltering: total - length-cutoff reads."""
    total_reads = 0.0
    filtered_reads = 0.0
    for line in open(log_file):
        m = TOTAL_READS_RE.search(line)
        if m:
            total_reads = float(m.group(1))
        m = FILTERED_READS_RE.search(line)
        if m:
            filtered_reads = float(m.group(1))
    return total_reads - filtered_reads


def get_inferexperiment_strandedness(infer_file, stranded_threshold, unstranded_threshold):
    """Upstream getInferexperimentStrandedness + calculateStrandedness."""
    forward = reverse = unstranded = 0.0
    for line in open(infer_file):
        m = INFER_RE[0].search(line)
        if m:
            unstranded = float(m.group(1)) * 100
        m = INFER_RE[1].search(line)
        if m:
            forward = float(m.group(1)) * 100
        m = INFER_RE[2].search(line)
        if m:
            reverse = float(m.group(1)) * 100
        m = INFER_RE[3].search(line)
        if m:
            forward = float(m.group(1)) * 100
        m = INFER_RE[4].search(line)
        if m:
            reverse = float(m.group(1)) * 100

    total = forward + reverse + unstranded
    stranded_total = forward + reverse

    inferred = "undetermined"
    if stranded_total > 0:
        forward_prop = forward / stranded_total
        reverse_prop = reverse / stranded_total
        diff = abs(forward_prop - reverse_prop)
        if forward_prop >= stranded_threshold:
            inferred = "forward"
        elif reverse_prop >= stranded_threshold:
            inferred = "reverse"
        elif diff <= unstranded_threshold:
            inferred = "unstranded"

    return {
        "inferred_strandedness": inferred,
        "forwardFragments": (forward / total) * 100 if total else 0.0,
        "reverseFragments": (reverse / total) * 100 if total else 0.0,
        "unstrandedFragments": (unstranded / total) * 100 if total else 0.0,
    }


def inference_certainty(analysis):
    """Upstream inferenceCertainty: inferred direction's share of the stranded pool."""
    if not analysis:
        return None
    fwd = analysis["forwardFragments"]
    rev = analysis["reverseFragments"]
    stranded = fwd + rev
    if stranded == 0:
        return None
    s = analysis["inferred_strandedness"]
    if s == "forward":
        return round_one_decimal((fwd / stranded) * 100)
    if s == "reverse":
        return round_one_decimal((rev / stranded) * 100)
    return None


def round_one_decimal(v):
    """Upstream roundOneDecimal: Math.round(v * 10) / 10.0."""
    if v is None:
        return None
    return math.floor(v * 10 + 0.5) / 10.0


def load_multiqc_asset(asset_path):
    """Upstream loadMultiqcAsset: YAML parse, drop top-level '_' keys."""
    with open(asset_path) as f:
        parsed = yaml.safe_load(f)
    return {k: v for k, v in parsed.items() if not k.startswith("_")}


def multiqc_sample_merge_yaml_pattern(sample_id, read):
    """Upstream multiqcSampleMergeYamlPattern: escape regex metachars, single-quote."""
    esc = REGEX_ESCAPE_RE.sub(lambda m: "\\" + m.group(0), sample_id)
    esc = esc.replace("'", "''")
    return f"    - type: regex\n      pattern: '(?<=^{esc})_{read}$'"


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out-dir", required=True)
    parser.add_argument("--aligner", default="star_salmon", help="aligner subdir under out_dir (upstream: params.aligner)")
    parser.add_argument("--samples", required=True, help="comma-joined sample ids")
    parser.add_argument("--reads-dir", required=True)
    parser.add_argument("--strandedness", required=True, choices=["forward", "reverse", "unstranded"])
    parser.add_argument("--min-trimmed-reads", type=float, required=True)
    parser.add_argument("--min-mapped-reads", type=float, required=True)
    parser.add_argument("--stranded-threshold", type=float, required=True)
    parser.add_argument("--unstranded-threshold", type=float, required=True)
    parser.add_argument("--status-header", required=True)
    parser.add_argument("--summary-asset", required=True)
    parser.add_argument("--composition-asset", required=True)
    args = parser.parse_args()

    samples = sorted(s for s in args.samples.split(",") if s)
    os.makedirs(args.out_dir, exist_ok=True)

    # -- fail_trimmed_samples_mqc.tsv --------------------------------------
    # The R2 trimming report is parsed (upstream: trim_log[-1] for PE).
    trimmed_rows = []
    for sid in samples:
        report = os.path.join(
            os.path.dirname(args.out_dir), "trimgalore", f"{sid}_trimmed_2.fastq.gz_trimming_report.txt"
        )
        if not os.path.exists(report):
            report = os.path.join(
                os.path.dirname(args.out_dir), "trimgalore", f"{sid}_trimmed_1.fastq.gz_trimming_report.txt"
            )
        if not os.path.exists(report):
            continue
        n = get_trimgalore_reads_after_filtering(report)
        if n <= args.min_trimmed_reads:
            trimmed_rows.append(f"{sid}\t{n:g}\n")
    with open(os.path.join(args.out_dir, "fail_trimmed_samples_mqc.tsv"), "w") as f:
        f.write("Sample\tReads after trimming\n")
        for row in trimmed_rows:
            f.write(row)

    # -- fail_mapped_samples_mqc.tsv ---------------------------------------
    # Parent header lines (sample_status_header.txt) + column row, then one
    # row per failing sample (upstream: skip=status_header_lines on merge).
    mapped_rows = []
    for sid in samples:
        log = os.path.join(os.path.dirname(args.out_dir), args.aligner, "log", f"{sid}.Log.final.out")
        star_log = os.path.exists(log)
        if not star_log:
            # bowtie2_salmon: upstream ALIGN_BOWTIE2 emits {id}.bowtie2.log
            # ("N% overall alignment rate") and the multiqc_rnaseq subworkflow
            # parses it with getBowtie2PercentMapped into the same
            # fail_mapped table. The header line below stays the hardcoded
            # STAR text — that is an upstream quirk we preserve verbatim.
            log = os.path.join(os.path.dirname(args.out_dir), args.aligner, "log", f"{sid}.bowtie2.log")
        if not os.path.exists(log):
            continue
        percent = None
        for line in open(log):
            m = UNIQUELY_MAPPED_RE.search(line)
            if not m and not star_log:
                m = BOWTIE2_OVERALL_ALIGNMENT_RE.search(line)
            if m:
                percent = float(m.group(1))
        if percent is not None and percent < args.min_mapped_reads:
            mapped_rows.append(f"{sid}\t{percent}\n")
    with open(args.status_header) as f:
        header_text = f.read()
    with open(os.path.join(args.out_dir, "fail_mapped_samples_mqc.tsv"), "w") as f:
        f.write(header_text)
        f.write("Sample\tSTAR uniquely mapped reads (%)\n")
        for row in mapped_rows:
            f.write(row)

    # -- strand checks -----------------------------------------------------
    summary_static = load_multiqc_asset(args.summary_asset)
    composition_static = load_multiqc_asset(args.composition_asset)
    composition_static.update(
        {k: summary_static[k] for k in ("parent_id", "parent_name", "parent_description")}
    )
    header_keys = list(summary_static["headers"].keys())

    rows = []
    for sid in samples:
        infer_file = os.path.join(
            os.path.dirname(args.out_dir), args.aligner, "rseqc", "infer_experiment", f"{sid}.infer_experiment.txt"
        )
        if not os.path.exists(infer_file):
            continue
        rseqc = get_inferexperiment_strandedness(
            infer_file, args.stranded_threshold, args.unstranded_threshold
        )
        provided = args.strandedness
        status = "pass" if provided == rseqc["inferred_strandedness"] else "fail"
        rows.append((sid, provided, status, None, rseqc))

    # Summary table (upstream strandSummaryCells + strandCheckSummaryYaml)
    data = {}
    for sid, provided, status, salmon, rseqc in rows:
        raw = {
            "provided": provided,
            "salmon_inferred": (salmon or {}).get("inferred_strandedness", "-"),
            "salmon_pct": inference_certainty(salmon),
            "salmon_s": round_one_decimal((salmon or {}).get("forwardFragments")),
            "salmon_a": round_one_decimal((salmon or {}).get("reverseFragments")),
            "salmon_u": round_one_decimal((salmon or {}).get("unstrandedFragments")),
            "rseqc_inferred": rseqc["inferred_strandedness"],
            "rseqc_pct": inference_certainty(rseqc),
            "rseqc_s": round_one_decimal(rseqc["forwardFragments"]),
            "rseqc_a": round_one_decimal(rseqc["reverseFragments"]),
            "rseqc_u": round_one_decimal(rseqc["unstrandedFragments"]),
            "status": status,
        }
        unknown = set(raw) - set(header_keys)
        if unknown:
            raise SystemExit(f"strand_check_summary.yaml headers do not declare columns: {sorted(unknown)}")
        cells = {k: raw[k] for k in header_keys if raw[k] is not None}
        data[sid] = cells
    summary_json = dict(summary_static)
    summary_json["data"] = data
    with open(os.path.join(args.out_dir, "strand_check_summary_mqc.json"), "w") as f:
        json.dump(summary_json, f, indent=4)

    # Composition bargraph (upstream strandCheckCompositionYaml)
    rseqc_data = {}
    for sid, _provided, _status, _salmon, rseqc in rows:
        rseqc_data[sid] = {
            "Sense": round_one_decimal(rseqc["forwardFragments"]),
            "Antisense": round_one_decimal(rseqc["reverseFragments"]),
            "Unstranded": round_one_decimal(rseqc["unstrandedFragments"]),
        }
    composition_json = dict(composition_static)
    composition_json["data"] = rseqc_data
    with open(os.path.join(args.out_dir, "strand_check_composition_mqc.json"), "w") as f:
        json.dump(composition_json, f, indent=4)

    # -- name_replacement.txt ----------------------------------------------
    # Upstream multiqcNameReplacements: fastq simpleName -> <id>_1 / <id>_2,
    # skipped when the simpleName already equals the sample id.
    mappings = []
    for sid in samples:
        r1 = os.path.join(args.reads_dir, f"{sid}_R1.fastq.gz")
        r2 = os.path.join(args.reads_dir, f"{sid}_R2.fastq.gz")
        simple1 = os.path.basename(r1)[: -len(".fastq.gz")]
        simple2 = os.path.basename(r2)[: -len(".fastq.gz")]
        if simple1 != sid:
            mappings.append(f"{simple1}\t{sid}_1")
            mappings.append(f"{simple2}\t{sid}_2")
    with open(os.path.join(args.out_dir, "name_replacement.txt"), "w") as f:
        for m in mappings:
            f.write(m + "\n")

    # -- multiqc_sample_merge.yml ------------------------------------------
    # Upstream multiqcSampleMergeYaml: PE sample ids, fixed-length lookbehind
    # so ids ending in _1 / _2 are not wrongly collapsed.
    r1_lines = "\n".join(multiqc_sample_merge_yaml_pattern(s, 1) for s in samples)
    r2_lines = "\n".join(multiqc_sample_merge_yaml_pattern(s, 2) for s in samples)
    with open(os.path.join(args.out_dir, "multiqc_sample_merge.yml"), "w") as f:
        if not samples:
            f.write("table_sample_merge: {}\n")
        else:
            f.write("table_sample_merge:\n")
            f.write('  "Read 1":\n')
            f.write(r1_lines + "\n")
            f.write('  "Read 2":\n')
            f.write(r2_lines + "\n")

    # -- nf_core_rnaseq_software_mqc_versions.yml --------------------------
    # Static version manifest matching the envs/*.yaml pins (upstream merged
    # mode collates per-process versions.yml at runtime instead).
    versions = {
        "FASTQC": {"fastqc": "0.12.1"},
        "TRIMGALORE": {"trim-galore": "2.1.0"},
        "FQ_LINT": {"fq": "0.12.0"},
        "FQ_LINT_AFTER_TRIMMING": {"fq": "0.12.0"},
        "STAR_ALIGN": {"star": "2.7.11b", "samtools": "1.21", "htslib": "1.21", "gawk": "5.1.0"},
        "SAMTOOLS_SORT": {"samtools": "1.23.1", "htslib": "1.23.1"},
        "SAMTOOLS_INDEX": {"samtools": "1.23.1", "htslib": "1.23.1"},
        "SAMTOOLS_STATS": {"samtools": "1.23.1", "htslib": "1.23.1"},
        "SAMTOOLS_FLAGSTAT": {"samtools": "1.23.1", "htslib": "1.23.1"},
        "SAMTOOLS_IDXSTATS": {"samtools": "1.23.1", "htslib": "1.23.1"},
        "PICARD_MARKDUPLICATES": {"picard": "3.4.0"},
        "HISAT2_ALIGN": {"hisat2": "2.2.1", "samtools": "1.20"},
        "BBMAP_BBSPLIT": {"bbmap": "39.18"},
        "SORTMERNA": {"sortmerna": "4.3.7"},
        "BOWTIE2_ALIGN": {"bowtie2": "2.5.4", "htslib": "1.21", "samtools": "1.21"},
        "RSEM_CALCULATEEXPRESSION": {"rsem": "1.3.3"},
        "UMITOOLS_EXTRACT": {"umi_tools": "1.1.6", "pysam": "0.22.0"},
        "UMITOOLS_DEDUP": {"umi_tools": "1.1.6", "pysam": "0.22.0"},
        "UMITOOLS_PREPAREFORRSEM": {"umi_tools": "1.1.6", "pysam": "0.22.0"},
        "UMICOLLAPSE": {"umicollapse": "1.1.0"},
        "SUBREAD_FEATURECOUNTS": {"subread": "2.0.6"},
        "CUSTOM_MULTIQCCUSTOMBIOTYPE": {"python": "3.12.12"},
        "RSEQC_BAMSTAT": {"rseqc": "5.0.4", "r-base": "4.3"},
        "RSEQC_INFEREXPERIMENT": {"rseqc": "5.0.4", "r-base": "4.3"},
        "RSEQC_INNERDISTANCE": {"rseqc": "5.0.4", "r-base": "4.3"},
        "RSEQC_JUNCTIONANNOTATION": {"rseqc": "5.0.4", "r-base": "4.3"},
        "RSEQC_JUNCTIONSATURATION": {"rseqc": "5.0.4", "r-base": "4.3"},
        "RSEQC_READDISTRIBUTION": {"rseqc": "5.0.4", "r-base": "4.3"},
        "RSEQC_READDUPLICATION": {"rseqc": "5.0.4", "r-base": "4.3"},
        "DUPRADAR": {"bioconductor-dupradar": "1.38.0"},
        "SAMTOOLS_SORT_QUALIMAP": {"samtools": "1.23.1", "htslib": "1.23.1"},
        "QUALIMAP_RNASEQ": {"qualimap": "2.3"},
        "BEDTOOLS_GENOMECOV_FW": {"bedtools": "2.31.1"},
        "BEDTOOLS_GENOMECOV_REV": {"bedtools": "2.31.1"},
        "BEDTOOLS_GENOMECOV_COMBINED": {"bedtools": "2.31.1"},
        "UCSC_BEDCLIP": {"ucsc-bedclip": "377"},
        "UCSC_BEDGRAPHTOBIGWIG": {"ucsc-bedgraphtobigwig": "469"},
        "MULTIQC": {"multiqc": "1.33"},
    }
    with open(os.path.join(args.out_dir, "nf_core_rnaseq_software_mqc_versions.yml"), "w") as f:
        f.write("# Software versions for this run (pinned in envs/*.yaml)\n")
        yaml.safe_dump(versions, f, sort_keys=False)


if __name__ == "__main__":
    sys.exit(main())
