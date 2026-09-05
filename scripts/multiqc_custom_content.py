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
#   workflow_summary_mqc.yaml           "Workflow Summary" params section
#   methods_description_mqc.yaml        "Methods Description" section
#
# Behaviour mirrors the upstream Groovy functions:
#   getTrimGaloreReadsAfterFiltering (R2 report: total - length-cutoff reads)
#   getInferexperimentStrandedness + calculateStrandedness
#   classifyStrand (provided = config strandedness; salmon is null in the port)
#   strandSummaryCells / strandCheckSummaryYaml / strandCheckCompositionYaml
#   multiqcNameReplacements / multiqcSampleMergeYaml / loadMultiqcAsset
#   paramsSummaryMap / paramsSummaryMultiqc (workflow_summary_mqc.yaml)
#   methodsDescriptionText (methods_description_mqc.yaml)
#
# Deviation from upstream: the merged-mode software versions file is static
# (tools are version-pinned in envs/*.yaml) instead of runtime-collated from
# per-process versions.yml files. workflow_summary_mqc.yaml ports
# paramsSummaryMap/paramsSummaryMultiqc against the upstream schema defaults
# (nextflow_schema.json, inlined below), with the 'Core Nextflow options'
# group adapted to the oxo-flow engine (engine version, invocation, launch
# dir) and a trailing group for the port-specific strandedness / chrom_sizes
# keys. methods_description_mqc.yaml renders the upstream template against
# the oxo-flow engine (the port metadata has no pipeline DOI yet).

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
# Params summary (port of upstream paramsSummaryMap / paramsSummaryMultiqc /
# methodsDescriptionText, nf-core/rnaseq nextflow_schema.json inlined)
# ---------------------------------------------------------------------------

N_A_SPAN = '<span style="color:#999999;">N/A</a>'

# Upstream methods_description_mqc.template.yaml (nf-core/rnaseq 3.26.0),
# verbatim including the ${workflow.*} / ${tool_*} placeholders.
METHODS_TEMPLATE = """<h4>Methods</h4>
<p>Data was processed using nf-core/rnaseq v${workflow.manifest.version} ${doi_text} of the nf-core collection of workflows (<a href="https://doi.org/10.1038/s41587-020-0439-x">Ewels <em>et al.</em>, 2020</a>), utilising reproducible software environments from the Bioconda (<a href="https://doi.org/10.1038/s41592-018-0046-7">Grüning <em>et al.</em>, 2018</a>) and Biocontainers (<a href="https://doi.org/10.1093/bioinformatics/btx192">da Veiga Leprevost <em>et al.</em>, 2017</a>) projects.</p>
<p>The pipeline was executed with Nextflow v${workflow.nextflow.version} (<a href="https://doi.org/10.1038/nbt.3820">Di Tommaso <em>et al.</em>, 2017</a>) with the following command:</p>
<pre><code>${workflow.commandLine}</code></pre>
<p>${tool_citations}</p>
<h4>References</h4>
<ul>
  <li>Di Tommaso, P., Chatzou, M., Floden, E. W., Barja, P. P., Palumbo, E., &amp; Notredame, C. (2017). Nextflow enables reproducible computational workflows. Nature Biotechnology, 35(4), 316-319. doi: <a href="https://doi.org/10.1038/nbt.3820">10.1038/nbt.3820</a></li>
  <li>Ewels, P. A., Peltzer, A., Fillinger, S., Patel, H., Alneberg, J., Wilm, A., Garcia, M. U., Di Tommaso, P., &amp; Nahnsen, S. (2020). The nf-core framework for community-curated bioinformatics pipelines. Nature Biotechnology, 38(3), 276-278. doi: <a href="https://doi.org/10.1038/s41587-020-0439-x">10.1038/s41587-020-0439-x</a></li>
  <li>Grüning, B., Dale, R., Sjödin, A., Chapman, B. A., Rowe, J., Tomkins-Tinch, C. H., Valieris, R., Köster, J., &amp; Bioconda Team. (2018). Bioconda: sustainable and comprehensive software distribution for the life sciences. Nature Methods, 15(7), 475–476. doi: <a href="https://doi.org/10.1038/s41592-018-0046-7">10.1038/s41592-018-0046-7</a></li>
  <li>da Veiga Leprevost, F., Grüning, B. A., Alves Aflitos, S., Röst, H. L., Uszkoreit, J., Barsnes, H., Vaudel, M., Moreno, P., Gatto, L., Weber, J., Bai, M., Jimenez, R. C., Sachsenberg, T., Pfeuffer, J., Vera Alvarez, R., Griss, J., Nesvizhskii, A. I., &amp; Perez-Riverol, Y. (2017). BioContainers: an open-source and community-driven framework for software standardization. Bioinformatics (Oxford, England), 33(16), 2580–2582. doi: <a href="https://doi.org/10.1093/bioinformatics/btx192">10.1093/bioinformatics/btx192</a></li>
  ${tool_bibliography}
</ul>
<div class="alert alert-info">
  <h5>Notes:</h5>
  <ul>
    ${nodoi_text}
    <li>The command above does not include parameters contained in any configs or profiles that may have been used. Ensure the config file is also uploaded with your publication!</li>
    <li>You should also cite all software used within this run. Check the "Software Versions" of this report to get version information.</li>
  </ul>
</div>"""


def parse_bool(text):
    """CLI config values arrive as strings; bool flags render lowercase."""
    return text == "true"


def groovy_to_string(value):
    """Render a value the way Groovy's GString interpolation would."""
    if isinstance(value, bool):
        return "true" if value else "false"
    if isinstance(value, float) and value.is_integer():
        return str(int(value))
    return str(value)


SCHEMA_GROUPS = [
    ("Input/output options", [
        ("input", "string", None),
        ("outdir", "string", None),
        ("email", "string", None),
        ("multiqc_title", "string", None),
    ]),
    ("Reference genome options", [
        ("genome", "string", None),
        ("fasta", "string", None),
        ("gtf", "string", None),
        ("gff", "string", None),
        ("gene_bed", "string", None),
        ("transcript_fasta", "string", None),
        ("additional_fasta", "string", None),
        ("splicesites", "string", None),
        ("star_index", "string", None),
        ("hisat2_index", "string", None),
        ("rsem_index", "string", None),
        ("salmon_index", "string", None),
        ("kallisto_index", "string", None),
        ("bowtie2_index", "string", None),
        ("hisat2_build_memory", "string", "200.GB"),
        ("gencode", "boolean", None),
        ("prokaryotic", "boolean", None),
        ("gffread_transcript_fasta", "boolean", None),
        ("gtf_extra_attributes", "string", "gene_name"),
        ("gtf_group_features", "string", "gene_id"),
        ("featurecounts_group_type", "string", "gene_biotype"),
        ("featurecounts_feature_type", "string", "exon"),
        ("igenomes_ignore", "boolean", None),
        ("arm", "boolean", None),
        ("igenomes_base", "string", "s3://ngi-igenomes/igenomes/"),
    ]),
    ("Read trimming options", [
        ("trimmer", "string", "trimgalore"),
        ("extra_trimgalore_args", "string", None),
        ("extra_fastp_args", "string", None),
        ("min_trimmed_reads", "integer", 10000),
    ]),
    ("Read filtering options", [
        ("bbsplit_fasta_list", "string", None),
        ("bbsplit_index", "string", None),
        ("sortmerna_index", "string", None),
        ("remove_ribo_rna", "boolean", None),
        ("ribo_removal_tool", "string", "sortmerna"),
        ("use_gpu_ribodetector", "boolean", None),
        ("ribo_database_manifest", "string", "${projectDir}/workflows/rnaseq/assets/rrna-db-defaults.txt"),
    ]),
    ("UMI options", [
        ("with_umi", "boolean", None),
        ("umi_dedup_tool", "string", "umitools"),
        ("umitools_extract_method", "string", "string"),
        ("umitools_bc_pattern", "string", None),
        ("umitools_bc_pattern2", "string", None),
        ("umi_discard_read", "integer", None),
        ("umitools_umi_separator", "string", None),
        ("umitools_grouping_method", "string", "directional"),
        ("umitools_dedup_stats", "boolean", None),
        ("umitools_dedup_primary_only", "boolean", None),
    ]),
    ("Alignment options", [
        ("aligner", "string", "star_salmon"),
        ("use_sentieon_star", "boolean", None),
        ("use_parabricks_star", "boolean", None),
        ("gpu_container_options", "string", None),
        ("pseudo_aligner", "string", None),
        ("pseudo_aligner_kmer_size", "integer", 31),
        ("bam_csi_index", "boolean", None),
        ("star_ignore_sjdbgtf", "boolean", None),
        ("salmon_quant_libtype", "string", None),
        ("min_mapped_reads", "number", 5),
        ("seq_center", "string", None),
        ("seq_platform", "string", None),
        ("stringtie_ignore_gtf", "boolean", None),
        ("extra_star_align_args", "string", None),
        ("extra_bowtie2_align_args", "string", None),
        ("extra_salmon_quant_args", "string", None),
        ("extra_kallisto_quant_args", "string", None),
        ("kallisto_quant_fraglen", "integer", 200),
        ("kallisto_quant_fraglen_sd", "integer", 200),
        ("stranded_threshold", "number", 0.8),
        ("unstranded_threshold", "number", 0.1),
    ]),
    ("Optional outputs", [
        ("save_merged_fastq", "boolean", None),
        ("save_umi_intermeds", "boolean", None),
        ("save_non_ribo_reads", "boolean", None),
        ("save_bbsplit_reads", "boolean", None),
        ("save_reference", "boolean", None),
        ("save_trimmed", "boolean", None),
        ("save_align_intermeds", "boolean", None),
        ("save_unaligned", "boolean", None),
        ("save_kraken_assignments", "boolean", None),
        ("save_kraken_unassigned", "boolean", None),
    ]),
    ("Quality Control", [
        ("extra_fqlint_args", "string", "--disable-validator P001"),
        ("deseq2_vst", "boolean", "true"),
        ("rseqc_modules", "string", "bam_stat,inner_distance,infer_experiment,junction_annotation,junction_saturation,read_distribution,read_duplication"),
        ("contaminant_screening", "string", None),
        ("contaminant_screening_input", "string", "unmapped"),
        ("kraken_db", "string", None),
        ("bracken_precision", "string", "S"),
        ("sylph_db", "string", None),
        ("sylph_taxonomy", "string", None),
    ]),
    ("Process skipping options", [
        ("skip_gtf_filter", "boolean", None),
        ("skip_gtf_transcript_filter", "boolean", None),
        ("skip_bbsplit", "boolean", True),
        ("skip_umi_extract", "boolean", None),
        ("skip_linting", "boolean", None),
        ("skip_trimming", "boolean", None),
        ("skip_alignment", "boolean", None),
        ("skip_pseudo_alignment", "boolean", None),
        ("skip_quantification_merge", "boolean", None),
        ("skip_markduplicates", "boolean", None),
        ("skip_bigwig", "boolean", None),
        ("skip_stringtie", "boolean", None),
        ("skip_fastqc", "boolean", None),
        ("use_rustqc", "boolean", False),
        ("skip_preseq", "boolean", True),
        ("skip_dupradar", "boolean", None),
        ("skip_qualimap", "boolean", None),
        ("skip_rseqc", "boolean", None),
        ("skip_biotype_qc", "boolean", None),
        ("skip_deseq2_qc", "boolean", None),
        ("skip_multiqc", "boolean", None),
        ("skip_qc", "boolean", None),
    ]),
    ("Institutional config options", [
        ("custom_config_version", "string", "master"),
        ("custom_config_base", "string", "https://raw.githubusercontent.com/nf-core/configs/master"),
        ("config_profile_name", "string", None),
        ("config_profile_description", "string", None),
        ("config_profile_contact", "string", None),
        ("config_profile_url", "string", None),
    ]),
    ("Generic options", [
        ("version", "boolean", None),
        ("publish_dir_mode", "string", "copy"),
        ("email_on_fail", "string", None),
        ("plaintext_email", "boolean", None),
        ("max_multiqc_email_size", "string", "25.MB"),
        ("monochrome_logs", "boolean", None),
        ("multiqc_config", "string", None),
        ("multiqc_logo", "string", None),
        ("multiqc_methods_description", "string", None),
        ("validate_params", "boolean", True),
        ("pipelines_testdata_base_path", "string", "https://raw.githubusercontent.com/nf-core/test-datasets/7f1614baeb0ddf66e60be78c3d9fa55440465ac8/"),
        ("trace_report_suffix", "string", None),
        ("help", "['boolean', 'string']", None),
        ("help_full", "boolean", None),
        ("show_hidden", "boolean", None),
    ]),
    ("Port options", [
        ("strandedness", "string", None),
        ("chrom_sizes", "string", None),
    ]),
]


def params_summary_map(config, engine_version, launch_dir, command_line):
    """Port of paramsSummaryMap: keep a param when the port config has the key
    AND (no schema default → value non-empty/non-false; else value != default).
    The 'Core Nextflow options' group is adapted to the oxo-flow engine:
    engine version, launch dir, invocation command line (upstream reads
    workflow.revision/runName/containerEngine/... which do not exist here).
    Config values arrive as CLI strings: "true"/"false" parse to bool so the
    false-vs-default comparison matches Groovy semantics."""
    core = {
        "version": engine_version,
        "commandLine": command_line,
        "launchDir": launch_dir,
        "projectDir": os.path.dirname(os.path.abspath(__file__)) or ".",
    }
    params_summary = {"Core Nextflow options": core}
    for group_name, group_params in SCHEMA_GROUPS:
        sub_params = {}
        for param, param_type, schema_value in group_params:
            if param not in config:
                continue
            params_value = config[param]
            if isinstance(params_value, str) and param_type in ("boolean", "['boolean', 'string']"):
                if params_value in ("true", "false"):
                    params_value = parse_bool(params_value)
            elif param_type in ("integer", "number") and isinstance(params_value, str):
                # CLI numerics arrive as strings; coerce so falsiness and
                # rendering match upstream Groovy (0 → N/A span at render).
                try:
                    params_value = int(params_value)
                except ValueError:
                    try:
                        params_value = float(params_value)
                    except ValueError:
                        pass
            params_text = groovy_to_string(params_value)
            if schema_value is not None:
                schema_text = groovy_to_string(schema_value)
                # $projectDir / ${projectDir} substitution quirk (upstream)
                if param_type == "string" and (
                    "$projectDir" in schema_text or "${projectDir}" in schema_text
                ):
                    sub_string = schema_text.replace("$projectDir", "").replace("${projectDir}", "")
                    if sub_string and params_value.startswith(sub_string):
                        schema_text = params_value
                # $params.outdir / ${params.outdir} substitution quirk (upstream)
                if param_type == "string" and (
                    "$params.outdir" in schema_text or "${params.outdir}" in schema_text
                ):
                    sub_string = schema_text.replace("$params.outdir", "").replace("${params.outdir}", "")
                    outdir = groovy_to_string(config.get("outdir", ""))
                    if outdir + sub_string == params_value:
                        schema_text = params_value
                if params_value != schema_text:
                    sub_params[param] = params_value
            else:
                # No schema default: keep unless Groovy-falsy ""/null/false.
                # Identity checks keep numeric 0 (upstream 0 != false); the
                # render-time ?: quirk then shows 0 as the N/A span.
                if params_value != "" and params_value is not None and params_value is not False:
                    sub_params[param] = params_value
        # Groovy renders booleans as lowercase true/false in GStrings.
        for param, value in list(sub_params.items()):
            if isinstance(value, bool):
                sub_params[param] = groovy_to_string(value)
        params_summary[group_name] = sub_params
    return params_summary


def workflow_summary_yaml(summary_params):
    """Port of paramsSummaryMultiqc: HTML dl rows per group, N/A span for
    Groovy-falsy values, params sorted alphabetically within each group."""
    summary_section = ""
    for group, group_params in summary_params.items():
        if not group_params:
            continue
        summary_section += f'    <p style="font-size:110%"><b>{group}</b></p>\n'
        summary_section += '    <dl class="dl-horizontal">\n'
        for param in sorted(group_params):
            rendered = group_params[param]
            if rendered is False:
                rendered = "false"
            elif not rendered:
                rendered = N_A_SPAN
            summary_section += f"        <dt>{param}</dt><dd><samp>{rendered}</samp></dd>\n"
        summary_section += "    </dl>\n"
    yaml_text = "id: 'nf-core-rnaseq-summary'\n"
    yaml_text += "description: ' - this information is collected when the pipeline is started.'\n"
    yaml_text += "section_name: 'nf-core/rnaseq Workflow Summary'\n"
    yaml_text += "section_href: 'https://github.com/nf-core/rnaseq'\n"
    yaml_text += "plot_type: 'html'\n"
    yaml_text += "data: |\n"
    yaml_text += summary_section
    return yaml_text


def methods_description_yaml(engine_version, command_line):
    """Port of methodsDescriptionText: the upstream template rendered against
    the oxo-flow engine (pipeline version 3.26.0, engine version, invocation).
    The port metadata has no pipeline DOI yet, so doi_text is empty and the
    no-DOI note is kept. tool_citations / tool_bibliography stay placeholders
    (upstream fills them from the process directives, which the port does not
    capture)."""
    doi_text = ""
    nodoi_text = '        <li>If you used nf-core/rnaseq for your analysis please cite it as above and reference the pipeline version 3.26.0.</li>\n'
    tool_citations = ""
    tool_bibliography = ""
    data = METHODS_TEMPLATE
    data = data.replace("${doi_text}", doi_text)
    data = data.replace("${workflow.manifest.version}", "3.26.0")
    data = data.replace("${workflow.nextflow.version}", engine_version)
    data = data.replace("${workflow.commandLine}", command_line)
    data = data.replace("${tool_citations}", tool_citations)
    data = data.replace("${tool_bibliography}", tool_bibliography)
    data = data.replace("${nodoi_text}", nodoi_text)
    yaml_text = "id: 'nf-core-rnaseq-methods-description'\n"
    yaml_text += "description: \"Suggested text and references to use when describing pipeline usage within the methods section of a publication.\"\n"
    yaml_text += "section_name: 'nf-core/rnaseq Methods Description'\n"
    yaml_text += "section_href: 'https://github.com/nf-core/rnaseq'\n"
    yaml_text += "plot_type: 'html'\n"
    yaml_text += "data: |\n"
    for line in data.rstrip("\n").split("\n"):
        yaml_text += "  " + line + "\n"
    return yaml_text


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
    parser.add_argument("--engine-version", required=True, help="oxo-flow engine version (upstream: workflow.nextflow.version)")
    parser.add_argument("--engine-command", default="oxo-flow run main.oxoflow", help="invocation command line (upstream: workflow.commandLine)")
    parser.add_argument("--config", action="append", default=[], metavar="KEY=VALUE", help="pipeline config key=value pairs for the Workflow Summary section (repeatable)")
    args = parser.parse_args()

    samples = sorted(s for s in args.samples.split(",") if s)
    os.makedirs(args.out_dir, exist_ok=True)

    # -- workflow_summary_mqc.yaml -----------------------------------------
    # Port config keys relevant to the upstream schema groups. Upstream
    # renders every params entry against the schema; here the CLI passes the
    # pipeline-relevant subset as key=value strings.
    config = dict(kv.split("=", 1) for kv in args.config if "=" in kv)
    summary_params = params_summary_map(
        config, args.engine_version, os.path.abspath(os.path.join(args.out_dir, os.pardir)), args.engine_command
    )
    with open(os.path.join(args.out_dir, "workflow_summary_mqc.yaml"), "w") as f:
        f.write(workflow_summary_yaml(summary_params))

    # -- methods_description_mqc.yaml --------------------------------------
    with open(os.path.join(args.out_dir, "methods_description_mqc.yaml"), "w") as f:
        f.write(methods_description_yaml(args.engine_version, args.engine_command))

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
