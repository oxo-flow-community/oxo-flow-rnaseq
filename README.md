# oxo-flow-rnaseq — RNA-seq: alignment, quantification and QC

> ★ Verified · ⇄ Official port of [`nf-core/rnaseq`](https://github.com/nf-core/rnaseq) @ `3.26.0` — same tools, same versions, same commands. Part of the [oxo-flow-community catalog](https://oxo-flow-community.github.io/).

Run a complete bulk RNA-seq analysis on paired-end reads: FastQC and fq lint
raw-read QC, TrimGalore adapter/quality trimming, STAR alignment, Picard
duplicate marking, Salmon alignment-mode quantification with tximport-merged
gene/transcript count tables and SummarizedExperiment R objects, StringTie
reference-guided transcript assembly and quantification, featureCounts gene
counts, RSeQC / dupRadar / Qualimap QC, DESeq2 sample-level QC (PCA, sample
distances, size factors), strand-specific bigWig tracks, and one final
MultiQC report with the nf-core/rnaseq custom content. The pipeline follows
the upstream default path (aligner `star_salmon`, trimmer `trimgalore`),
with every tool and command pinned to nf-core/rnaseq 3.26.0 — from raw
FASTQs to a single report with gene counts, transcript-level quantification
and per-sample QC.

## Installation

### 1. Install oxo-flow

This workflow requires **oxo-flow >= 0.12.0**. Release binary (recommended):

```bash
curl -fL -o oxo-flow.tar.gz https://github.com/Traitome/oxo-flow/releases/latest/download/oxo-flow-latest-x86_64-unknown-linux-gnu.tar.gz
tar xzf oxo-flow.tar.gz && sudo mv oxo-flow /usr/local/bin/
```

Alternatively via conda:

```bash
conda install -c bioconda oxo-flow-cli
```

(The conda package may lag behind releases; binaries for other platforms are
on the [releases page](https://github.com/Traitome/oxo-flow/releases).)

### 2. Get this workflow

```bash
git clone https://github.com/oxo-flow-community/oxo-flow-rnaseq.git
cd oxo-flow-rnaseq
```

### 3. Requirements

**Reference data (genome preparation — the upstream `PREPARE_GENOME` artifact
prep — is not ported; you must provide these):**

- `fasta` — reference genome (uncompressed FASTA)
- `gtf` — gene annotation (GTF)
- `transcript_fasta` — transcriptome FASTA (upstream `--transcript_fasta`);
  used by Salmon alignment-mode quant and StringTie
- `gene_bed` — 12-column BED of the same genes (RSeQC input)
- `chrom_sizes` — UCSC chrom.sizes file

The branch index builders ARE ported: the STAR index is built by the
`[[references]]` builder, and the HISAT2 / RSEM / Salmon indexes by
when-gated builder rules — all from the shipped fixtures when the
corresponding `config.<tool>_index` key is empty. For real data, set
`config.star_index` / `config.hisat2_index` / `config.rsem_index` /
`config.salmon_index` to your own index directories (they are symlinked in
instead of rebuilt).

**Reads:** `reads_dir/<sample>_R1.fastq.gz` and `reads_dir/<sample>_R2.fastq.gz`
(paired-end only). The sample cohort is declared in `[[sample_groups]]`.

The repository ships tiny synthetic fixtures for all of the above
(`test/fixtures/`), so the default config validates and dry-runs cleanly.

**Compute:** per-rule resources reproduce the upstream process labels
(`[rules.resources]`). STAR alignment is the most demanding rule — 12 CPUs /
72 GB; the majority of rules run on 6 CPUs / 36 GB or 1 CPU / 6 GB. Plan
node capacity accordingly.

**Tool delivery:** conda environments with pinned versions. Every rule
declares `[rules.environment] conda = "envs/<tool>.yaml"` (20 environments
in `envs/`), with each package version pinned exactly to the upstream
nf-core/rnaseq 3.26.0 module environment (e.g. `star=2.7.11b`,
`salmon=1.10.3`, `multiqc=1.33`). No containers are used — you need conda or
mamba to create and activate these environments.

## Usage

```bash
# Point OXO at your oxo-flow binary (>= 0.12.0)
export OXO=oxo-flow

# 1. Validate and lint the workflow
"$OXO" validate main.oxoflow
"$OXO" lint main.oxoflow

# 2. Preview the execution plan (fixtures are the default config)
"$OXO" dry-run main.oxoflow --samples first:1

# 3. Run with your own data: set reads_dir, fasta, gtf, transcript_fasta,
#    gene_bed, chrom_sizes, star_index and the sample list, then
"$OXO" run main.oxoflow

# Acceptance test (validate + lint + dry-run + debug):
./test/run.sh
```

The pipeline follows the upstream default path:

```
FASTQ ──> fq lint ──> FastQC (raw) ──> TrimGalore (+FastQC trim) ──> fq lint (trimmed)
      ──> STAR (star_salmon quantifier args) ──> samtools sort + index + stats
      ──> Salmon quant (alignment mode on the toTranscriptome BAM, tx2gene
           mapping from the GTF) ──> tximport merged gene/transcript tables
           ──> SummarizedExperiment RDS objects (gene + transcript)
      ──> Picard MarkDuplicates ──> samtools index + stats (markdup)
      ──> StringTie (reference-guided assembly + gene abundance)
      ──> featureCounts (gene counts) ──> biotype MultiQC table
      ──> RSeQC (bam_stat, infer_experiment, inner_distance, junction_annotation,
                 junction_saturation, read_distribution, read_duplication)
      ──> dupRadar ──> Qualimap
      ──> DESeq2 QC (PCA, sample distances, size factors)
      ──> bigWig tracks (strand-specific for stranded libraries, combined always)
      ──> MultiQC report (with the upstream custom content: fail_trimmed /
           fail_mapped tables, strandedness checks, salmon, DESeq2, sample
           merge, versions)
```

### Configuration

All upstream `params.*` of the ported path are exposed as `[config]` keys with
the upstream defaults:

| config key | upstream | default | notes |
|---|---|---|---|
| `out_dir` | `--outdir` | `results` | |
| `reads_dir` | (samplesheet) | `test/fixtures/raw` | input FASTQs |
| `strandedness` | samplesheet column | `unstranded` | `forward` / `reverse` / `unstranded`; `auto` not ported (pipeline-level setting) |
| `aligner` | `--aligner` | `star_salmon` | only `star_salmon` is ported |
| `transcript_fasta` | `--transcript_fasta` | `test/fixtures/reference/transcripts.fa` | Salmon alignment-mode quant + StringTie |
| `salmon_quant_libtype` | `--salmon_quant_libtype` | (empty) | empty = derive from `strandedness` (forward → ISF, reverse → ISR, else IU); set e.g. `A` for auto-detection |
| `min_trimmed_reads` | `--min_trimmed_reads` | `10000` | used for the MultiQC fail_trimmed table only (per-sample drop filter not ported) |
| `min_mapped_reads` | `--min_mapped_reads` | `5` | MultiQC fail_mapped table |
| `stranded_threshold` | `--stranded_threshold` | `0.8` | RSeQC strand classification |
| `unstranded_threshold` | `--unstranded_threshold` | `0.1` | RSeQC strand classification |
| `featurecounts_group_type` | `--featurecounts_group_type` | `gene_biotype` | |
| `featurecounts_feature_type` | `--featurecounts_feature_type` | `exon` | |
| `extra_fqlint_args` | `--extra_fqlint_args` | `--disable-validator P001` | |
| `skip_fastqc` | `--skip_fastqc` | `false` | |
| `skip_linting` | `--skip_linting` | `false` | |
| `skip_trimming` | `--skip_trimming` | `false` | breaks the downstream chain (see fidelity table) |
| `skip_markduplicates` | `--skip_markduplicates` | `false` | breaks the downstream chain |
| `skip_qc` | `--skip_qc` | `false` | featureCounts / RSeQC / dupRadar / Qualimap / DESeq2 QC |
| `skip_bigwig` | `--skip_bigwig` | `false` | |
| `skip_stringtie` | `--skip_stringtie` | `false` | |
| `skip_quantification_merge` | `--skip_quantification_merge` | `false` | skips the cross-sample tximport merge + SummarizedExperiment rules; leaves the MultiQC DESeq2 inputs missing (see fidelity table) |
| `skip_deseq2_qc` | `--skip_deseq2_qc` | `false` | leaves the MultiQC DESeq2 inputs missing (see fidelity table) |
| `deseq2_vst` | `--deseq2_vst` | `true` | variance stabilizing transformation for DESeq2 QC |
| `save_trimmed` / `save_align_intermeds` | `--save_trimmed` / `--save_align_intermeds` | `false` | accepted for parity; trimmed FASTQs and intermediate BAMs are kept at results/ paths regardless |

### Outputs

`results/` mirrors the upstream `outdir/` layout under `results/<aligner>/`
(aligner = `star_salmon`):

- `fastqc/raw/`, `fastqc/trim/` — FastQC HTML + zip for raw and trimmed reads
- `fq_lint/raw/`, `fq_lint/trimmed/` — fq lint reports
- `trimgalore/` — trimmed FASTQs and trimming reports
- `star_salmon/` — `<id>.Aligned.out.bam`, `<id>.Aligned.toTranscriptome.out.bam`,
  `<id>.sorted.bam`(+`.bai`), `<id>.markdup.sorted.bam`(+`.bai`),
  `<id>/quant.sf`, `<id>/quant.genes.sf`, `<id>/logs/salmon_quant.log`,
  `<id>_meta_info.json`, `<id>_lib_format_counts.json` (Salmon alignment mode),
  `salmon.merged.*` (tximport tables), `salmon.merged.gene/transcript.SummarizedExperiment.rds`,
  `salmon.merged.tx2gene.tsv`
- `star_salmon/stringtie/` — `<id>.transcripts.gtf`, `<id>.gene.abundance.txt`,
  `<id>.coverage.gtf`, `<id>.ballgown/` (moved alongside but not declared as
  a rule output)
- `star_salmon/deseq2_qc/` — `deseq2.dds.RData`, `deseq2.pca.vals.txt`,
  `deseq2.plots.pdf`, `deseq2.sample.dists.txt`, `size_factors/`,
  `star_salmon.pca.vals_mqc.tsv`, `star_salmon.sample.dists_mqc.tsv`,
  `R_sessionInfo.log`
- `star_salmon/log/` — STAR logs (`*.Log.final.out`, `*.Log.out`,
  `*.Log.progress.out`, `*.SJ.out.tab`)
- `star_salmon/samtools_stats/` — `*.sorted.bam.{stats,flagstat,idxstats}` and
  `*.markdup.sorted.bam.{stats,flagstat,idxstats}`
- `star_salmon/picard_metrics/` — `*.markdup.sorted.metrics.txt`
- `star_salmon/featurecounts/` — `*.featureCounts.tsv`(+`.summary`),
  `*.biotype_counts_mqc.tsv`, `*.biotype_counts_rrna_mqc.tsv`
- `star_salmon/rseqc/` — per-module `bam_stat/`, `infer_experiment/`,
  `inner_distance/{txt,pdf,rscript}/`, `junction_annotation/{pdf,bed,xls,log,rscript}/`,
  `junction_saturation/{pdf,rscript}/`, `read_distribution/`,
  `read_duplication/{pdf,xls,rscript}/`
- `star_salmon/dupradar/` — `scatter_plot/`, `box_plot/`, `histogram/`,
  `gene_data/`, `intercepts_slope/`, `multiqc/`, `<id>.R_sessionInfo.log`
- `star_salmon/qualimap/<id>/` — Qualimap RNA-seq QC results
- `star_salmon/bigwig/` — `<id>.forward.bigWig`, `<id>.reverse.bigWig`,
  `<id>.bigWig` (FW/REV only for stranded libraries)
- `multiqc/` — custom content files and
  `multiqc/star_salmon/multiqc_report.html`

## Source

Upstream: [`nf-core/rnaseq`](https://github.com/nf-core/rnaseq) @ `3.26.0`
(commit `e7ca46272c8f9d5ceee3f71759f4ba551d3217a4`), licensed MIT. Created
2026-08-15; this workflow may lag behind upstream releases. See `NOTICE.md`
for upstream attribution and the licensing of copied files.

## Fidelity

Commands mirror the upstream modules byte-for-byte under default parameters
(flag-for-flag, including upstream quirks such as `samtools stats` receiving
the `.bai` as a positional argument and RSeQC's stdout redirections). Upstream
process labels are reproduced as `[rules.resources]`. Every tool is pinned to
the exact upstream conda version (see `envs/`).

Known, documented deviations:

| # | upstream (3.26.0) | port | reason |
|---|---|---|---|
| 1 | Per-sample strandedness from the samplesheet (`auto` supported) | Pipeline-level `config.strandedness`, three explicit values only | oxo-flow has one config per run; `auto` needs a Salmon inference branch |
| 2 | `PREPARE_GENOME` prepares the reference artifacts (fasta, gtf, transcript_fasta, gene_bed, chrom_sizes) and builds the branch indexes (STAR / HISAT2 / RSEM / Salmon) | The artifacts are inputs; the index BUILDERS are ported: STAR via the `[[references]]` builder, HISAT2 / RSEM / Salmon via when-gated builder rules that build from the shipped fixtures when the config key is empty and symlink a user-supplied directory otherwise | The artifact prep stays excluded as infra; index building is now part of the run (see the fidelity rows below) |
| 3 | Non-default branches: `star_rsem`, `hisat2`, `--with_umi`, `--pseudo_aligner salmon` | Ported — see rows 16-23 for their deviations | The `bowtie2_salmon` aligner, the `kallisto` pseudo-aligner and the RSEM `as_quantification` mode remain excluded |
| 4 | SALMON_QUANT (alignment mode) + CUSTOM_TX2GENE + TXIMETA_TXIMPORT + SUMMARIZEDEXPERIMENT_* — the default-path quantification chain | Ported as `quantification::salmon_quant` / `tx2gene` / `tximport` / `summarizedexperiment` | The upstream 4-process chain is mirrored as 4 rules; tx2gene runs on the first sample's quant dir (upstream `.first()`); the SE process runs twice (gene + transcript) inside one rule with the upstream `--assay_names` values |
| 5 | `min_trimmed_reads` gate drops failing samples from the downstream chain | Only the MultiQC fail_trimmed table is produced | The filter is data-dependent per-sample state (n of trimmed reads), not expressible as a static DAG |
| 6 | `skip_trimming` / `skip_markduplicates` rewire the downstream inputs (QC runs on raw / sorted BAM) | `skip_trimming=true` / `skip_markduplicates=true` break the downstream chain (trimmed reads / markdup BAM are rule inputs) | oxo-flow inputs are static paths; use the defaults |
| 7 | `save_trimmed` / `save_align_intermeds` control publication; intermediates live in workdir | Trimmed FASTQs and intermediate BAMs are always kept at `results/` paths (they double as run checkpoints) | oxo-flow re-executes from declared outputs |
| 8 | RSeQC PDFs are published upstream: `*.pdf` outputs of RSEQC_JUNCTIONANNOTATION (`splicing_events_pie.pdf`, `splicing_junction_pie.pdf`), RSEQC_JUNCTIONSATURATION (`junctionSaturation_plot.pdf`), read_duplication and inner_distance — plus two zero-byte touch placeholders (`junction.pdf`, `events.pdf`) | The same PDFs are kept under `junction_annotation/pdf/`, `junction_saturation/pdf/`, `read_duplication/pdf/`, `inner_distance/pdf/` with `<id>.`-prefixed names (e.g. `<id>.junction_events.pdf`); the zero-byte `junction.pdf` / `events.pdf` touch placeholders are not produced | Layout only — the published artifact set is the same; the touch placeholders are upstream artifacts MultiQC ignores |
| 9 | `BEDTOOLS_GENOMECOV_FW/REV` swap their prefixes between forward and reverse libraries | `genomecov_fw` always emits `<id>.forward` (strand `+`), `genomecov_rev` always `<id>.reverse` (strand `-`) | With pipeline-level strandedness both rules never run together; the published artifact set is identical |
| 10 | `workflow_summary_mqc.yaml` and `methods_description_mqc.yaml` MultiQC sections (Nextflow-param rendered) | Not generated | Nextflow-specific param rendering |
| 11 | Merged-mode software versions are runtime-collated from per-process `versions.yml` | Static `nf_core_rnaseq_software_mqc_versions.yml` pinned to the env versions | Tools are pinned in `envs/*.yaml`; there are no per-process version captures in oxo-flow |
| 12 | `CUSTOM_MULTIQCCUSTOMBIOTYPE` supports `--max_biotypes` via `ext.args` | Fixed at the upstream default `100` | The upstream pipeline never sets it |
| 13 | STRINGTIE_STRINGTIE (default path, runs on the markdup BAM with `-G gtf -e`) | Ported as `quantification::stringtie` (`--fr`/`--rf` from strandedness like upstream) | The `<id>.ballgown/` directory is moved into `results/` but is not declared as a rule output (upstream emits it) |
| 14 | DESEQ2_QC (default path, runs on `salmon.merged.gene_counts_length_scaled.tsv` with `--id_col 1 --sample_suffix '' --count_col 3`, `--vst TRUE` by default) | Ported as `quantification::deseq2_qc` with the upstream header sed (label `star_salmon`) | Blind design (`design=~1`, as upstream); upstream's sample-name group decomposition (Group columns for PCA-by-group plots) is not ported — coldata is the sample IDs only; the `star_salmon.*_mqc.tsv` tables are kept in `results/` (upstream feeds them to MultiQC without publishing). Like `skip_qc` for the other QC files, `skip_deseq2_qc=true` / `skip_quantification_merge=true` leave the MultiQC rule's DESeq2 inputs missing — use the defaults |
| 15 | UMI extraction (`umitools`), BBSplit, SortMeRNA/Bowtie2 rRNA removal | Ported as when-gated rules (off by default, same gates as upstream: `with_umi` / `!skip_bbsplit` / `remove_ribo_rna` + `ribo_removal_tool`) | The four trimmed-read variants each feed the aligners, quantification and MultiQC exactly like upstream; `cat_fastq` (multi-fastq-sample branch) remains excluded |
| 16 | UMI transcriptome intermediates are unpublished Nextflow work-dir files (`{id}.bam`, `{id}.sorted.bam`, `{id}.filtered.bam` from `bam_dedup_umi`'s SAMTOOLS_SORT / UMITOOLS_PREPAREFORRSEM) | Stable canonical names: `{sample}.transcriptome.sorted.bam` → `{sample}.umi_dedup.transcriptome.sorted.bam` → `{sample}.umi_dedup.transcriptome.bam` → `{sample}.umi_dedup.transcriptome.filtered.bam` | oxo-flow has no work dirs; every intermediate is a declared output. Published names are unchanged upstream (logs, stats, prepared BAM) |
| 17 | UMI dedup outputs are tool-specific upstream (`{prefix}.dedup.bam` from UMITOOLS_DEDUP, `{prefix}.UMICollapse.bam` from UMICOLLAPSE) | All four umitools variants and the umicollapse variant write the shared path `{sample}.markdup.sorted.bam` (exclusive when-gates; downstream rules resolve one path) | Duplicate-output exclusive-gate idiom — same published artifact set per config; the `.log` / `_UMICollapse.log` logs keep their tool-specific names |
| 18 | Transcriptome-side BAM stats (`samtools_stats` for `{prefix}.umi_dedup.transcriptome.sorted.bam`) | Only the dedup-side stats are ported (`{aligner}/samtools_stats/{sample}.umi_dedup.transcriptome.sorted.bam.{stats,flagstat,idxstats}`); the coordinate-sorted index + sort-side stats are not | Sort-side stats and the index are unpublished upstream unless `--save_umi_intermeds`; dedup-side stats publish unconditionally. MultiQC excludes the transcriptome stats upstream too (`bam_dedup_umi` never mixes them into `multiqc_files`) — the port mirrors that |
| 19 | `RSEM_PREPAREREFERENCE` emits `transcripts.fa` next to the index | The `rsem_index` builder does not emit it | Nothing in the RSEM chain consumes `transcripts.fa`; the align-mode RSEM input is the toTranscriptome BAM |
| 20 | `STAR_ALIGN` passes no `--limitBAMsortRAM` | The port adds `--limitBAMsortRAM $(( effective_memory_mb * 1000000 ))` | Without it STAR's 50 GB default sort-RAM cap can fail on small hosts; the value is derived from the rule's memory like every other engine resource |
| 21 | `HISAT2_EXTRACTSPLICESITES` names the splice-site file after the GTF (`{gtf.baseName}.splice_sites.txt`) | Fixed canonical path `reference/genes.splice_sites.txt` | The port's hisat2 index builder and align rules consume it; the align command's `--rna-strandness` is rendered via a shell branch (FR forward / RF reverse / omitted unstranded — same values as the upstream `meta.strandedness` branch) |
| 22 | `SALMON_QUANT` (alignment mode) runs without `--no-version-check` | The port adds `--no-version-check` | Pre-existing port-wide deviation kept for consistency across all salmon quant rules (bam / umi / pseudo) |
| 23 | RSEM tximport reads the flat per-sample `*.isoforms.results` files; `DESEQ2_QC_RSEM` passes `--id_col 1 --sample_suffix '' --count_col 3` via the rsem deseq2 config | `tx2gene_rsem` stages the first sample's `isoforms.results` into a flat dir (same first-sample semantics as the salmon tx2gene); `deseq2_qc_rsem` passes the three args explicitly | The args equal the port script defaults but are passed explicitly for parity; the flat staging preserves the upstream first-sample `.first()` semantics |

## Not ported (metadata `excluded`)

`bowtie2_salmon` alignment (upstream `--aligner`); `kallisto` pseudo-alignment
(upstream `--pseudo_aligner`); the RSEM `as_quantification` mode; `PREPARE_GENOME`
reference-artifact prep (fasta, gtf, transcript_fasta, gene_bed, chrom_sizes are
inputs; the branch index builders ARE ported — see row 2); per-sample
`min_trimmed_reads` filtering (data-dependent per-sample state; only the
MultiQC fail_trimmed table is produced); `cat_fastq` (only active for samples
with more than one fastq pair); `auto` strandedness (per-sample inference from
the samplesheet); the DESeq2 QC sample-name group decomposition (Group columns
for PCA-by-group plots); the Nextflow-param-rendered MultiQC sections
(`workflow_summary_mqc.yaml` / `methods_description_mqc.yaml`).

## Test

```bash
bash test/run.sh
```

Runs `validate` + `lint` + `dry-run` (plus DAG-ordering and wildcard-expansion
regression checks) against the default fixture config; exits non-zero on any
failure.

## License

The port itself is Apache-2.0 (see `LICENSE`). The upstream nf-core/rnaseq
pipeline is MIT (see `LICENSE.upstream`); module scripts and MultiQC assets
copied from upstream keep their MIT headers/attribution (see `NOTICE.md`).
