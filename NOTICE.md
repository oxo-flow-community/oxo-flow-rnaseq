oxo-flow-rnaseq
Copyright (c) 2026 oxo-flow-community

This pipeline is a port of nf-core/rnaseq
(https://github.com/nf-core/rnaseq), version 3.26.0, authored by
The nf-core/rnaseq team.

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

    http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.

---------------------------------------------------------------------
Upstream license

This port is derived from nf-core/rnaseq under the MIT license.
The upstream LICENSE is included verbatim in this repository at
LICENSE.upstream (fetched from the upstream repository at the ported
commit e7ca46272c8f9d5ceee3f71759f4ba551d3217a4, tag 3.26.0).
(Apache-2.0 §4(d): attribution notices from the Source form must be
retained.)
---------------------------------------------------------------------

Copied files that retain upstream MIT headers and attribution:

- scripts/mqc_features_stat.py — written by Senthilkumar Panneerselvam,
  adapted for nf-core/modules by Jonathan Manning (MIT)
  (modules/nf-core/custom/multiqccustombiotype/templates/mqc_features_stat.py)
- scripts/dupradar.r — from modules/nf-core/dupradar/templates/dupradar.r
  (MIT)
- scripts/gtf2bed.pl — from modules/nf-core/ea-utils/gtf2bed/templates/gtf2bed.pl
  (MIT, Copyright (c) 2011 Erik Aronesty); template variables replaced by
  command-line arguments (ARGV[0] = GTF, ARGV[1] = output prefix)
- assets/multiqc/* — from subworkflows/local/multiqc_rnaseq/assets/*
  (multiqc_config.yml, biotypes_header.txt, sample_status_header.txt,
  strand_check_summary.yaml, strand_check_composition.yaml)

Ported subworkflows/modules (commands and file naming transcribed from the
upstream sources below; no files copied verbatim):

- subworkflows/nf-core/bam_dedup_umi (UMI-tools / UMICollapse genome and
  transcriptome dedup chains, name-sort, prepare-for-rsem)
- subworkflows/nf-core/bam_sort_stats_samtools (transcriptome name/coordinate
  sorts)
- subworkflows/nf-core/bam_dedup_stats_samtools_umitools and
  bam_dedup_stats_samtools_umicollapse (dedup-side stats)
- subworkflows/nf-core/quantify_rsem (RSEM index, calculateexpression,
  custom/rsemmergecounts)
- subworkflows/nf-core/quantify_pseudo_alignment (salmon/index + salmon/quant
  pseudo mode)
- subworkflows/nf-core/quant_tximport_summarizedexperiment (tximport /
  tximeta, summarizedexperiment) — extended with the rsem and pseudo variants
- modules/nf-core/hisat2/{align,index,extractsplicesites},
  bowtie2/align, bbmap/bbsplit, sortmerna/sortmerna,
  umitools/{extract,dedup,prepareforrsem}, umicollapse/umicollapse,
  rsem/{preparereference,calculateexpression}, custom/rsemmergecounts
- workflows/rnaseq/main.nf branch wiring for the star_rsem / hisat2 /
  with_umi / pseudo_aligner configurations, and conf/modules/{
  umi_dedup,quantify_bam_salmon,quantify_rsem,quantify_pseudo_alignment,
  deseq2_qc,prepare_genome,align_hisat2,align_star}.config resource labels
  and args
