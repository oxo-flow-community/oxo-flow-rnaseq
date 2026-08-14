# STAR index (placeholder)

The STAR index is an upstream input artifact (PREPARE_GENOME is not
ported). Build one for your genome before running, e.g.:

    STAR --runMode genomeGenerate --genomeDir star_index \
        --genomeFastaFiles genome.fa --sjdbGTFfile genes.gtf \
        --runThreadN 8

and point `config.star_index` at the directory containing the
`SA`, `SAindex`, `genomeParameters.txt` etc. files.
