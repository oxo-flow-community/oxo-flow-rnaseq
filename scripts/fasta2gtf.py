#!/usr/bin/env python3

# Ported from the nf-core/rnaseq 3.26.0 CUSTOM_CATADDITIONALFASTA template
# (modules/nf-core/custom/catadditionalfasta/templates/fasta2gtf.py).

# Written by Pranathi Vemuri, later modified by Jonathan Manning and released under the MIT license.

import argparse
import logging
import os
from collections.abc import Iterator
from itertools import groupby


def setup_logging() -> logging.Logger:
    """Configure logging for the script."""
    logging.basicConfig(format="%(name)s - %(asctime)s %(levelname)s: %(message)s")
    logger = logging.getLogger(__file__)
    logger.setLevel(logging.INFO)
    return logger


def parse_fasta(fasta_file: str) -> Iterator[tuple[str, str]]:
    """Parse a fasta file and yield tuples of header and sequence.

    Fasta iterator from https://www.biostars.org/p/710/#120760
    """
    with open(fasta_file) as file_handle:
        fasta_iter = (x[1] for x in groupby(file_handle, lambda line: line[0] == ">"))
        for header in fasta_iter:
            header_str = next(header)[1:].strip()
            sequence = "".join(s.strip() for s in next(fasta_iter))
            yield (header_str, sequence)


def fasta_to_gtf(fasta: str, output_file: str, biotype: str) -> None:
    """Read a fasta file and create a GTF file."""
    fasta_iter = parse_fasta(fasta)
    lines = []

    for header, sequence in fasta_iter:
        seq_name = header.split()[0].replace(" ", "_")
        line = generate_gtf_line(seq_name, len(sequence), biotype)
        lines.append(line)

    with open(output_file, "w") as file_handle:
        file_handle.writelines(lines)


def generate_gtf_line(name: str, length: int, biotype: str) -> str:
    """Generate a single GTF line given sequence name, length, and biotype."""
    biotype_attr = f' {biotype} "transgene";' if biotype else ""
    attributes = f'exon_id "{name}.1"; exon_number "1";{biotype_attr} gene_id "{name}_gene"; gene_name "{name}_gene"; gene_source "custom"; transcript_id "{name}_gene"; transcript_name "{name}_gene";\n'
    return f"{name}\ttransgene\texon\t1\t{length}\t.\t+\t.\t{attributes}"


def main() -> None:
    logger = setup_logging()
    logger.info("Starting fasta to GTF conversion.")

    parser = argparse.ArgumentParser()
    parser.add_argument("additional_fasta", help="FASTA whose sequences become transgene GTF lines")
    parser.add_argument("--out", default="transgenes.gtf", help="output GTF of transgene lines")
    parser.add_argument("--biotype", default="", help="biotype attribute value for the transgene lines")
    args = parser.parse_args()

    # Add fasta lines to GTF
    add_name = os.path.splitext(os.path.basename(args.additional_fasta))[0]
    fasta_to_gtf(args.additional_fasta, args.out, args.biotype)

    logger.info("Conversion completed successfully.")


if __name__ == "__main__":
    main()
