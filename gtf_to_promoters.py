#!/usr/bin/env python3
"""
Extract promoter regions from a GTF or BED12 file into a promoters BED file.

For each gene:
  + strand: promoter = [start - upstream, start + downstream]
  - strand: promoter = [end - downstream, end + upstream]

Coordinates are clamped to a minimum of 1.

Gene boundaries are determined automatically depending on input format:

  GTF input:
    - If the file contains "gene" feature-type rows, those are used directly
      (matches Ensembl GTFs from recent releases).
    - Otherwise (e.g. older Ensembl/iGenomes GTFs that only go down to
      "transcript" level), gene extent is inferred by taking the min(start)
      and max(end) across all "transcript" rows sharing the same gene_id.

  BED12 input (e.g. iGenomes genes.bed, one row per transcript):
    - gene_id is derived from the transcript name (column 4) by stripping a
      trailing ".<number>" transcript suffix (e.g. AT1G01010.1 -> AT1G01010).
    - Gene extent is the min(chromStart) / max(chromEnd) across all
      transcripts sharing that gene_id.
    - BED chromStart (0-based) is converted to 1-based internally so the
      output follows the same convention as the GTF path.

Format is auto-detected from the file extension (.bed vs .gtf/.gtf.gz),
or can be forced with --format.

Usage:
    ./gtf_to_promoters.py -i genes.gtf.gz -o promoters.bed
    ./gtf_to_promoters.py -i genes.gtf -u 1000 -d 200 -o promoters.bed
    ./gtf_to_promoters.py -i genes.bed -o promoters.bed
    ./gtf_to_promoters.py -i genes.txt --format bed -o promoters.bed
"""

import argparse
import gzip
import re
import sys
from collections import OrderedDict

GENE_ID_RE = re.compile(r'gene_id "([^"]+)"')
BED_TRANSCRIPT_SUFFIX_RE = re.compile(r'^(.*)\.\d+$')

CHR_RENAME = {
    "Pt": "Cp",  # plastid -> Cp, matches original script's chrPt -> chrCp
}


def open_maybe_gz(path):
    if path.endswith(".gz"):
        return gzip.open(path, "rt")
    return open(path, "r")


def detect_format(path, forced_format):
    if forced_format:
        return forced_format
    lower = path.lower()
    if lower.endswith(".bed") or lower.endswith(".bed.gz"):
        return "bed"
    return "gtf"


def format_chrom(raw_chr, add_chr_prefix):
    chrom = raw_chr
    if add_chr_prefix:
        chrom = f"chr{chrom}"
        for old, new in CHR_RENAME.items():
            if chrom == f"chr{old}":
                chrom = f"chr{new}"
    return chrom


def get_gene_id_gtf(attr_field, line):
    m = GENE_ID_RE.search(attr_field)
    if not m:
        sys.exit(f"Could not obtain gene ID from line: {line.strip()}")
    return m.group(1)


def get_gene_id_bed(transcript_name):
    m = BED_TRANSCRIPT_SUFFIX_RE.match(transcript_name)
    if m:
        return m.group(1)
    return transcript_name


def write_promoter(out, chrom, start, end, strand, gene_id, upstream, downstream):
    if strand == "+":
        left = start - upstream
        right = start + downstream
    elif strand == "-":
        left = end - downstream
        right = end + upstream
    else:
        sys.exit(f"Strand not specified for gene {gene_id}")

    if left < 1:
        left = 1

    out.write("\t".join([chrom, str(left), str(right), gene_id, "0", strand]) + "\n")


def process_gtf(args):
    # Pass 1: detect whether "gene" rows exist at all
    has_gene_rows = False
    with open_maybe_gz(args.infile) as infh:
        for line in infh:
            if line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) >= 3 and fields[2] == "gene":
                has_gene_rows = True
                break

    n_genes = 0
    n_written = 0

    if has_gene_rows:
        with open_maybe_gz(args.infile) as infh, open(args.outfile, "w") as out:
            for line in infh:
                if line.startswith("#"):
                    continue
                fields = line.rstrip("\n").split("\t")
                if len(fields) < 9 or fields[2] != "gene":
                    continue

                n_genes += 1
                gene_id = get_gene_id_gtf(fields[8], line)
                chrom = format_chrom(fields[0], args.chr_prefix)
                write_promoter(out, chrom, int(fields[3]), int(fields[4]), fields[6],
                                gene_id, args.upstream, args.downstream)
                n_written += 1

        print(f"[GTF mode] Used 'gene' feature rows. Processed {n_genes} genes, "
              f"wrote {n_written} promoters to {args.outfile}", file=sys.stderr)
    else:
        genes = OrderedDict()  # gene_id -> [chrom, min_start, max_end, strand]
        with open_maybe_gz(args.infile) as infh:
            for line in infh:
                if line.startswith("#"):
                    continue
                fields = line.rstrip("\n").split("\t")
                if len(fields) < 9 or fields[2] != "transcript":
                    continue

                gene_id = get_gene_id_gtf(fields[8], line)
                chrom = format_chrom(fields[0], args.chr_prefix)
                start = int(fields[3])
                end = int(fields[4])
                strand = fields[6]

                if gene_id not in genes:
                    genes[gene_id] = [chrom, start, end, strand]
                else:
                    g = genes[gene_id]
                    g[1] = min(g[1], start)
                    g[2] = max(g[2], end)

        with open(args.outfile, "w") as out:
            for gene_id, (chrom, start, end, strand) in genes.items():
                n_genes += 1
                write_promoter(out, chrom, start, end, strand, gene_id, args.upstream, args.downstream)
                n_written += 1

        print(f"[GTF mode] No 'gene' rows found — aggregated gene extent from 'transcript' rows. "
              f"Processed {n_genes} genes, wrote {n_written} promoters to {args.outfile}", file=sys.stderr)


def process_bed(args):
    genes = OrderedDict()  # gene_id -> [chrom, min_start(1-based), max_end, strand]
    n_transcripts = 0

    with open_maybe_gz(args.infile) as infh:
        for line in infh:
            line = line.rstrip("\n")
            if not line or line.startswith(("#", "track", "browser")):
                continue
            fields = line.split("\t")
            if len(fields) < 6:
                continue

            n_transcripts += 1

            raw_chrom, chrom_start, chrom_end, name, _score, strand = fields[:6]
            chrom = format_chrom(raw_chrom, args.chr_prefix)
            start = int(chrom_start) + 1  # BED 0-based -> 1-based, matches GTF convention
            end = int(chrom_end)
            gene_id = get_gene_id_bed(name)

            if gene_id not in genes:
                genes[gene_id] = [chrom, start, end, strand]
            else:
                g = genes[gene_id]
                g[1] = min(g[1], start)
                g[2] = max(g[2], end)

    n_genes = 0
    n_written = 0
    with open(args.outfile, "w") as out:
        for gene_id, (chrom, start, end, strand) in genes.items():
            n_genes += 1
            write_promoter(out, chrom, start, end, strand, gene_id, args.upstream, args.downstream)
            n_written += 1

    print(f"[BED mode] Read {n_transcripts} transcript rows, aggregated into {n_genes} genes, "
          f"wrote {n_written} promoters to {args.outfile}", file=sys.stderr)


def main():
    ap = argparse.ArgumentParser(description="Extract gene promoters from a GTF or BED12 file into a BED file.")
    ap.add_argument("-i", "--infile", required=True,
                     help="Input file: GTF (.gtf/.gtf.gz) or BED12 (.bed), e.g. iGenomes genes.gtf or genes.bed")
    ap.add_argument("-o", "--outfile", required=True, help="Output BED file")
    ap.add_argument("-u", "--upstream", type=int, default=500, help="Bases upstream of TSS (default: 500)")
    ap.add_argument("-d", "--downstream", type=int, default=500, help="Bases downstream of TSS (default: 500)")
    ap.add_argument("--chr-prefix", action="store_true",
                     help="Add 'chr' prefix to chromosome names (and rename Pt -> Cp)")
    ap.add_argument("--format", choices=["gtf", "bed"], default=None,
                     help="Force input format instead of auto-detecting from file extension")
    args = ap.parse_args()

    fmt = detect_format(args.infile, args.format)
    if fmt == "bed":
        process_bed(args)
    else:
        process_gtf(args)


if __name__ == "__main__":
    main()
