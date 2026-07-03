# gtf_to_promoters

Extract gene promoter regions from a GTF or BED12 file into a BED file, based on strand-aware windows around the transcription start site (TSS).

For each gene:

| Strand | Promoter window |
|--------|------------------|
| `+`    | `[start - upstream, start + downstream]` |
| `-`    | `[end - downstream, end + upstream]` |

By default, `upstream = 500` and `downstream = 500` (both configurable via `-u`/`-d`, see [Options](#options)).

Output coordinates are clamped to a minimum of `1`.

This is a Python reimplementation of a Perl script originally written for Ensembl-style GTFs, extended to also handle GTFs that lack `gene`-level rows and to support BED12 input directly.

## Why

Some GTFs — notably older Ensembl releases and several of the genome annotation files bundled in [iGenomes](https://ewels.github.io/AWS-iGenomes/) (e.g. TAIR10 for *Arabidopsis thaliana*) — don't include explicit `gene` feature rows, only `transcript`, `exon`, `CDS`, etc. This script detects that automatically and falls back to inferring gene boundaries from `transcript` rows. It can also work directly from a `genes.bed` (BED12) file, which is a common alternative format for the same annotation.

## Requirements

- Python 3.6+
- No external dependencies (standard library only)

## Usage

```bash
chmod +x gtf_to_promoters.py

# GTF input (plain or gzipped)
./gtf_to_promoters.py -i genes.gtf -o promoters.bed
./gtf_to_promoters.py -i genes.gtf.gz -o promoters.bed

# BED12 input (auto-detected from .bed extension)
./gtf_to_promoters.py -i genes.bed -o promoters.bed

# Custom promoter window
./gtf_to_promoters.py -i genes.gtf -u 1000 -d 200 -o promoters.bed

# Add 'chr' prefix to chromosome names (also renames Pt -> Cp for plastid)
./gtf_to_promoters.py -i genes.gtf --chr-prefix -o promoters.bed

# Force input format if the extension doesn't match (.txt, etc.)
./gtf_to_promoters.py -i genes.txt --format bed -o promoters.bed
```

### Options

| Flag | Description | Default |
|------|-------------|---------|
| `-i`, `--infile` | Input GTF (`.gtf`/`.gtf.gz`) or BED12 (`.bed`) file | required |
| `-o`, `--outfile` | Output promoters BED file | required |
| `-u`, `--upstream` | Bases upstream of TSS | `500` |
| `-d`, `--downstream` | Bases downstream of TSS | `500` |
| `--chr-prefix` | Add `chr` prefix to chromosome names; renames `Pt` → `Cp` | off |
| `--format {gtf,bed}` | Force input format instead of auto-detecting from extension | auto |

## How gene boundaries are determined

**GTF input:**
- If the file contains `gene` feature rows, those are used directly.
- Otherwise, gene extent is inferred by taking `min(start)` / `max(end)` across all `transcript` rows sharing the same `gene_id`.

**BED12 input** (e.g. an iGenomes `genes.bed`, one row per transcript):
- `gene_id` is derived from the transcript name (column 4) by stripping a trailing `.<number>` suffix (e.g. `AT1G01010.1` → `AT1G01010`).
- Gene extent is `min(chromStart)` / `max(chromEnd)` across all transcripts sharing that `gene_id`.
- BED's 0-based `chromStart` is converted to 1-based internally so output is consistent with the GTF path.

In both cases, the script prints a summary to stderr indicating which code path was used and how many genes/transcripts were processed — useful for a quick sanity check against `grep -oP 'gene_id "\K[^"]+' genes.gtf | sort -u | wc -l` or similar.

## Output format

Standard 6-column BED:

```
chrom   start   end   gene_id   score   strand
```

`score` is always `0` (placeholder, matches the original Perl script's convention).

## Example

```bash
$ ./gtf_to_promoters.py -i genes.gtf -o promoters.bed
[GTF mode] No 'gene' rows found — aggregated gene extent from 'transcript' rows. Processed 27655 genes, wrote 27655 promoters to promoters.bed

$ head -3 promoters.bed
1    3131    4131    AT1G01010    0    +
1    8237    9237    AT1G01020    0    -
1    11148   12148   AT1G01030    0    -
```

## Notes

- Chromosome naming is passed through as-is unless `--chr-prefix` is used. Check what's in your file before deciding:
  ```bash
  zcat -f genes.gtf | grep -v '^#' | cut -f1 | sort -u
  ```
- Coordinates follow the same (slightly non-standard) convention as the original Perl script: promoter start/end are computed directly from 1-based GTF-style coordinates rather than strictly BED half-open coordinates. This is consistent with how the original script's output has been used downstream, but worth knowing if you're feeding the output into strict BED-spec tools.

## License

MIT (or update to match your repo's license)
