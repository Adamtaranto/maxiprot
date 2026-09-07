# maxiprot

Pick the **best miniprot alignment per locus** from a miniprot GFF3 (with `##PAF` headers), score candidates with flexible, biology‑aware metrics, and emit a **standards‑compliant GFF3** (gene → mRNA → multi‑line CDS with a shared `ID`) plus a **TSV summary**.

- **Input**: miniprot GFF3 that includes `##PAF` header lines (the “PAF‑in‑GFF” lines that carry `AS:i`, `ms:i`, `np:i`, `fs:i`, `st:i`, `cg:Z`, `cs:Z`).
- **Output**:
  1) best‑per‑locus **GFF3** with valid hierarchy and shared CDS `ID` across exons;
  2) **TSV** summary (one row per winning candidate)
- **Selection**: longest or best‑scoring, with policies to **prefer intact** (no frameshift/stop) where available, while **never discarding pseudogenes** unless you ask to.
- **Locus clustering**: group non‑overlapping hits on the **same chromosome & strand** into loci using a tunable max gap (nt).

---

## Why this exists

When you map a family of orthologous proteins to a new genome with **miniprot**, multiple references often hit the same locus—sometimes with different extents (domain differences), frameshifts, or internal stops. `maxiprot` consolidates these alternatives, **scores** them using interpretable metrics (identity, coverage, length, miniprot scores), and outputs a **single, best** annotation per locus—while keeping pseudogenes discoverable.

---

## Quick start

```bash
# Sequence to annotate
GENOME="genome.fa"

# Multifasta set of orthologous proteins
REF_PROT="refs.faa"

# Map reference proteins to genome with miniprot
# Require 80% sequence identity and 60% sequence coverage
miniprot -t 4 --gff -N 50 --outs 0.80 --outc 0.60 $GENOME $REF_PROT > miniprot.gff3

# Run maxiprot on miniprot GFF3 file; outputs to files.
maxiprot filter miniprot.gff3 --out-gff3 best_per_locus.gff3 --out-tsv best.tsv --strict

# Extract translated CDS for all non-pseudogenes
maxiprot extract best_per_locus.gff3 -g $GENOME --extract protein --exclude-pseudogenes > best_intact.faa

```

Read from **stdin**:

```bash
miniprot --gff $GENOME $REF_PROT | maxiprot filter - \
  --out-gff3 best.gff3 --out-tsv best.tsv
```

Write gff on **stdout**:

```bash
maxiprot filter miniprot.gff3 1> best_per_locus.gff3
```

---

## Installation

Requirements: Python 3.10+, biopython, numpy, pandas, pyfaidx and rich.

Installation options:

1) Create a Conda environment with the latest version of `maxiprot` and its dependencies using the `environment.yml` config file in this repo.

```bash
# Create maxiprot env
conda env create -f environment.yml
# Activate new env
conda activate maxiprot
```

2) `pip install` the latest development version directly from this repo.

```bash
python -m pip install git+https://github.com/Adamtaranto/maxiprot.git
```

3) Clone and install an editable version of `maxiprot` if you want to edit the source code.

```bash
git clone https://github.com/Adamtaranto/maxiprot.git && cd maxiprot && pip install -e '.[dev]'
```

---

## Inputs and assumptions

- The GFF3 **must** contain miniprot’s `##PAF` header lines; they provide per‑alignment fields used in scoring (e.g., `AS:i`, `ms:i`, `np:i`, `fs:i`, `st:i`, `cg:Z`, `cs:Z`).
- Coverage gating is computed on the **reference (query) protein**: `(qend - qstart) / qlen`.
- Identity is computed from `cs:Z:` blocks as **identical AAs ÷ aligned AAs**.
- Output GFF3 follows the [GFF3 specification](https://github.com/The-Sequence-Ontology/Specifications/blob/master/gff3.md): multi‑line CDS **share the same `ID`**, and each CDS belongs to an `mRNA`, which belongs to a `gene`.

---

## Outputs

### 1) GFF3 (best per locus)

- One **gene** per locus, one **mRNA** child, and **CDS** parts.
- Every line of a multi‑exon CDS **shares the same `ID=`** (GFF3 requirement).
- If the original mRNA/CDS features are present in the input for the winner, they’re reused; otherwise a minimal hierarchy is synthesized from the PAF coordinates.

### 2) TSV (one row per winner)

Columns include:

- locus ID (`chrom:strand:index`), `tname`, `tstart`, `tend`, `strand`
- reference `qname`, `qlen`
- `cov_aa`, `pid_aa`, `positives`, `ms`, `AS`, `score_raw`
- `cds_aa_len`, frameshifts `fs`, in‑frame stops `st`, `status` (`intact`|`pseudogene`)
- `mrna_id`, `gff_cds_nt_len`, `passes` (gates), `mapq`

---

## Scoring modes

`--score-mode` controls how candidates are ranked **within a locus** (after gates). You can also **gate** with `--min-cov` (default 0.60) and `--min-pid` (default 0.30).

Gates always **exclude** failing candidates from selection. If no candidate at a locus passes, the locus is dropped under `--strict`; otherwise its single best-scoring candidate is emitted as a fallback, flagged `passes=False` in the TSV.

Common ingredients:

- **`pid_aa`**: identical AA fraction from `cs:Z:`
- **`cov_aa`**: query coverage `(qend - qstart) / qlen`
- **`positives`**: similar AA fraction from `np:i / aligned_AA`
- **`ms_per_qlen`**: `ms:i / qlen`
- **`AS_per_qlen`**: `AS:i / qlen`
- **Length metric** (`--length-metric frac|aa`): `len_frac = aligned_AA/qlen` or `cds_aa_len` (absolute AA)

### Modes

1. **`ms_cov_pos`** *(default; robust “all‑rounder”)*
   `ms_per_qlen × cov_aa × (0.5 + 0.5×positives)`
   Balances raw alignment quality, coverage, and substitution quality.

2. **`AS`**
   `AS_per_qlen` (miniprot’s full score per qlen). Mirrors mapper ranking.

3. **`ms`**
   `ms_per_qlen` (match score per qlen). Slightly more optimistic than `AS`.

4. **`pid_cov`**
   `pid_aa × cov_aa` (transparent, biology‑first).

5. **`pid_cov_len`**
   `pid_aa^w_pid × cov_aa^w_cov × length^w_len`
   Choose length via `--length-metric`. Set weights with `--w-pid --w-cov --w-len`.

6. **`length`**
   Just the length metric (fractional or absolute).

7. **`linear`**
   Weighted sum:
   `w_pid×pid + w_cov×cov + w_pos×positives + w_ms×ms_per_qlen + w_AS×AS_per_qlen + w_len×length`

8. **`geom`**
   Weighted geometric product (weights act as exponents) over the ingredients with **positive weights**. Rewards **balanced** candidates; any near‑zero factor drags the score down. Factors are clamped to ≥1e‑9 before exponentiation so zero/negative ingredients can't produce NaN; requesting `geom` with all weights zero is an error.

**Missing tags**: if a record lacks `cs:Z`, identity falls back to `nmatch/alen` and aligned length to `qend-qstart`; a missing `ms:i`/`AS:i` falls back along `AS → ms → nmatch` (warned once); a missing `np:i` makes `positives` 0. Note `len_frac` (gap-excluded aligned AA fraction) is deliberately distinct from `cov_aa` (span-based).

---

## Selection policies

`--selection-mode`:

- `best` *(default)*: best score (after gates)
- `prefer_intact`: pick best **intact** (fs=0 & st=0) if any exist, else best overall (pseudogenes still selectable)
- `longest`: longest AA; tie‑break on score
- `longest_prefer_intact`: longest among intact if available; else longest overall


> Note: By default, pseudogenes are not dropped; they simply compete under your policy.

---

## Locus clustering

Hits are clustered into loci by **target** (`tname`) and **strand**. A new locus starts when the next hit begins more than `--locus-pad` nt after the running end of the current locus. Locus labels have the form `tname:strand:index`.

- `--locus-pad` (alias `--locus-gap`) default: `5000`

Clustering is single-linkage chaining: a run of tandem gene copies each spaced closer than `--locus-pad` merges into **one** locus (and `--emit-mode best` then emits one gene for the whole array). Tune `--locus-pad` upward for fragmented assemblies or long introns; downward to split nearby tandem copies, or use `--emit-mode all_passing` to keep every passing copy.

---

## CLI reference (abridged)

```text
positional:
  gff                   miniprot GFF3 path or '-' for stdin (default: '-')

scoring:
  --score-mode {ms_cov_pos,AS,ms,pid_cov,pid_cov_len,length,linear,geom}
  --length-metric {frac,aa}           (default frac)
  --w-pid --w-cov --w-len             (default 1.0)
  --w-pos --w-ms --w-AS               (default 0.0 for linear/geom)

gates:
  --min-cov 0.60
  --min-pid 0.30
  --strict                           (drop loci with no passing candidates)

locus:
  --locus-pad 5000                   (alias: --locus-gap)

selection:
  --selection-mode {best,prefer_intact,longest,longest_prefer_intact}

emission:
  --emit-mode {best,all_passing}     (default best: one winner per locus)
  --max-per-locus N                  (cap for all_passing)

output:
  --id-prefix PBM
  --out-gff3 -                       (stdout by default)
  --out-tsv FILE                     (not written by default; '-' = stdout,
                                      only when --out-gff3 is a file)
  --log-level {DEBUG,INFO,WARNING,ERROR,CRITICAL}

```

---

## Example recipes

### 1) Identity × coverage, but reward being long

```bash
maxiprot filter input.gff3 \
  --score-mode pid_cov_len --w-pid 1 --w-cov 1 --w-len 1 --length-metric aa \
  --selection-mode best \
  --out-gff3 best.gff3 --out-tsv best.tsv
```

### 2) Balanced, prefer intact if possible

```bash
maxiprot filter input.gff3 \
  --score-mode ms_cov_pos \
  --selection-mode prefer_intact \
  --min-cov 0.60 --min-pid 0.80 \
  --locus-pad 6000 \
  --out-gff3 best.gff3 --out-tsv best.tsv
```

### 3) Pure longest, within your gates

```bash
maxiprot filter input.gff3 \
  --selection-mode longest \
  --score-mode length --length-metric aa \
  --min-cov 0.60 --min-pid 0.80 \
  --out-gff3 longest.gff3 --out-tsv longest.tsv
```

### 4) Strict balance via geometric mean

```bash
maxiprot filter input.gff3 \
  --score-mode geom --w-pid 1 --w-cov 1 --w-len 1 --length-metric frac \
  --selection-mode best \
  --out-gff3 best_geom.gff3 --out-tsv best_geom.tsv
```

### 5) Mirror miniprot’s own ranking

```bash
maxiprot filter input.gff3 \
  --score-mode AS \
  --selection-mode best \
  --out-gff3 best_AS.gff3 --out-tsv best_AS.tsv
```

### 6) Prefer intact, but if none are intact, choose the longest

```bash
maxiprot filter input.gff3 \
  --selection-mode longest_prefer_intact \
  --score-mode pid_cov_len --w-pid 1 --w-cov 1 --w-len 0.5 --length-metric frac \
  --out-gff3 best_pref_long.gff3 --out-tsv best_pref_long.tsv
```

---

## Protein extraction companion: `maxiprot extract`

Extract protein sequences from the **GFF3 output of** `maxiprot` by pulling and translating the annotated CDS features.

**Features**

- Input: GFF3 from `maxiprot` (file or stdin)
- Output: **FASTA to stdout** by default (or --out-faa PATH)
- Genome FASTA may be **plain or bgzip-compressed** (.bgz); .fai index auto-created with pyfaidx if missing
- Choose **NCBI translation table**: `--transl-table` (default `1`; all NCBI tables supported via Biopython)
- **Exclude pseudogenes**: --exclude-pseudogenes skips proteins with internal `*`
- **Max annotations per contig** filter: `--max-annos-per-contig N` (e.g., `1` to enforce exclusivity)
- Warnings for: CDS length not divisible by 3; non-ATGC bases in CDS

**Usage**

Basic pipeline (miniprot → maxiprot → proteins):

```bash
miniprot --gff genome.fa ref_proteins.faa  > miniprot.gff3

# Pick best per locus; capture GFF to file, TSV to file.
maxiprot filter miniprot.gff3 --out-gff3 best.gff3 --out-tsv best.tsv

# Extract proteins from the resulting GFF into FASTA
maxiprot extract best.gff3 -g genome.fa > best.faa
```

Streaming (no intermediate files):

```bash
miniprot --gff genome.fa ref_proteins.faa \
| maxiprot filter - \
| maxiprot extract -g genome.fa > best.faa
```

Exclude pseudogenes and use bacterial table (11):

```bash
maxiprot extract best.gff3 -g genome.fa \
  --exclude-pseudogenes --transl-table 11 > proteins.faa
```

Skip contigs with more than one annotation (keep only uniquely annotated contigs):

```bash
maxiprot extract best.gff3 -g genome.fa \
  --max-annos-per-contig 1 > proteins_unique.faa
```

**Notes**

If `*` appears **internally** in the translation and `--exclude-pseudogenes` is **not** set, the sequence is still output and `*` is preserved so downstream tools can detect pseudogenes.

---
## Logging

`--log-level` controls verbosity. For each locus, the log prints:

- number of candidates (DEBUG)
- the **selected** candidate per locus with score/cov/pid/length and pass status

Logs go to **stderr** only; data outputs (GFF3, TSV) never share a stream with logs.

---

## Limitations & notes

- Requires `##PAF` header lines; if your miniprot run didn’t include them, re‑run with GFF output that preserves `PAF` headers.
- `positives` relies on miniprot’s `np:i` (AA positives) and aligned AA length.

---

## FAQ

**Q: I want to keep pseudogenes in the output.**
A: That’s the default—pseudogenes are scored and can win. If you **prefer intact**, use `--selection-mode prefer_intact`.

**Q: My reference set mixes short/long variants (extra domains).**
A: Use fractional length (`--length-metric frac`) so long references don’t dominate by sheer size.

**Q: How do I merge far‑apart exons?**
A: Increase `--locus-pad`. It merges non‑overlapping hits on the same strand within that distance into one locus.

---

## Citation

When using the workflow described in this README please cite the `maxiprot` GitHub repo directly as well as `miniprot`.

- Li, H. (2023) Protein-to-genome alignment with miniprot. *Bioinformatics*, **39**, btad014 [[PMID: 36648328]][mp-pmid].

---

## Changelog

See [CHANGELOG.md](CHANGELOG.md).
