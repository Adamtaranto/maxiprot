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
# Make it executable if you saved it as ./maxiprot
chmod +x maxiprot

# Run on a miniprot GFF3 file
./maxiprot input.gff3 \
  --out-gff3 best_per_locus.gff3 \
  --out-tsv  best_per_locus.tsv
```

Read from **stdin**:

```bash
cat input.gff3 | ./maxiprot - \
  --out-gff3 best.gff3 --out-tsv best.tsv
```

---

## Installation

Requirements: Python 3.8+ and `pandas`.

```bash
python -m pip install pandas
# Then drop maxiprot somewhere on $PATH, e.g.:
install -m 0755 maxiprot /usr/local/bin/maxiprot
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
   Weighted geometric product (weights act as exponents). Rewards **balanced** candidates; any near‑zero factor drags the score down.

---

## Selection policies

`--selection-mode`:

- `best` *(default)*: best score (after gates)
- `prefer_intact`: pick best **intact** (fs=0 & st=0) if any exist, else best overall (pseudogenes still selectable)
- `longest`: longest AA; tie‑break on score
- `longest_prefer_intact`: longest among intact if available; else longest overall

> By default, pseudogenes are not dropped; they simply compete under your policy.

---

## Locus clustering

Hits are clustered into loci by **target** (`tname`) and **strand**. A new locus starts when the next hit begins more than `--locus-pad` nt after the running end of the current locus.

- `--locus-pad` (alias `--locus-gap`) default: `5000`

Tune this upward for fragmented assemblies or long introns; downward to split nearby tandem copies.

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

output:
  --id-prefix PBM
  --out-gff3 best_per_locus.gff3
  --out-tsv  best_per_locus.tsv
  --log-level INFO
```

---

## Example recipes

### 1) Balanced, prefer intact if possible

```bash
./maxiprot input.gff3 \
  --score-mode ms_cov_pos \
  --selection-mode prefer_intact \
  --min-cov 0.60 --min-pid 0.30 \
  --locus-pad 6000 \
  --out-gff3 best.gff3 --out-tsv best.tsv
```

### 2) Identity × coverage, but reward being long

```bash
./maxiprot input.gff3 \
  --score-mode pid_cov_len --w-pid 1 --w-cov 1 --w-len 1 --length-metric aa \
  --selection-mode best \
  --out-gff3 best.gff3 --out-tsv best.tsv
```

### 3) Pure longest, within your gates

```bash
./maxiprot input.gff3 \
  --selection-mode longest \
  --score-mode length --length-metric aa \
  --min-cov 0.60 --min-pid 0.30 \
  --out-gff3 longest.gff3 --out-tsv longest.tsv
```

### 4) Strict balance via geometric mean

```bash
./maxiprot input.gff3 \
  --score-mode geom --w-pid 1 --w-cov 1 --w-len 1 --length-metric frac \
  --selection-mode best \
  --out-gff3 best_geom.gff3 --out-tsv best_geom.tsv
```

### 5) Mirror miniprot’s own ranking

```bash
./maxiprot input.gff3 \
  --score-mode AS \
  --selection-mode best \
  --out-gff3 best_AS.gff3 --out-tsv best_AS.tsv
```

### 6) Prefer intact, but if none are intact, choose the longest

```bash
./maxiprot input.gff3 \
  --selection-mode longest_prefer_intact \
  --score-mode pid_cov_len --w-pid 1 --w-cov 1 --w-len 0.5 --length-metric frac \
  --out-gff3 best_pref_long.gff3 --out-tsv best_pref_long.tsv
```

---

## Logging

`--log-level` controls verbosity. For each locus, the log prints:

- number of candidates
- per candidate summary: **score**, **cov**, **pid**, **lenAA**, **positives**, **ms**, **AS**, **fs**, **st**, pass/fail, intact/pseudogene  
- the **selected** candidate and the global settings used

This is designed for traceability and to justify selections downstream.

---

## Performance tips

- miniprot GFFs can be large; `pandas` handles the PAF header extraction efficiently, and body feature parsing only inspects `mRNA`/`CDS` lines.
- If you only need the winner’s synthesized features (not original CDS parts), runtime is essentially linear in number of `##PAF` lines.

---

## Limitations & notes

- Requires `##PAF` header lines; if your miniprot run didn’t include them, re‑run with GFF output that preserves `PAF` headers.
- `positives` relies on miniprot’s `np:i` (AA positives) and aligned AA length.
- Synthesized features (when no matching mRNA/CDS exists) are minimal but standards‑compliant.

---

## Reproducibility

- Version string is embedded and printed at startup.  
- Pin `pandas` in your environment if you want byte‑for‑byte identical TSV ordering across machines.

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

- Li, H. (2023) Protein-to-genome alignment with miniprot. *Bioinformatics*, **39**, btad014 [[PMID: 36648328]][mp-pmid].

---

## Changelog (short)

- **0.4.0** — renamed to `maxiprot`, NumPy-docs, type hints, improved logging, flexible scoring/selection, correct multi‑line CDS IDs.
