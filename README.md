# Disease-regulatory flux map

A Nextflow pipeline that builds a **disease-regulatory flux map**: it splits each gene's
disease association into an own (cis) part and a network (trans) part, and attributes the
trans part to the specific upstream regulators that deliver it.

The pipeline is **data-agnostic** — it works on any cohort by changing configuration, not
code. This README describes the inputs you provide; pipeline stages are documented as they
are added.

---

## Inputs

The pipeline runs **per tissue** and needs three matched reference files. All are
tab-separated tables and may be gzipped (`.gz`). Column names and the delimiter are
configurable; the defaults are shown below.

### 1. Expression matrix — `expression.tsv[.gz]`

A gene-by-sample matrix of **analysis-ready** expression (already normalised and
covariate-adjusted — the pipeline uses the values as given). The first column is the gene
id; the remaining columns are the samples.

```
gene_id     S1       S2       S3
ENSG0001    0.42    -1.10     0.05
ENSG0002   -0.88     0.31     1.22
ENSG0003    0.10     0.07    -0.44
```

All genes are kept: every expressed gene can be a **target** in the regulatory network.

### 2. Genotype matrix — `genotype.tsv[.gz]`

A variant-by-sample matrix of **alternate-allele dosages in [0, 2]**. The leading
column(s) identify the variant (see *Variant identifiers* below); the remaining columns are
the samples — the **same samples** as in the expression matrix.

```
rsID      S1     S2     S3
rs1001    0.0    1.0    2.0
rs1002    1.0    0.0    1.0
```

### 3. eQTL map — `eqtl.tsv[.gz]`

A long table of **cis-eQTL associations** linking variants to genes — one row per
gene–variant pair, with effect size, standard error and p-value. This is what marks a gene
as a possible **regulator** (only genes with a cis-eQTL can regulate others).

```
gene_id     rsID      beta     se      pvalue
ENSG0001    rs1001    0.31    0.06    2.1e-7
ENSG0001    rs1002   -0.12    0.05    3.0e-3
ENSG0002    rs1002    0.27    0.06    8.0e-6
```

### Variant identifiers

The genotype and eQTL files must identify variants the **same way**:

- if **both** files have an `rsID` column, variants are matched by rsID;
- otherwise **both** must have `chromosome`, `position`, `ref`, `alt` columns, and variants
  are matched by `chr:pos:ref:alt`.

### How the files connect

| key | links |
|-----|-------|
| sample columns | expression ↔ genotype |
| `gene_id`      | expression ↔ eQTL |
| variant id     | genotype ↔ eQTL |

---

## Running the pipeline

Requires Nextflow ≥ 24.04.2 (and Java 11+).

Provide a samplesheet CSV with a header and one row per tissue:

```
tissue,expression,genotype,eqtl
AOR,/data/AOR.expr.tsv.gz,/data/AOR.geno.tsv.gz,/data/AOR.eqtl.tsv.gz
LIV,/data/LIV.expr.tsv.gz,/data/LIV.geno.tsv.gz,/data/LIV.eqtl.tsv.gz
```

Then run:

```
nextflow run flux.nf --samplesheet samples.csv --outdir results
```

Outputs are written per tissue under `results/harmonised/`. The column names and field
delimiter are set in `nextflow.config` (or overridden with `--<param>`), so a new dataset
needs no code changes.

---

## First step: harmonisation

The first stage aligns the inputs so the rest of the pipeline can rely on them
(`bin/harmonise_inputs.py`):

- keeps **all** expression genes (every gene is a potential network target);
- reduces the **genotype to the variants present in the eQTL map** (the only variants the
  method uses);
- restricts both expression and genotype to the **samples they share**, in one consistent
  order;
- writes a QC report of how many samples, genes and variants were kept.

It stops with a clear error if the files don't share samples, or if no eQTL pair survives
the matching.
