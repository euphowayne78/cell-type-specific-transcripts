# Cell-Type Specific Transcript Identification Pipeline

A computational pipeline for identifying **cell-type specific cytosolic transcripts** in human tissues, designed to support ADAR-based therapeutic mRNA technology.

## Background

Therapeutic mRNAs can be engineered with an upstream UAG stop codon in a sensor module. When this sensor base-pairs with a cell-type-specific endogenous transcript, endogenous ADAR edits the UAG → UIG, enabling **cell-type-restricted translation**. This pipeline identifies candidate trigger transcripts that are:

1. Highly and reliably expressed in the target cell type
2. Minimally expressed in all other cell types
3. Localized to the cytosol

---

## Project Structure

```
Specific_Transcript/
├── notebooks/                         # Jupyter notebooks (run in order)
│   ├── 01_data_acquisition.ipynb      # Download & cache pseudobulk data from CELLxGENE
│   ├── 02_preprocessing.ipynb         # QC, normalization, filtering
│   ├── 03_specificity_analysis.ipynb  # Compute Tau, fold-change, AUROC, detection rates
│   ├── 04_visualization.ipynb         # Heatmaps, dot plots, volcano plots
│   └── 05_validation.ipynb            # Cross-reference HPA, GTEx, literature
├── src/                               # Reusable Python modules
│   ├── data_loader.py                 # CELLxGENE Census API wrapper & pseudobulk aggregation
│   ├── preprocessing.py               # Normalization, QC, batch correction
│   ├── specificity.py                 # Tau, fold-change, AUROC, detection rate computation
│   ├── filters.py                     # Biotype, localization, annotation quality filters
│   ├── visualization.py               # Heatmap, dot plot, expression profile generators
│   └── query.py                       # End-user query interface logic
├── app/
│   └── streamlit_app.py               # Interactive Streamlit dashboard
├── data/
│   ├── raw/                           # Cached census queries (parquet/h5ad)
│   ├── processed/                     # Normalized pseudobulk matrices
│   └── references/                    # Housekeeping gene lists, lncATLAS data, etc.
├── results/                           # Output tables, heatmaps, ranked transcript lists
├── tests/                             # Unit tests for src/ modules
├── requirements.txt
└── README.md
```

---

## Quick Start

### 1. Set up the environment

```bash
python -m venv venv
source venv/bin/activate        # Windows: venv\Scripts\activate
pip install -r requirements.txt
```

### 2. Run the notebooks in order

Launch Jupyter and run each notebook sequentially:

```bash
jupyter notebook
```

| Notebook | Description | Est. Runtime |
|---|---|---|
| `01_data_acquisition.ipynb` | Downloads pseudobulk data from CELLxGENE Census | 20–60 min (first run only) |
| `02_preprocessing.ipynb` | QC, normalization, and filtering | < 5 min |
| `03_specificity_analysis.ipynb` | Specificity scoring for all genes and cell types | < 5 min |
| `04_visualization.ipynb` | Heatmaps, dot plots, expression profiles | < 5 min |
| `05_validation.ipynb` | Cross-validation against HPA, GTEx, known markers | Variable |

Results are cached to `data/processed/` after the first run — subsequent runs load from cache and are near-instant.

### 3. Launch the Streamlit dashboard (optional)

```bash
streamlit run app/streamlit_app.py
```

---

## Data Source

### Primary: CZ CELLxGENE Census
The most curated, standardized aggregation of human scRNA-seq data, covering 50M+ cells across 600+ cell type annotations. Accessed via the `cellxgene-census` Python API — no full download required.

**Filters applied at query time:**
- Adult human primary tissue only
- Healthy / normal donors (disease == 'normal')
- Minimum 500 cells per cell type
- Minimum 2 independent datasets per cell type
- Excludes: fetal/embryonic tissue, cell lines, organoids, tumor tissue, activated immune cells

### Supplementary (validation)
- **Human Protein Atlas (HPA)** — protein-level validation
- **GTEx** — bulk tissue expression cross-validation
- **lncATLAS** — cytoplasmic localization of lncRNAs
- **UniProt** — annotation quality tiers

---

## Specificity Metrics

All metrics are computed on a **pseudobulk matrix** (~25,000 genes × 100+ cell types). Pseudobulk means that instead of working with individual cells, each cell type is represented by a single aggregated profile — the mean CPM (counts per million) across all cells of that type. This reduces noise from single-cell dropout and makes computation tractable on a laptop.

---

### Primary: Tau Specificity Index

Tau is a well-validated summary statistic (Yanai et al., 2005) that captures how concentrated a gene's expression is across cell types, in a single number between 0 and 1.

**Formula:**

```
τ = Σᵢ(1 - x̂ᵢ) / (n - 1)
```

**Variables defined precisely:**

- **`i`** — index over cell types. If you have 100 cell types, `i` runs from 1 to 100.
- **`xᵢ`** — the mean CPM expression of the gene **in cell type `i`**. This is a single number per cell type — the average expression across all cells of that type (pseudobulk).
- **`max(x)`** — the maximum of `xᵢ` across **all** cell types for this gene. This is *not* the max within one cell type's cells; it is asking: *which cell type has the highest mean expression of this gene?* That highest value becomes the denominator for normalisation.
- **`x̂ᵢ = xᵢ / max(x)`** — the normalised expression of the gene in cell type `i`. Because you divide by the maximum, `x̂ᵢ` always lies between 0 and 1. The cell type with the highest expression gets `x̂ᵢ = 1`; all others get a value between 0 and 1.
- **`(1 - x̂ᵢ)`** — how far *below* the maximum is this cell type? A cell type with zero expression contributes 1; the top cell type contributes 0.
- **`n`** — total number of cell types included in the analysis.
- **`n - 1`** — normalisation factor so that the sum cannot exceed 1.

**Worked example with 5 cell types:**

Suppose gene ALB has these mean CPM values:

| Cell type | Mean CPM (xᵢ) | x̂ᵢ = xᵢ / 800 | (1 - x̂ᵢ) |
|---|---|---|---|
| Hepatocyte | 800 | 1.00 | 0.00 |
| T cell | 2 | 0.003 | 0.997 |
| Neuron | 0 | 0.000 | 1.000 |
| Cardiomyocyte | 1 | 0.001 | 0.999 |
| Keratinocyte | 3 | 0.004 | 0.996 |

```
τ = (0.00 + 0.997 + 1.000 + 0.999 + 0.996) / (5 - 1)
  = 3.992 / 4
  = 0.998
```

ALB is nearly perfectly specific to hepatocytes — τ ≈ 1.0.

Now consider a housekeeping gene like ACTB expressed at similar levels everywhere:

| Cell type | Mean CPM (xᵢ) | x̂ᵢ | (1 - x̂ᵢ) |
|---|---|---|---|
| Hepatocyte | 500 | 1.00 | 0.00 |
| T cell | 490 | 0.98 | 0.02 |
| Neuron | 510 | 1.00* | — |
| Cardiomyocyte | 495 | 0.99 | 0.01 |
| Keratinocyte | 505 | 1.00* | — |

*(max = 510, from neuron)*

```
τ = (0.02 + 0.04 + 0.00 + 0.03 + 0.01) / 4 ≈ 0.025
```

ACTB is essentially ubiquitous — τ ≈ 0.

**Interpretation:**
- **τ = 0** — gene is expressed equally in all cell types (e.g. GAPDH, ACTB)
- **τ = 1** — gene is expressed exclusively in one cell type
- **τ > 0.85** — high specificity, used as the primary filter in this pipeline
- **τ > 0.90** — highest confidence; use when false-positive risk must be minimised

**Why log-transform before computing Tau:**
Tau is computed on **log2(1 + mean CPM)** values, not raw CPM. Without log-transformation, a single extreme outlier cell type — say, a gene expressed at 10,000 CPM in one type and 10 CPM in all others — would produce a near-perfect τ ≈ 1 even if the absolute off-target expression (10 CPM) is biologically non-negligible. Log-transformation compresses the dynamic range, so τ reflects the *pattern* of specificity rather than being dominated by one extreme value.

Tau serves as a **pre-filter**: only genes passing the Tau threshold are evaluated for the supporting metrics, which makes the pipeline efficient.

---

### Supporting Metrics

Tau alone is necessary but not sufficient. A gene could have high Tau yet still be a poor therapeutic trigger — for example, if it is expressed in only 5% of target cells (unreliable), or at levels too low to drive ADAR editing. The five supporting metrics address these failure modes independently.

---

#### 1. On-Target Detection Rate — threshold > 25%

**What it measures:** The fraction of individual cells of the target type that express the gene (raw count > 0).

**Why it matters:** Single-cell RNA-seq has 10–30% technical dropout — transcripts present in a cell are not always captured. A gene detected in only 10% of target cells may be a technical artifact, not a reliable transcript. Requiring > 25% detection ensures the gene is robustly expressed across the population, not driven by a small outlier subpopulation.

**Why 25% specifically:** Going higher (e.g. > 50%) would discard legitimate tissue-resident markers with moderate but consistent expression. Going lower would admit noisy, dropout-dominated signals. 25% is a pragmatic midpoint supported by the scRNA-seq literature.

---

#### 2. On-Target Mean Expression — threshold > 1.0 log-normalized CPM

**What it measures:** The average expression level across all cells of the target type, in log2(1 + CPM) units.

**Why it matters:** ADAR editing efficiency is concentration-dependent. The therapeutic sensor RNA must compete with RNA secondary structure and RNA-binding proteins for access to the trigger transcript. A transcript expressed at very low levels will not reliably drive sufficient ADAR editing to activate translation. The > 1.0 log-CPM threshold corresponds roughly to > 1 CPM on the linear scale — a level associated with biologically meaningful, consistently detectable expression in bulk and single-cell data.

---

#### 3. Off-Target Detection Rate — threshold < 5% in any non-target cell type

**What it measures:** The maximum detection rate of the gene across all cell types *other* than the target. The worst-case off-target cell type must fall below 5%.

**Why it matters:** This is the **safety threshold**. Even if a transcript is highly expressed in the target cell, if it is also expressed in 20% of liver cells, 15% of kidney cells, or 10% of cardiomyocytes, the therapeutic mRNA would be translated in those tissues too — potentially causing toxicity. The < 5% threshold sets a strict ceiling on leaky off-target expression.

**Note:** For safety-critical therapeutic applications, this threshold should be tightened to < 1%. The 5% starting point is appropriate for initial candidate discovery; hits should be re-evaluated with stricter thresholds before advancing to experimental validation.

---

#### 4. Fold-Change vs. Next Highest Cell Type — threshold > 10x

**What it measures:** The ratio of mean CPM in the target cell type to the mean CPM in the **single highest-expressing non-target cell type** — i.e. how much more is the gene expressed in the target compared to its closest competitor.

**Formula:**

```
FC = (CPM_target + 0.1) / (CPM_next_highest + 0.1)
```

- **`CPM_target`** — mean CPM of the gene in the target cell type
- **`CPM_next_highest`** — mean CPM in whichever non-target cell type expresses it the most (the "worst-case competitor")
- **`0.1`** — pseudocount added to both sides to prevent division by zero and to dampen artificially inflated fold-changes when the denominator is near-zero (e.g. 0.001 CPM)

**Why the next-highest, not the average?** The therapeutic risk comes from the single most leaky cell type. If expression in even one off-target type is high, the therapeutic mRNA could be activated there. Using the maximum off-target value is the conservative, safety-first choice.

**Worked example:**

Suppose a candidate gene has these mean CPM values across 5 cell types:

| Cell type | Mean CPM | Role |
|---|---|---|
| Hepatocyte | 500 | Target |
| Cholangiocyte | 45 | Non-target (next highest) |
| T cell | 1.2 | Non-target |
| Neuron | 0.3 | Non-target |
| Keratinocyte | 0.8 | Non-target |

The next-highest non-target is cholangiocyte at 45 CPM.

```
FC = (500 + 0.1) / (45 + 0.1)
   = 500.1 / 45.1
   ≈ 11.1x   ✓  PASSES (> 10x)
```

Now suppose a second gene looks similar but has higher leakage into cholangiocytes:

| Cell type | Mean CPM | Role |
|---|---|---|
| Hepatocyte | 500 | Target |
| Cholangiocyte | 170 | Non-target (next highest) |
| T cell | 1.0 | Non-target |
| Neuron | 0.5 | Non-target |
| Keratinocyte | 0.6 | Non-target |

```
FC = (500 + 0.1) / (170 + 0.1)
   = 500.1 / 170.1
   ≈ 2.9x   ✗  FAILS (< 10x)
```

This gene might still have a high Tau (hepatocyte is clearly dominant), but the cholangiocyte signal is too close — a therapeutic mRNA using this trigger would risk being translated in cholangiocytes as well.

**Why the pseudocount matters:** Without it, a gene expressed at 500 CPM in the target and 0.001 CPM in the next-highest would yield FC = 500,000x — a meaningless number driven by sequencing noise at near-zero expression. The pseudocount of 0.1 caps this: FC = 500.1 / 0.101 ≈ 4,950x, still very high but no longer dominated by noise. In practice, any gene with a next-highest value below ~1 CPM is already well-covered by the off-target detection rate filter.

**Why 10x specifically:** ADAR editing rate in a given cell is approximately proportional to the local concentration of the trigger transcript (assuming the sensor mRNA is delivered at the same dose to all cell types). A 10x fold-difference in transcript abundance therefore translates to approximately a 10x difference in editing rate between the target cell and the worst-case off-target cell. This gives a practical therapeutic window: editing in the target cell is an order of magnitude more efficient than in any other cell type.

Note that sensor *occupancy* (how much of the time the sensor is actually duplexed with a trigger transcript) depends on the **absolute** concentration of the trigger relative to the duplex Kd — not on the fold-difference between cell types. Two cell types could have a 10x fold-difference in transcript levels yet both sit far below the Kd, meaning the sensor is mostly unbound in both. The fold-change threshold is therefore a proxy for relative editing efficiency, not a direct measure of occupancy or editing probability. Absolute expression level is addressed separately by the mean expression threshold (metric 2).

The 5–10x minimum is an empirical starting point drawn from published ADAR sensor literature; 10x is used here as a conservative threshold. Higher fold-change is always preferred.

**Biological basis:** ADAR editing efficiency scales with the local concentration of the guide RNA–target RNA duplex. Literature models of ADAR-based logic gates suggest a minimum 5–10x expression differential is required for reliable on/off switching. 10x is used here as a conservative starting point.

---

#### 5. AUROC (One-vs-Rest) — threshold > 0.9

**What it measures:** The Area Under the Receiver Operating Characteristic curve, treating each gene's expression profile as a binary classifier: can the gene's pseudobulk expression value alone correctly distinguish the target cell type from all other cell types?

**How it is computed:** For each (gene, cell type) pair, a one-vs-rest classification problem is constructed:
- The **n** cell-type mean CPM values serve as classifier scores
- The target cell type is the single **positive class** (label = 1)
- All other n−1 cell types are the **negative class** (label = 0)
- The ROC curve is drawn by sweeping a threshold from high to low CPM, at each step recording the True Positive Rate (TPR) and False Positive Rate (FPR)
- AUROC is the area under that curve

In this pseudobulk setting with exactly **one positive sample**, AUROC has a simple intuitive meaning:

> **AUROC = the fraction of non-target cell types that the target cell type outranks in expression.**

If the target cell type has the highest CPM of all cell types → AUROC = 1.0.
If one non-target cell type has higher CPM → AUROC = (n−2)/(n−1).

**Worked example with 5 cell types:**

*Gene A — hepatocyte-specific, target = Hepatocyte:*

| Rank | Cell type | Mean CPM | Class |
|---|---|---|---|
| 1 | Hepatocyte | 500 | Positive ✓ |
| 2 | Keratinocyte | 3 | Negative |
| 3 | T cell | 2 | Negative |
| 4 | Cardiomyocyte | 1 | Negative |
| 5 | Neuron | 0 | Negative |

Sweeping the threshold from high to low:

| Threshold crossed | Cell type included | TPR | FPR | Point on ROC |
|---|---|---|---|---|
| > 500 | (none) | 0.00 | 0.00 | (0.00, 0.00) |
| ↓ 500 | Hepatocyte (positive) | 1.00 | 0.00 | (0.00, 1.00) |
| ↓ 3 | Keratinocyte (negative) | 1.00 | 0.25 | (0.25, 1.00) |
| ↓ 2 | T cell (negative) | 1.00 | 0.50 | (0.50, 1.00) |
| ↓ 1 | Cardiomyocyte (negative) | 1.00 | 0.75 | (0.75, 1.00) |
| ↓ 0 | Neuron (negative) | 1.00 | 1.00 | (1.00, 1.00) |

AUROC = area of the rectangle = 1.0 × 1.0 = **1.0** ✓ PASSES

The target (Hepatocyte) was ranked #1 — it was captured before any false positive.

---

*Gene B — same target, but T cell has higher expression than Hepatocyte:*

| Rank | Cell type | Mean CPM | Class |
|---|---|---|---|
| 1 | T cell | 600 | Negative |
| 2 | Hepatocyte | 500 | Positive ✓ |
| 3 | Keratinocyte | 3 | Negative |
| 4 | Cardiomyocyte | 1 | Negative |
| 5 | Neuron | 0 | Negative |

| Threshold crossed | Cell type included | TPR | FPR | Point on ROC |
|---|---|---|---|---|
| > 600 | (none) | 0.00 | 0.00 | (0.00, 0.00) |
| ↓ 600 | T cell (negative) | 0.00 | 0.25 | (0.25, 0.00) |
| ↓ 500 | Hepatocyte (positive) | 1.00 | 0.25 | (0.25, 1.00) |
| ↓ 3 | Keratinocyte (negative) | 1.00 | 0.50 | (0.50, 1.00) |
| ... | ... | ... | ... | ... |

AUROC = (0.25 × 0) + (0.75 × 1.0) = **0.75** ✗ FAILS (< 0.9)

One false positive (T cell) was captured before the positive (Hepatocyte) — the classifier had to "pay" a 25% FPR before reaching TPR = 1.

---

**What AUROC catches that fold-change misses:**

With 100 cell types, the AUROC > 0.9 threshold means the target must be ranked in the **top 10% of all cell types** by expression. Fold-change only checks the single worst competitor; AUROC checks the target's rank across all cell types simultaneously. A gene could pass fold-change (10x over the next-highest) yet have AUROC < 0.9 if many other cell types cluster just below the next-highest, making the gene a poor global discriminator.

**Caveat:** Because this AUROC is computed on pseudobulk (one CPM value per cell type, not per individual cell), it is an approximation of the true per-cell AUROC. With only one positive sample, the metric reduces to a ranking statistic. It is used here as a consistency check alongside the other metrics, not as the primary criterion.

---

### Combined Scoring

A transcript is reported as a **high-confidence hit** (`all_pass = True`) only if it passes **all five** supporting metrics in addition to the Tau threshold. The full results table ranks candidates by:

1. `all_pass` (True first)
2. Tau (descending)
3. Fold-change (descending)

This ranking surfaces the most specific, highest-expressing, best-separated candidates at the top of each cell-type report.

---

## Output Files

After running all notebooks, the following files are generated:

| File | Description |
|---|---|
| `data/processed/pseudobulk_mean_cpm.parquet` | Mean CPM per gene per cell type |
| `data/processed/pseudobulk_detection_rate.parquet` | Fraction of cells expressing each gene |
| `data/processed/gene_id_map.parquet` | Ensembl ID ↔ gene symbol mapping |
| `data/processed/cell_type_summary.parquet` | Cell type counts and dataset representation |
| `results/specificity_scores.csv` | Full ranked table of all gene–cell type pairs |
| `results/heatmaps/` | Expression heatmaps per cell type |
| `results/cell_type_reports/` | Per-cell-type ranked transcript lists |

---

## Transcript Filtering

### Biotype priority
1. **Protein-coding mRNAs** — primary focus, cytosolic by default
2. **lncRNAs** — included only if cytoplasmic localization confirmed via lncATLAS (RCI > 0)
3. **Excluded:** rRNA, tRNA, snRNA, snoRNA, miRNA, mitochondrial transcripts, pseudogenes

### Additional exclusions
- Housekeeping genes (Eisenberg & Levanon, 2013)
- Ribosomal protein genes
- Mitochondrially-encoded genes
- Sex chromosome genes (unless sex-specific analysis)
- Transcripts < 200 nt

### Annotation quality tiers (ranked, not hard-filtered)
| Tier | Description |
|---|---|
| Tier 1 | UniProt reviewed (Swiss-Prot), well-characterized protein |
| Tier 2 | UniProt unreviewed (TrEMBL), ENSEMBL annotated |
| Tier 3 | Novel/uncharacterized with consistent expression |

---

## Cell Types Covered

The pipeline queries 100+ cell types spanning all major tissue categories:

- **Immune:** T cells (CD4+/CD8+ naive/memory/effector, Treg), B cells, NK cells, monocytes, macrophages, dendritic cells, neutrophils, mast cells
- **Epithelial:** alveolar type I/II, enterocytes, goblet cells, hepatocytes, renal tubular cells, keratinocytes, thyrocytes
- **Stromal:** fibroblasts, pericytes, smooth muscle, adipocytes, chondrocytes, osteoblasts
- **Endothelial:** arterial, venous, capillary, lymphatic, sinusoidal
- **Neural:** excitatory/inhibitory neurons, astrocytes, oligodendrocytes, microglia, Schwann cells
- **Muscle:** cardiomyocytes, skeletal myocytes
- **Other:** erythrocytes, megakaryocytes, podocytes, beta cells, acinar cells, melanocytes

A minimum of 100–150 cell types is used to avoid false positives from an incomplete competitor set.

---

## Sanity Checks

The pipeline validates itself against known marker genes:

| Gene | Expected cell type |
|---|---|
| ALB | Hepatocyte |
| CD3D | T cell |
| SNAP25 | Neuron |
| TNNT2 | Cardiomyocyte |
| KRT14 | Keratinocyte |
| CD79A | B cell |
| PECAM1 | Endothelial cell |

---

## Computational Requirements

The pseudobulk approach makes this **laptop-friendly**:

- Pseudobulk matrix: ~25,000 genes × 150 cell types = ~3.75M values (trivially small)
- All specificity computations run on the aggregated matrix in seconds
- Raw single-cell data is streamed one cell type at a time and discarded after aggregation
- Results are cached as parquet files to avoid repeated API calls

For environments without sufficient local resources, notebooks are Google Colab compatible.

---

## Running Tests

```bash
pytest tests/
```

---

## Dependencies

See [requirements.txt](requirements.txt). Key packages:

- `cellxgene-census` — primary data access
- `scanpy` / `anndata` — single-cell analysis framework
- `pandas`, `numpy`, `scipy`, `scikit-learn` — data processing and statistics
- `matplotlib`, `seaborn`, `plotly` — visualization
- `streamlit` — interactive dashboard

Python 3.10+ required.

---

## Development Phases

- [x] **Phase 1:** Data acquisition & preprocessing
- [x] **Phase 2:** Specificity analysis (Tau, fold-change, AUROC, detection rate)
- [x] **Phase 3:** Visualization & validation
- [x] **Phase 4:** Streamlit query interface
- [ ] **Phase 5 (planned):** ADAR sensor design integration — identify optimal target regions, predict RNA secondary structure, score editing efficiency, output sensor sequences

---

## References

- Yanai et al. (2005). Genome-wide midrange transcription profiles reveal expression level relationships in human tissue specification. *Bioinformatics*.
- Kryuchkova-Mostacci & Robinson-Rechavi (2017). A benchmark of gene expression tissue-specificity metrics. *Briefings in Bioinformatics*.
- Eisenberg & Levanon (2013). Human housekeeping genes, revisited. *Trends in Genetics*.
