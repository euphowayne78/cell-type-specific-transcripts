# CLAUDE.md — Cell-Type Specific Transcript Identification Pipeline

## Project Overview

This project identifies **cell-type specific cytosolic transcripts** in human tissues for use as
triggers in ADAR-based therapeutic mRNA technology. The core application: engineer therapeutic mRNAs
with an upstream UAG stop codon in a sensor module that, upon base-pairing with a cell-type-specific
endogenous transcript, gets edited (UAG -> UIG) by endogenous ADAR, enabling cell-type-restricted
translation.

**Critical requirement:** identified transcripts must be (1) highly and reliably expressed in the
target cell type, (2) minimally expressed in all other cell types, and (3) localized to the cytosol.

## Data Source

### Primary: CZ CELLxGENE Census (https://cellxgene.cziscience.com/)
- The most curated, standardized aggregation of human scRNA-seq datasets
- Uses Cell Ontology (CL) for standardized cell type annotations
- Accessible via Python API: `cellxgene-census` (TileDB-backed, no full download needed)
- Supports pseudobulk queries (aggregate expression per cell type) — laptop-friendly
- Covers 600+ cell type annotations, 50M+ cells across hundreds of studies

### Supplementary (for validation/filtering)
- **Human Protein Atlas (HPA)** — protein-level validation of transcript specificity
- **GTEx** — bulk tissue-level expression for cross-validation
- **lncATLAS** — subcellular localization of lncRNAs (filter for cytoplasmic)
- **UniProt** — annotation quality tiers for protein-coding genes

### Data Filters (applied at download/query time)
- **Include:** primary adult human tissue, normal/healthy donors only
- **Exclude:**
  - Tumor / cancer cells and tissues
  - Fetal / embryonic tissues
  - Cell lines (immortalized, iPSC-derived unless validated)
  - Disease-state samples (inflammatory, infected, autoimmune)
  - Activated / stimulated immune cells (use resting-state only)
  - Organoids (unless no primary tissue alternative exists)

## Cell Type Strategy

### Minimum 100-150 cell types for comprehensive specificity assessment
A short list (10-20) will produce **false positives** — transcripts that appear specific only because
competitor cell types were omitted. Use the Cell Ontology hierarchy from CELLxGENE.

### Representative cell type categories (non-exhaustive):
- **Immune:** T cells (CD4+ naive, CD4+ memory, CD8+ naive, CD8+ effector, Treg, gamma-delta),
  B cells (naive, memory, germinal center), plasma cells, NK cells, monocytes (classical,
  non-classical), macrophages (alveolar, Kupffer, tissue-resident), dendritic cells (cDC1, cDC2,
  pDC), neutrophils, eosinophils, basophils, mast cells
- **Epithelial:** alveolar type I/II, intestinal enterocytes, goblet cells, Paneth cells, club cells,
  ciliated cells, keratinocytes (basal, suprabasal), hepatocytes, cholangiocytes, renal tubular
  (proximal, distal, collecting duct), urothelial cells, thyrocytes
- **Stromal:** fibroblasts (multiple tissue origins), myofibroblasts, pericytes, smooth muscle cells,
  adipocytes (white, brown), chondrocytes, osteoblasts, osteoclasts
- **Endothelial:** arterial, venous, capillary, lymphatic, sinusoidal (liver)
- **Neural:** neurons (excitatory, inhibitory subtypes), astrocytes, oligodendrocytes, OPCs,
  microglia, Schwann cells, satellite glial cells
- **Muscle:** cardiomyocytes, skeletal myocytes, smooth muscle
- **Other specialized:** erythrocytes, megakaryocytes, platelets, podocytes, mesangial cells,
  melanocytes, Sertoli cells, Leydig cells, beta cells (pancreatic), alpha cells, acinar cells,
  ductal cells (pancreatic)

### Cell type quality thresholds
- Minimum 500 cells per cell type in the census (for statistical robustness)
- Require representation from >= 2 independent studies per cell type (batch robustness)

## Specificity Metrics

### Primary: Tau Specificity Index
- Range: 0 (ubiquitous) to 1 (perfectly cell-type specific)
- **Threshold: Tau > 0.85** (stringent; use > 0.9 for highest-confidence hits)
- Well-validated in literature (Yanai et al., 2005; Kryuchkova-Mostacci & Robinson-Rechavi, 2017)
- Computed on pseudobulk mean expression across cell types

### Supporting metrics (all must pass):
| Metric | Threshold | Rationale |
|--------|-----------|-----------|
| On-target detection rate | > 25% of cells | Ensures transcript is reliably present, not driven by outlier cells |
| On-target mean expression | > 1.0 log-normalized CPM | Meaningful expression level for ADAR sensor pairing |
| Off-target detection rate | < 5% in ANY non-target type | Minimizes leaky therapeutic activation |
| Fold-change vs next highest | > 10x | Clear separation from nearest competitor cell type |
| AUROC (one-vs-rest) | > 0.9 | Discriminative power as a cell type classifier |

### Why these thresholds:
- **25% detection rate:** scRNA-seq has ~10-30% dropout; transcripts detected in >25% of cells are
  robustly expressed. Going higher (e.g., >50%) loses tissue-resident markers with moderate
  expression. Going lower picks up noise.
- **<5% off-target:** for a therapeutic, even 5% leaky expression may be too high — this is a
  starting point. For safety-critical applications, tighten to <1%.
- **10x fold-change:** ADAR editing efficiency is concentration-dependent. A 10x expression
  differential provides a reasonable therapeutic window. Literature suggests 5-10x is minimum
  for functional ADAR sensor discrimination.

## Transcript Filtering

### Biotype priority:
1. **Protein-coding mRNAs** — primary focus, cytosolic by default
2. **lncRNAs** — include ONLY if cytoplasmic localization confirmed via lncATLAS (RCI > 0)
3. **Exclude:** rRNA, tRNA, snRNA, snoRNA, miRNA, mitochondrial transcripts, pseudogenes

### Annotation quality tiers (do NOT hard-filter, rank instead):
- **Tier 1:** UniProt reviewed (Swiss-Prot), well-characterized protein, known function
- **Tier 2:** UniProt unreviewed (TrEMBL), predicted protein, ENSEMBL annotated
- **Tier 3:** Novel/uncharacterized transcripts with consistent expression patterns
- Present all tiers in results; let the user decide based on risk tolerance

### Additional filters:
- Exclude housekeeping genes (use list from Eisenberg & Levanon, 2013)
- Exclude ribosomal protein genes
- Exclude mitochondrially-encoded genes
- Exclude genes on sex chromosomes (unless sex-specific analysis desired)
- Minimum transcript length > 200 nt (practical for ADAR sensor design)

## Project Architecture

```
Specific_Transcript/
├── CLAUDE.md                          # This file
├── requirements.txt                   # Python dependencies
├── .claude/
│   └── launch.json                    # Dev server config (Streamlit)
├── notebooks/                         # Jupyter notebooks for development/exploration
│   ├── 01_data_acquisition.ipynb      # Download and cache pseudobulk data from CELLxGENE
│   ├── 02_preprocessing.ipynb         # QC, normalization, filtering
│   ├── 03_specificity_analysis.ipynb  # Compute Tau, fold-change, AUROC, detection rates
│   ├── 04_visualization.ipynb         # Heatmaps, dot plots, volcano plots
│   └── 05_validation.ipynb            # Cross-reference HPA, GTEx, literature
├── src/                               # Reusable modules
│   ├── __init__.py
│   ├── data_loader.py                 # CELLxGENE Census API wrapper
│   ├── preprocessing.py               # Normalization, QC, batch correction
│   ├── specificity.py                 # Tau, fold-change, AUROC, detection rate computation
│   ├── filters.py                     # Biotype, localization, annotation quality filters
│   ├── visualization.py               # Heatmap, dot plot, expression profile generators
│   └── query.py                       # End-user query interface logic
├── app/                               # Streamlit dashboard for end-user querying
│   └── streamlit_app.py
├── data/
│   ├── raw/                           # Cached census queries (parquet/h5ad)
│   ├── processed/                     # Normalized pseudobulk matrices
│   └── references/                    # Housekeeping gene lists, lncATLAS data, etc.
├── results/                           # Output tables, heatmaps, ranked transcript lists
│   ├── specificity_scores.csv
│   ├── heatmaps/
│   └── cell_type_reports/
└── tests/                             # Unit tests for src/ modules
    ├── test_specificity.py
    └── test_filters.py
```

## Interface Strategy

### Development phase: Jupyter Notebooks
- Use notebooks (01-05) for iterative development, exploration, visualization
- Each notebook is self-contained with clear inputs/outputs
- Notebooks call into `src/` modules for reusable logic

### End-user phase: Streamlit Dashboard
- User selects a cell type of interest from dropdown
- Dashboard displays: ranked transcript list, expression heatmap, specificity scores
- Filters: Tau threshold slider, biotype checkboxes, annotation tier selection
- Export: CSV download of filtered results
- Preferred over Jupyter for non-technical end users (no code execution needed)

## Computational Strategy

### Laptop-friendly approach (primary):
- Use **pseudobulk aggregation** via CELLxGENE Census API — queries return cell-type-level
  summary statistics, NOT individual cell matrices
- Pseudobulk matrix size: ~25,000 genes x 150 cell types = ~3.75M values (trivially small)
- All specificity computations run on this aggregated matrix — seconds, not hours
- Cache query results locally as parquet files to avoid repeated API calls

### If raw single-cell access is needed (e.g., detection rate computation):
- Query cell-type subsets via Census API (stream, don't download all)
- Process one cell type at a time, compute statistics, discard raw matrix
- Use `tiledbsoma` lazy loading — only reads data into memory as needed
- Chunked processing: iterate over cell types in batches of 10-20

### Testing strategy for code correctness:
- Use a **small test subset**: 5 well-known cell types (hepatocyte, T cell, neuron, cardiomyocyte,
  keratinocyte) with known marker genes as ground truth
- Verify pipeline recovers known markers: ALB (hepatocyte), CD3D (T cell), SNAP25 (neuron),
  TNNT2 (cardiomyocyte), KRT14 (keratinocyte)
- Include synthetic data tests: create mock expression matrices with planted specific genes
- Run full pipeline on small subset first, validate, then scale

### Google Colab fallback:
- If laptop is insufficient, notebooks are Colab-compatible
- Mount Google Drive for persistent data/results storage
- Use Colab's free GPU/TPU runtime for any heavy computation
- Keep `src/` modules importable from Colab via sys.path or pip install from GitHub

## Key Python Dependencies

```
cellxgene-census          # CELLxGENE Census API (primary data access)
tiledbsoma                # Underlying data access layer for Census
scanpy                    # scRNA-seq analysis framework
anndata                   # Data structure for single-cell data
pandas                    # DataFrames for pseudobulk and results
numpy                     # Numerical computation
scipy                     # Statistical tests
scikit-learn              # AUROC computation, clustering
matplotlib                # Plotting
seaborn                   # Heatmaps, statistical visualization
plotly                    # Interactive plots (for Streamlit)
streamlit                 # End-user dashboard
pyarrow                   # Parquet file I/O
requests                  # API calls to HPA, UniProt
```

## Coding Conventions

- Python 3.10+
- Type hints on all function signatures in `src/`
- Docstrings (NumPy style) on all public functions in `src/`
- Notebooks: markdown cells explaining each step, rationale, and expected output
- Use `logging` module, not print statements, in `src/`
- Pandas operations: prefer vectorized over iterrows
- All file paths via `pathlib.Path`, not string concatenation
- Cache expensive API calls to `data/raw/` as parquet
- Results reproducibility: set random seeds, log Census version/date

## ADAR Sensor Design Considerations (Context for Transcript Selection)

- Sensor requires ~20-40 nt complementary region to the trigger transcript
- Trigger transcript must be **abundant in cytosol** at the time of therapeutic mRNA delivery
- Avoid transcripts with high secondary structure in the target region (impedes base-pairing)
- Avoid transcripts with known RNA-binding protein (RBP) occupancy at target site
- Prefer transcripts with long, accessible 3'UTR or CDS regions
- A-to-I editing by ADAR has sequence context preferences (5' neighbor: U>A>C>G)
- Final output should include transcript sequence + suggested sensor target regions (future phase)

## Development Phases

### Phase 1: Data Acquisition & Preprocessing
- Set up CELLxGENE Census API access
- Query pseudobulk expression for 100+ cell types, adult human, healthy tissue
- Normalize and QC the aggregated matrix
- Cache locally

### Phase 2: Specificity Analysis
- Compute Tau, fold-change, detection rate, AUROC for all genes across all cell types
- Apply biotype and localization filters
- Generate ranked transcript lists per cell type

### Phase 3: Visualization & Validation
- Expression heatmaps (genes x cell types)
- Dot plots (expression level + detection rate)
- Cross-validate top hits against HPA, GTEx, published marker gene lists
- Sanity check with known cell type markers

### Phase 4: Query Interface
- Build Streamlit app for end-user exploration
- Cell type selector, threshold controls, export functionality
- Interactive heatmap with filtering

### Phase 5 (Future): Sensor Design Integration
- For top transcript hits, identify optimal ADAR sensor target regions
- Predict RNA secondary structure (RNAfold) at candidate binding sites
- Score ADAR editing efficiency based on sequence context
- Output sensor sequences ready for experimental validation
