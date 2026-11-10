# Gut microbiome across MASLD, cirrhosis, and HCC in a Saudi cohort (16S rRNA, V4)

Analysis code for:

> Bahitham W. *et al.* Gut microbiome alterations across MASLD, cirrhosis, and hepatocellular carcinoma in a Saudi Arabian cohort. *Microbiology Spectrum* (under revision).

Cross-sectional study of 205 participants: healthy controls (n = 70), MASLD (n = 59), cirrhosis (n = 59), and HCC (n = 17). Fecal 16S rRNA gene sequencing (V4, 515F/806R; Illumina MiSeq 2 × 300 bp).

## Data availability

| Data | Where |
|---|---|
| Raw sequencing reads | NCBI SRA, BioProject **PRJNA1459468** |
| Participant-level clinical metadata | Not public (ethics approval KAIMRC NRJ22J/236/09; KAUST 24IBEC012). Available from the corresponding author on reasonable request. |
| Analysis code and summary outputs | This repository |

Because the clinical metadata cannot be shared, scripts that use clinical covariates will not run end-to-end without it. Column names and units expected by the scripts are listed in `metadata_template.tsv` (header and one example row; no participant data).

## Repository layout

The scripts use relative paths and expect two sibling folders, `16S_reprocess/` (data; not included) and `16S_liver_disease/` (code):

```
.
├── 16S_reprocess/                 # NOT included: QIIME 2 outputs + sample metadata
│   ├── qiime2_analysis_v2/        #   exported tables (alpha, beta, taxonomy, feature tables)
│   └── sample_metadata.tsv        #   clinical metadata (restricted)
├── qiime2/                        # QIIME 2 commands (primer removal → taxonomy → diversity)
└── 16S_liver_disease/
    ├── Figure_1/  … Figure_8/     # one folder per main figure
    ├── Supp_Figure_S1/ … S3/      # supplementary figures
    ├── Supp_Analyses/             # sensitivity analyses + alpha effect sizes
    ├── Main_Tables/               # Tables 1–2
    └── Supp_Tables/               # Supplementary Tables S1–S6
```

Each script sets its working directory to its own folder and writes outputs there.

## Sequence processing (QIIME 2 2024.10.1)

| Step | Settings |
|---|---|
| Primer removal | q2-cutadapt, anchored 5′ primers (515F `GTGYCAGCMGCCGCGGTAA`, 806R `GGACTACNVGGGTWTCTAAT`); untrimmed pairs discarded |
| Denoising | DADA2 1.30.0; truncation 230 (F) / 200 (R) bp; maxEE 2; consensus chimera removal |
| Filtering | ASV length 250–256 bp; ASVs without a phylum assignment and mitochondrial, chloroplast and eukaryotic ASVs removed (archaea retained) → 3,868 ASVs |
| Taxonomy | Naïve Bayes classifier (confidence 0.7), SILVA 138 (99%), human-stool-weighted |
| Phylogeny | MAFFT + FastTree |
| Rarefaction | 13,000 reads/sample (all 205 retained); sensitivity at 5,000–20,000 reads |

## Scripts → manuscript items

| Script | Folder | Produces |
|---|---|---|
| `Figure1_Study.R` | `Figure_1` | Figure 1 (study design) |
| `Figure2_Alpha_Diversity.R` | `Figure_2` | Figure 2 |
| `Alpha_Effect_Sizes.R` | `Supp_Analyses` | Alpha effect sizes, KW/BH, minimum detectable effect (Supplementary Table S2C) |
| `Figure3_BetaDiversity_FINAL.R` | `Figure_3` | Figure 3 (PCoA, PERMANOVA, PERMDISP) |
| `Figure4_Confounder_PERMANOVA.R` | `Figure_4` | Figure 4; Table 2 inputs |
| `Figure5_TaxonomicComposition_PhylumGenus.R` | `Figure_5` | Figure 5; F/B ratio with bootstrap CI |
| `Figure6_ANCOMBC2_Heatmap.R` | `Figure_6` | Figure 6 (ANCOM-BC2, unadjusted and adjusted) |
| `Figure7_Spearman_Heatmap.R` | `Figure_7` | Figure 7; Supplementary Figure S4 (partial correlations) |
| `Dysbiosis_Index_PCA.R` | `Figure_8` | Dysbiosis Index (CLR → PCA → distance to healthy centroid) |
| `Figure8_Dysbiosis_Index.R` | `Figure_8` | Figure 8; leave-one-out, split-half, within-group and medication checks |
| `FigS1_Rarefaction.R` | `Supp_Figure_S1` | Supplementary Figure S1 |
| `Supplementary_Figure_S2_DESeq2_Concordance.R` | `Supp_Figure_S2` | Supplementary Figure S2 (DESeq2 concordance) |
| `Supplementary_Figure_S3_Variance_Decomposition.R` | `Supp_Figure_S3` | Supplementary Figure S3 (variance decomposition, ω²) |
| `Supplementary_Analyses_Revision.R` | `Supp_Analyses` | Sensitivity analyses: W*d test, balanced subsampling, rarefaction depth, patients-only, viral-etiology exclusion, covariate-adjusted DA |
| `Generate_Main_Tables.R` | `Main_Tables` | Tables 1 and 2 |
| `Generate_Supplementary_Tables.R` | `Supp_Tables` | Supplementary Tables S1–S6 (+ S4 Excel workbook) |

### Run order

1. QIIME 2 commands in `qiime2/`
2. `Figure_2` … `Figure_8`, `Supp_Figure_S1`, `Supp_Figure_S2`
3. `Supp_Analyses/Supplementary_Analyses_Revision.R` and `Supp_Analyses/Alpha_Effect_Sizes.R`
4. `Supp_Figure_S3`
5. `Main_Tables`, then `Supp_Tables`

## Software

R 4.5.1; phyloseq 1.54.0; vegan 2.7.5; ANCOMBC 2.12.0; DESeq2 1.50.2; ggplot2 4.0.2. Random seeds are set in each script.

Multiple testing: Benjamini–Hochberg FDR throughout (q < 0.05).

Reporting follows the STORMS checklist (Supplementary Table S6).

## License

Code: MIT (see `LICENSE`).
