"""
Random Forest classifier for functional gene vs. pseudogene prediction,
incorporating SNP enrichment features (1000 Genomes, gnomAD, pangenome,
dbSNP) AND promoter motif features (type 1/2/3 scores + motif counts)
alongside conservation/annotation scores.

FAMILY-WISE (leave-one-family-out) PREDICTION
---------------------------------------------
Genes and pseudogenes within a family (e.g. RNU1-*, RNU1-*P) are close
homologs, so a random train/test split leaks information: a pseudogene's
nearest relative is almost always in the training set. Instead, every
probability here is out-of-family:

    for each family F:
        train on all labeled, non-ambiguous genes NOT in F
        predict every gene in F (labeled and ambiguous alike)

So RNU1 probabilities come from a model that never saw any RNU1 gene or
pseudogene. This is analogous to out-of-fold prediction, except the folds
are biological families rather than random subsets.

Pipeline:
  1. Load the base conservation/annotation table (combined_gene_data.csv)
     and DEDUPLICATE it by gene (fixing VTRNA gene group)
  2. Load the SNP enrichment table and deduplicate by gene
  3. Load the promoter motif table and deduplicate by gene
  4. Merge all three on gene name -> results/combined_gene_data_with_snp_promoters.csv
     (asserting exactly one row per gene afterwards)
  5. Use the existing gene group column as the family
  6. Leave-one-family-out training/prediction for ALL genes
  7. Evaluate out-of-family performance (overall + per family)
  8. Report ambiguous gene predictions (also out-of-family)
  9. Plot probability distributions and enrichment/promoter-colored scatters
"""

import os
import re
import pandas as pd
import numpy as np
from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import classification_report, roc_auc_score, accuracy_score
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib.backends.backend_pdf import PdfPages

# ---------------------------------------------------------------------------
# Config
# ---------------------------------------------------------------------------

BASE_DATA_PATH = 'results/combined_gene_data.csv'
SNP_ENRICHMENT_PATH = 'data/snp_intergenic_genes_pseudogenes.csv'
PROMOTER_MOTIF_PATH = 'results/mcast_summary_grouped_annotated_with_snp.csv'
MERGED_DATA_PATH = 'results/combined_gene_data_with_snp_promoters.csv'
AMBIGUOUS_PREDICTIONS_PATH = 'results/ambiguous_gene_predictions.csv'
FAMILY_METRICS_PATH = 'results/rf_family_wise_metrics.csv'

FEATURE_IMPORTANCE_PATH = 'results/rf_feature_importance.csv'
FEATURE_IMPORTANCE_PER_FAMILY_PATH = 'results/rf_feature_importance_per_family.csv'

# Where dropped duplicate rows are written out for review
DUPLICATES_DIR = 'results/duplicates'

# Where to save plots. Each plot is saved as its own PNG (individually
# reviewable) and all plots are also collected into one multi-page PDF.
PLOTS_DIR = 'results/plots'
PLOTS_PDF_PATH = 'results/plots/all_plots.pdf'

# Gene group handling. Genes whose symbol starts with "VTRNA" (vault RNAs)
# must be in the VTRNA gene group, NOT the tRNA group. Adjust the label if
# your table spells the group differently (e.g. 'vtRNA').
GENE_GROUP_COL = 'Gene_group'
VTRNA_GROUP_LABEL = 'VTRNA'

# Families are taken directly from the existing gene group column
# (RNU1, RNU2, RNU4, RNU5, RNU6, RNU7, RNU11, RNU12, RNU4ATAC, RNU6ATAC,
# RN7SK, RN7SL, TRNA, VTRNA). Each group = one family, genes and
# pseudogenes together. No re-classification is done here.
FAMILY_COL = GENE_GROUP_COL

# Minimum number of training genes of EACH class required to fit a model for
# a held-out family. Families that would leave too little training data are
# reported and left as NaN rather than predicted from a degenerate model.
MIN_TRAIN_PER_CLASS = 5

RF_PARAMS = dict(n_estimators=100, random_state=42, n_jobs=-1)

# ---------------------------------------------------------------------------
# tRNA downsampling (training only)
# ---------------------------------------------------------------------------
# ~500/532 functional ncRNAs are tRNAs, so unrestricted they dominate every
# model trained on non-tRNA-held-out folds. We cap how many tRNAs enter the
# TRAINING pool; all tRNAs are still scored and evaluated.
#   'isotype'   -> cap per amino acid           (e.g. Ala, Gly, ...)
#   'anticodon' -> cap per amino acid+anticodon (e.g. Ala-AGC, Ala-CGC, ...)
TRNA_GROUP_LABEL = 'TRNA'        # value of Gene_group for tRNAs
TRNA_CAP_LEVEL = 'isotype'     # 'isotype' or 'anticodon'
TRNA_MAX_PER_GROUP = 3           # genes kept per group, per class
TRNA_CAP_PSEUDOGENES = True      # also cap tRNA pseudogenes (keeps classes balanced)
TRNA_SAMPLE_SEED = 42

# Conservation / annotation features from the original combined_gene_data.csv
CONSERVATION_FEATURES = [
    #'PhastCons30_median',
    'PhyloP100_median',
    #'PhyloP447_median',
    #'GTEX_max',
    'ENCODE_max',
]

# SNP enrichment features pulled in from snp_intergenic_genes_pseudogenes.csv
# (enrichment = flank-normalized SNP density relative to gene-body density)
ENRICHMENT_FEATURES = [
    #'enrichment_intergenic_1000genomes',
    'enrichment_intergenic_dbsnp',
    #'enrichment_intergenic_pangenome',
    #'enrichment_intergenic_gnomad'
]

# Promoter motif features pulled in from promoter_motifs.csv
# (type1/2/3 refer to distinct promoter motif classes; score = motif match
# strength, motif_counts = number of matches found)
PROMOTER_FEATURES = [
    'type1_score',
    'type1_motif_counts',
    'type2_score',
    'type2_motif_counts',
    'type3_score',
    'type3_motif_counts',
]

# Features to use for training the Random Forest:
# Conservation + SNP enrichment + promoter motif features
FEATURES = CONSERVATION_FEATURES + ENRICHMENT_FEATURES + PROMOTER_FEATURES

# List of ambiguous genes (assuming you have them defined)
ambiguous_genes = [
    "RNU1-128P", "RNU1-27P", "RNU1-28P", "RNU1-40P", "RNU1-41P", "RNU1-95P", "RNU1-97P", "RNU1-107P", "RNU1-48P", "RNU1-86P", "RNU2-57P", 
    "RNU2-63P", "RNU4-77P", "RNU5D-2P", "RNU5E-9P", "RNU5F-2P", "RNU5E-4P", "RNU5A-8P", "RNU5E-6P", "RNU6-1076P", "RNU6-109P", "RNU6-1171P", 
    "RNU6-184P", "RNU6-785P", "RNU6-1071P", "RNU6-465P", "RNU6-1239P", "RNU6-1334P", "RNU6-1268P", "RNU6-749P", "RNU6-18P", "RNU6-303P", "RNU6-255P", 
    "RNU6-1314P", "RNU6-941P", "RNU6-521P", "RNU6-631P", "RNU6-1235P", "RNU6-816P", "RNU6-368P", "RNU6-1269P", "RNU6-312P", "RNU6-473P", "RNU6-704P", 
    "RNU6-1193P", "RNU6-1321P", "RNU6-121P", "RNU6-1206P", "RNU6-546P", "RNU6-674P", "RNU6-89P", "RNU6-1283P", "RNU6-856P", "RNU6-162P", "RNU6-1096P", 
    "RNU6-436P", "RNU6-1056P", "RNU6-498P", "RNU6-501P", "RNU6-887P", "RNU6-424P", "RNU6-670P", "RNU6-747P", "RNU6-726P", "RNU6-440P", "RNU6-781P", 
    "RNY3P16", "RNY4P16", "TRQ-CTG16-1", "TRUND-NNN10-1", "TRV-CAC11-1", "TRV-CAC11-2", "TRA-TGC9-1", "TRS-ACT1-1", "RN7SL818P", "RN7SL725P", 
    "RN7SL763P", "RN7SL199P", "RN7SL544P", "RN7SL656P", "RN7SL640P", "RN7SL722P", "RN7SL702P", "RN7SL331P", "RN7SL717P", "RN7SL311P", "RN7SL134P", 
    "RNU7-182P", "RNU7-56P", "RNU7-110P", "RNU7-76P", "RNU7-78P", "RNU7-183P", "RN7SKP282", "RN7SKP131", "RN7SKP179", "RNU1-125P", "RN7SL471P", 
    "RN7SL4P", "RN7SL5P", "RN7SL128P", "RN7SKP255", "RNU2-4P", "RNU4-63P", "RNU5B-5P", "RNU6-1194P", "RNU6-1059P", "RNU6-202P", "RNU6-693P", 
    "RNU6-751P", "RNU4ATAC15P", "RNU4ATAC5P", "RNU4ATAC13P", "RNU4ATAC4P", "RNU4ATAC9P", "RNU4ATAC18P", "RNU4ATAC11P", "RNU4ATAC2P", "RNU4ATAC10P", 
    "RNU6ATAC3P", "RNU6ATAC34P", "RNU6ATAC26P", "RNU6ATAC11P", "RNU6ATAC4P", "RNU6ATAC27P", "RNU6ATAC2P", "RNU6ATAC8P", "RNU6ATAC38P", "RNU6ATAC5P", 
    "RNU6ATAC33P", "RNU6ATAC30P", "RNU6ATAC29P", "RNU6ATAC6P", "RNU6ATAC7P", "RNU6ATAC37P", "RNU6ATAC20P", "RNU6ATAC41P", "RNU6ATAC10P", "RNU6ATAC21P", 
    "RNU6ATAC15P", "RNU6ATAC31P", "RNU6ATAC18P", "RNU6ATAC25P", "RNU6ATAC16P", "RNU6ATAC23P", "RNU11-3P", "RNU11-6P", "RNU11-4P", "RNU12-2P", 
    "TRL-TAA5-1", "TRQ-TTG10-1", "RN7SKP70", "RN7SKP253", "RN7SKP12", "RN7SKP123"
]

# ---------------------------------------------------------------------------
# Deduplication helpers
# ---------------------------------------------------------------------------

os.makedirs(DUPLICATES_DIR, exist_ok=True)


def is_vtrna(gene_series):
    """True for vault RNA genes (VTRNA1-1, VTRNA2-1, ...)."""
    return gene_series.astype(str).str.upper().str.startswith('VTRNA')


def looks_like_trna_group(group_series):
    """True where the gene group label says tRNA (but not VTRNA)."""
    g = group_series.astype(str).str.lower()
    return g.str.contains('trna') & ~g.str.contains('vtrna')


def dedupe_by_gene(df, key, label, value_cols=None):
    """Return df with exactly one row per `key`.

    1. Strip whitespace from the key.
    2. Drop fully identical rows.
    3. For rows that still share a key (conflicting duplicates), keep the
       best one:
         - for VTRNA genes, prefer a row whose gene group is NOT tRNA
         - otherwise prefer the row with the most non-missing values
           (ties -> first occurrence in the file)
    Dropped conflicting rows are written to results/duplicates/<label>.csv
    so they can be reviewed.
    """
    df = df.copy()
    df[key] = df[key].astype(str).str.strip()
    n_start = len(df)

    df = df.drop_duplicates()
    n_exact = n_start - len(df)

    dup_mask = df.duplicated(key, keep=False)
    n_conflict_rows = 0
    if dup_mask.any():
        df.loc[dup_mask].sort_values(key).to_csv(
            os.path.join(DUPLICATES_DIR, f'{label}_conflicting_duplicates.csv'),
            index=False
        )

        cols = value_cols or [c for c in df.columns if c != key]
        df['_n_nonnull'] = df[cols].notna().sum(axis=1)

        # VTRNA preference: rows with a correct (non-tRNA) group win
        df['_group_ok'] = 1
        if GENE_GROUP_COL in df.columns:
            bad = is_vtrna(df[key]) & looks_like_trna_group(df[GENE_GROUP_COL])
            df.loc[bad, '_group_ok'] = 0

        before = len(df)
        df = (df.sort_values(['_group_ok', '_n_nonnull'], ascending=False,
                             kind='stable')
                .drop_duplicates(key, keep='first')
                .drop(columns=['_n_nonnull', '_group_ok'])
                .sort_index())
        n_conflict_rows = before - len(df)

    print(f"[dedupe:{label}] {n_start} rows -> {len(df)} rows "
          f"({n_exact} exact duplicates, {n_conflict_rows} conflicting "
          f"duplicates removed)")
    return df


def fix_vtrna_gene_group(df, gene_col='Gene'):
    """Force VTRNA genes into the VTRNA gene group (never tRNA)."""
    if GENE_GROUP_COL not in df.columns:
        print(f"Warning: no '{GENE_GROUP_COL}' column found; skipping VTRNA fix.")
        return df
    mask = is_vtrna(df[gene_col])
    n_wrong = (mask & (df[GENE_GROUP_COL] != VTRNA_GROUP_LABEL)).sum()
    df.loc[mask, GENE_GROUP_COL] = VTRNA_GROUP_LABEL
    print(f"[gene_group] {mask.sum()} VTRNA genes set to '{VTRNA_GROUP_LABEL}' "
          f"({n_wrong} were previously something else)")
    return df

# ---------------------------------------------------------------------------
# tRNA parsing / downsampling helpers
# ---------------------------------------------------------------------------

# Nuclear:           TRA-AGC10-1      -> isotype A,   anticodon AGC
# Nuclear mito-like: NMTRQ-TTG10-1    -> isotype Q,   anticodon TTG
# Suppressor:        TRSUP-CTA2-1     -> isotype SUP, anticodon CTA
# Undetermined:      TRUND-NNN10-1    -> isotype UND, anticodon NNN
# Mitochondrial:     MT-TS1, MT-TQ    -> isotype S (anticodon proxy 'S1'), Q
_NUC_RE = re.compile(r'^(NMTR|TR)([A-Z]|SUP|UND)-([ACGTUN]{3})\d*-\d+P?$')
_MT_RE = re.compile(r'^MT-T([A-Z])(\d?)$')


def parse_trna(symbol):
    """Return (source, isotype, anticodon) for an HGNC tRNA symbol, or
    ('UNPARSED', 'UNPARSED', symbol) if it doesn't match a known pattern."""
    s = str(symbol).strip().upper()
    m = _NUC_RE.match(s)
    if m:
        source = 'NMT' if m.group(1) == 'NMTR' else 'NUC'
        return source, m.group(2), m.group(3)
    m = _MT_RE.match(s)
    if m:
        # MT genes have no anticodon in the name; the optional digit
        # (MT-TS1 vs MT-TS2, MT-TL1 vs MT-TL2) distinguishes isoacceptors.
        return 'MT', m.group(1), m.group(1) + m.group(2)
    return 'UNPARSED', 'UNPARSED', s


def cap_trna_training_pool(df, train_mask):
    """Return a copy of train_mask where tRNAs are downsampled to at most
    TRNA_MAX_PER_GROUP genes per (class, source, isotype[, anticodon]).

    Source (NUC / NMT / MT) is part of the key so the mitochondrial and
    nuclear-mitochondrial-like tRNAs stay represented. Non-tRNA genes are
    untouched. Selection is random but seeded (reproducible)."""
    is_trna = df['family'].str.upper() == TRNA_GROUP_LABEL.upper()
    cand = df[train_mask & is_trna].copy()
    if not TRNA_CAP_PSEUDOGENES:
        cand = cand[cand['label'] == 1]

    parsed = cand['Gene'].apply(parse_trna)
    cand['_source'] = parsed.str[0]
    cand['_isotype'] = parsed.str[1]
    cand['_anticodon'] = parsed.str[2]

    n_unparsed = (cand['_source'] == 'UNPARSED').sum()
    if n_unparsed:
        print(f"Warning: {n_unparsed} tRNA symbols didn't match a known "
              f"pattern and are pooled into one 'UNPARSED' group: "
              f"{cand.loc[cand['_source'] == 'UNPARSED', 'Gene'].tolist()[:10]}")

    group_cols = ['label', '_source', '_isotype']
    if TRNA_CAP_LEVEL == 'anticodon':
        group_cols.append('_anticodon')
    elif TRNA_CAP_LEVEL != 'isotype':
        raise ValueError("TRNA_CAP_LEVEL must be 'isotype' or 'anticodon'")

    kept = (cand.groupby(group_cols, group_keys=False)
                .apply(lambda g: g.sample(n=min(len(g), TRNA_MAX_PER_GROUP),
                                          random_state=TRNA_SAMPLE_SEED)))

    dropped_idx = cand.index.difference(kept.index)
    capped_mask = train_mask.copy()
    capped_mask.loc[dropped_idx] = False

    print(f"\n[tRNA cap] level={TRNA_CAP_LEVEL}, max {TRNA_MAX_PER_GROUP}/group")
    print(f"  tRNAs in training pool: {len(cand)} -> {len(kept)} "
          f"({int((kept['label'] == 1).sum())} functional, "
          f"{int((kept['label'] == 0).sum())} pseudogenes)")
    print(f"  groups represented: "
          f"{kept[group_cols[1:]].drop_duplicates().shape[0]}")
    return capped_mask

# ---------------------------------------------------------------------------
# Step 0a: Load base data (dedupe) and merge SNP enrichment data (dedupe)
# ---------------------------------------------------------------------------

print("Loading base conservation/annotation data...")
base_data = pd.read_csv(BASE_DATA_PATH)
base_data = dedupe_by_gene(base_data, key='Gene', label='base_combined_gene_data')
base_data = fix_vtrna_gene_group(base_data, gene_col='Gene')

print("Loading SNP enrichment data...")
snp_data = pd.read_csv(SNP_ENRICHMENT_PATH)

# Only pull in the enrichment-related columns from the SNP table so we don't
# clobber columns that already exist in the base table (e.g. gene_group).
snp_cols_to_merge = ['gene_name'] + [
    c for c in snp_data.columns
    if c.startswith('snp_density_')
    or c.startswith('enrichment_intergenic_')
]
snp_data_subset = dedupe_by_gene(
    snp_data[snp_cols_to_merge], key='gene_name', label='snp_enrichment'
)

# Left join so every gene in the base table is kept, even if it has no
# matching entry in the SNP enrichment table (those rows get NaN enrichment).
# validate='m:1' raises if the right-hand table somehow still has duplicates.
data = base_data.merge(
    snp_data_subset,
    left_on='Gene',
    right_on='gene_name',
    how='left',
    validate='m:1'
).drop(columns=['gene_name'])

# ---------------------------------------------------------------------------
# Step 0b: Merge in promoter motif data (type1/2/3 scores + motif counts)
# ---------------------------------------------------------------------------

print("Loading promoter motif data...")
promoter_data = pd.read_csv(PROMOTER_MOTIF_PATH)

# Only pull in the motif score/count columns so we don't duplicate columns
# already present (e.g. enrichment_intergenic_dbsnp is in both this table
# and the SNP enrichment table -- we keep the version already merged above).
promoter_cols_to_merge = ['gene_name'] + PROMOTER_FEATURES
promoter_data_subset = promoter_data[promoter_cols_to_merge].copy()

# The promoter motif columns come in with blank strings ("") where no motif
# of that type was found, so coerce to numeric BEFORE deduplicating (so that
# "most non-missing values" is judged on real numbers).
for col in PROMOTER_FEATURES:
    promoter_data_subset[col] = pd.to_numeric(promoter_data_subset[col], errors='coerce')

promoter_data_subset = dedupe_by_gene(
    promoter_data_subset, key='gene_name', label='promoter_motifs',
    value_cols=PROMOTER_FEATURES
)

data = data.merge(
    promoter_data_subset,
    left_on='Gene',
    right_on='gene_name',
    how='left',
    validate='m:1'
).drop(columns=['gene_name'])

n_missing_promoter = data[PROMOTER_FEATURES].isna().any(axis=1).sum()
if n_missing_promoter:
    print(f"Note: {n_missing_promoter} genes have missing promoter motif "
          f"values (no motif of that type detected, or no match in the "
          f"promoter table). Filling missing motif scores/counts with 0.")
    data[PROMOTER_FEATURES] = data[PROMOTER_FEATURES].fillna(0)

# Final safety net: one row per gene, and VTRNA group still correct
assert data['Gene'].is_unique, (
    f"{data['Gene'].duplicated().sum()} duplicated genes remain after merging"
)
data = fix_vtrna_gene_group(data, gene_col='Gene')

# ---------------------------------------------------------------------------
# Step 1: Labels, ambiguous flag, and gene families
# ---------------------------------------------------------------------------

data['label'] = data['Gene_Type'].apply(lambda x: 1 if x == 'Functional' else 0)
data['is_ambiguous'] = data['Gene'].isin(ambiguous_genes)
if data[FAMILY_COL].isna().any():
    n_na = data[FAMILY_COL].isna().sum()
    print(f"Warning: {n_na} genes have no '{FAMILY_COL}' value; they are "
          f"treated as their own family 'UNASSIGNED'.")
data['family'] = data[FAMILY_COL].fillna('UNASSIGNED').astype(str).str.strip()

print("\nTotal counts in full dataset:")
print(f"Functional Genes: {data[data['label'] == 1].shape[0]}")
print(f"Pseudogenes: {data[data['label'] == 0].shape[0]}")

family_sizes = (data.groupby('family')
                    .agg(n_genes=('Gene', 'size'),
                         n_functional=('label', 'sum'),
                         n_ambiguous=('is_ambiguous', 'sum'))
                    .sort_values('n_genes', ascending=False))
print(f"\n{len(family_sizes)} gene families detected. Largest families:")
print(family_sizes.head(15))

# Genes eligible to be used as TRAINING data: labeled and not ambiguous.
# Labeled, non-ambiguous genes: used for EVALUATION (all of them).
labeled_mask = ~data['is_ambiguous']

# Genes eligible for TRAINING: same, but with tRNAs downsampled for diversity.
train_pool_mask = cap_trna_training_pool(data, labeled_mask)

tp = data[train_pool_mask]
print(f"\nTraining pool after tRNA cap: {len(tp)} genes "
      f"({int(tp['label'].sum())} functional, "
      f"{int((tp['label'] == 0).sum())} pseudogenes)")
print(f"  tRNAs: {int((tp['family'].str.upper() == TRNA_GROUP_LABEL.upper()).sum())}")

# ---------------------------------------------------------------------------
# Step 2: Leave-one-family-out training + prediction for ALL genes
# ---------------------------------------------------------------------------
# For each family F, a fresh Random Forest is trained on every training-pool
# gene outside F, and used to predict every gene inside F. Each gene's
# probability therefore comes from a model that has never seen any member
# (gene or pseudogene) of its own family.

print("\nRunning leave-one-family-out training/prediction...")
data['functional_probability'] = np.nan
importances_per_family = []
fitted_families = []
skipped_families = []

for fam in family_sizes.index:
    in_fam = data['family'] == fam
    train_mask = train_pool_mask & ~in_fam

    y_tr = data.loc[train_mask, 'label']
    n_pos, n_neg = int(y_tr.sum()), int(len(y_tr) - y_tr.sum())
    if n_pos < MIN_TRAIN_PER_CLASS or n_neg < MIN_TRAIN_PER_CLASS:
        skipped_families.append((fam, n_pos, n_neg))
        continue

    model_f = RandomForestClassifier(**RF_PARAMS)
    model_f.fit(data.loc[train_mask, FEATURES], y_tr)
    data.loc[in_fam, 'functional_probability'] = model_f.predict_proba(
        data.loc[in_fam, FEATURES])[:, 1]
    importances_per_family.append(model_f.feature_importances_)
    fitted_families.append(fam)

if skipped_families:
    print("\nWarning: these families were NOT predicted (too little training "
          "data of one class once the family was held out):")
    for fam, n_pos, n_neg in skipped_families:
        print(f"  {fam}: training would have {n_pos} functional / {n_neg} pseudogenes")

n_scored = data['functional_probability'].notna().sum()
print(f"Scored {n_scored} / {len(data)} genes out-of-family.")

# Mean feature importance across all family-held-out models
# Feature importance summary across all leave-one-family-out models.
# Each fitted model contributes one importance vector; SD is the spread across
# those models (ddof=1) and SE = SD / sqrt(n_models).
imp_matrix = np.array(importances_per_family)          # (n_models, n_features)
n_models = imp_matrix.shape[0]
imp_sd = imp_matrix.std(axis=0, ddof=1) if n_models > 1 else np.zeros(len(FEATURES))
imp_se = imp_sd / np.sqrt(n_models)

importance_df = (pd.DataFrame({
                    'feature': FEATURES,
                    'importance_mean': imp_matrix.mean(axis=0),
                    'importance_sd': imp_sd,
                    'importance_se': imp_se,
                    'n_models': n_models,
                 })
                 .sort_values('importance_mean', ascending=False)
                 .reset_index(drop=True))

os.makedirs(os.path.dirname(FEATURE_IMPORTANCE_PATH), exist_ok=True)
importance_df.to_csv(FEATURE_IMPORTANCE_PATH, index=False)

# Raw per-model importances (one row per held-out family) for full transparency
(pd.DataFrame(imp_matrix, columns=FEATURES)
   .assign(held_out_family=fitted_families)
   [['held_out_family'] + FEATURES]
   .to_csv(FEATURE_IMPORTANCE_PER_FAMILY_PATH, index=False))

print(f"\nFeature Importances (mean ± SD / SE across {n_models} "
      f"leave-one-family-out models; saved to {FEATURE_IMPORTANCE_PATH}):")
print(importance_df.to_string(index=False))

# ---------------------------------------------------------------------------
# Step 3: Evaluate out-of-family performance on labeled, non-ambiguous genes
# ---------------------------------------------------------------------------

eval_df = data[labeled_mask & data['functional_probability'].notna()].copy()
eval_df['pred'] = (eval_df['functional_probability'] >= 0.5).astype(int)

print("\nOut-of-family evaluation "
      f"({len(eval_df)} labeled non-ambiguous genes):")
print("\nClassification Report:")
print(classification_report(eval_df['label'], eval_df['pred']))
print(f"Out-of-family ROC AUC Score: "
      f"{roc_auc_score(eval_df['label'], eval_df['functional_probability']):.2f}")

# Per-family breakdown (AUC only defined where a family has both classes)
rows = []
for fam, g in eval_df.groupby('family'):
    both = g['label'].nunique() == 2
    rows.append({
        'family': fam,
        'n_genes': len(g),
        'n_functional': int(g['label'].sum()),
        'n_pseudogenes': int(len(g) - g['label'].sum()),
        'accuracy': accuracy_score(g['label'], g['pred']),
        'roc_auc': roc_auc_score(g['label'], g['functional_probability']) if both else np.nan,
        'mean_functional_probability': g['functional_probability'].mean(),
    })
family_metrics = pd.DataFrame(rows).sort_values('n_genes', ascending=False)
os.makedirs(os.path.dirname(FAMILY_METRICS_PATH), exist_ok=True)
family_metrics.to_csv(FAMILY_METRICS_PATH, index=False)
print(f"\nPer-family out-of-family metrics (saved to {FAMILY_METRICS_PATH}):")
print(family_metrics.head(15).to_string(index=False))

# ---------------------------------------------------------------------------
# Step 4: Ambiguous genes (also predicted out-of-family)
# ---------------------------------------------------------------------------

ambiguous_data = data[data['is_ambiguous']].copy()

if ambiguous_data.empty:
    print("\nNo ambiguous genes found in the dataset.")
else:
    print("\nAmbiguous Genes Counts:")
    print(f"Functional Genes: {ambiguous_data[ambiguous_data['label'] == 1].shape[0]}")
    print(f"Pseudogenes: {ambiguous_data[ambiguous_data['label'] == 0].shape[0]}")

    ambiguous_output = pd.DataFrame({
        'HGNC Gene Symbol': ambiguous_data['Gene'],
        'Family': ambiguous_data['family'],
        'Current Classification': ambiguous_data['Gene_Type'],
        'Predicted Classification': np.where(
            ambiguous_data['functional_probability'].isna(), 'Not predicted',
            np.where(ambiguous_data['functional_probability'] >= 0.5,
                     'Functional', 'Pseudogene')
        ),
        'Random Forest Predicted Functional Probability': ambiguous_data['functional_probability'],
    })

    ambiguous_output.to_csv(AMBIGUOUS_PREDICTIONS_PATH, index=False)

    print("\nPredictions for ambiguous genes (each from a model that excluded "
          "the gene's own family):")
    print(ambiguous_output)

# ---------------------------------------------------------------------------
# Step 5: Save the full dataset with out-of-family probabilities
# ---------------------------------------------------------------------------

data.to_csv(MERGED_DATA_PATH, index=False)
print(f"\nSaved {MERGED_DATA_PATH} ({len(data)} unique genes) with family "
      f"and out-of-family functional_probability columns")

# For the train-vs-out-of-family comparison plot: a single model fit on the
# whole training pool, scored in-sample, shows how much optimism a
# non-family-aware evaluation would have had.
full_model = RandomForestClassifier(**RF_PARAMS)
full_model.fit(data.loc[train_pool_mask, FEATURES], data.loc[train_pool_mask, 'label'])
insample_probs = full_model.predict_proba(data.loc[eval_df.index, FEATURES])[:, 1]
