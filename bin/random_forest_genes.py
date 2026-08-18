"""
Random Forest classifier for functional gene vs. pseudogene prediction,
now incorporating SNP enrichment features (1000 Genomes, gnomAD,
pangenome, dbSNP) alongside conservation/annotation scores.

Pipeline:
  1. Load the base conservation/annotation table (combined_gene_data.csv)
  2. Load the SNP enrichment table (snp_enrichment_genes_pseudogenes.csv)
  3. Merge them on gene name -> results/combined_gene_data_with_snp.csv
  4. Train/evaluate a Random Forest using conservation + enrichment features
  5. Score ambiguous genes and the full dataset
  6. Plot probability distributions and an enrichment-colored scatter
"""

import pandas as pd
import numpy as np
from sklearn.ensemble import RandomForestClassifier
from sklearn.model_selection import train_test_split
from sklearn.metrics import classification_report, roc_auc_score
import matplotlib.pyplot as plt
import seaborn as sns

# ---------------------------------------------------------------------------
# Config
# ---------------------------------------------------------------------------

BASE_DATA_PATH = 'results/combined_gene_data.csv'
SNP_ENRICHMENT_PATH = 'data/snp_enrichment_genes_pseudogenes.csv'
MERGED_DATA_PATH = 'results/combined_gene_data_with_snp.csv'
AMBIGUOUS_PREDICTIONS_PATH = 'results/ambiguous_gene_predictions.csv'

# Conservation / annotation features from the original combined_gene_data.csv
CONSERVATION_FEATURES = [
    'PhastCons30_median',
    'PhyloP100_median',
    'PhyloP447_median',
    'GTEX_max',
    'ENCODE_max',
]

# SNP enrichment features pulled in from snp_enrichment_genes_pseudogenes.csv
# (enrichment = flank-normalized SNP density relative to gene-body density)
ENRICHMENT_FEATURES = [
    'enrichment_1000genomes',
    'enrichment_gnomad',
    'enrichment_pangenome',
    'enrichment_dbsnp',
]

# Features to use for training the Random Forest:
# Conservation + enrichment + expression (ENCODE) features
FEATURES = [
    'PhyloP100_median',
    'ENCODE_max',
    'enrichment_dbsnp',
]

# List of ambiguous genes (assuming you have them defined)
ambiguous_genes = ['RNU6-1189P', 'RNU6-82P', 'RNU2-2P', 'RNU1-27P', 'RNU1-28P', 'RNU5B-1',
                   'RNU5F-1', 'RNU5B-5P', 'RNU6-1194P', 'RNU6-1334P', 'RN7SKP70', 'RN7SL471P',
                    'TRL-TAA5-1', 'RNU6ATAC10P','MT-TF', 'MT-TL1', 'NMTRS-TGA3-1',
                    'TRG-CCC7-1', 'TRG-TCC2-3', 'TRG-TCC2-4', 'TRD-GTC2-4', 'TRD-GTC2-3', 'TRP-AGG5-1',
                    'TRE-TTC15-1', 'TRK-TTT8-1', 'TRA-AGC15-1', 'TRK-TTT14-1', 'TRP-GGG1-1', 'TRV-CAC1-3', 'TRD-GTC2-5', 'TRSUP-CTA2-1', 'TRE-CTC1-3', 'TRG-TCC2-2',
                    'TRQ-TTG9-1', 'TRL-CAA7-1', 'TRE-CTC1-4', 'TRV-CAC1-2', 'TRSUP-CTA3-1', 'TRD-GTC2-2', 'TRG-TCC2-5', 'TRUND-NNN5-1', 'TRQ-TTG8-1', 'TRR-CCT6-2', 'TRE-CTC14-1',
                    'TRL-CAG1-5', 'TRA-TGC10-1', 'TRL-CAG1-1', 'TRL-CAG1-4', 'TRQ-CTG1-3', 'TRL-AAG1-2', 'TRV-CAC13-1', 'TRE-CTC1-5', 'TRE-CTC1-2', 'TRS-AGA7-1', 'TRE-TTC14-1', 'TRC-GCA7-1', 'TRL-AAG1-1',
                    'TRE-CTC13-1', 'TRL-CAG1-2', 'TRD-GTC2-1', 'TRUND-NNN1-1', 'TRK-CTT9-1', 'TRE-CTC15-1', 'TRE-CTC12-1', 'TRE-CTC4-1', 'TRU-TCA2-1', 'TRD-GTC10-1',
                    'TRK-CTT6-1', 'TRE-CTC11-1', 'TRE-CTC9-1', 'TRE-TTC4-2', 'TRG-CCC4-1', 'TRN-GTT2-7', 'TRV-AAC1-3', 'TRN-GTT2-8', 'TRV-CAC8-1', 'TRG-GCC1-4',
                    'TRG-TCC4-1', 'TRL-CAG1-3', 'TRK-TTT2-1', 'TRA-AGC9-1', 'TRL-CAA5-1', 'TRK-CTT7-1', 'TRY-ATA1-1', 'TRK-CTT5-1', 'TRR-CCT8-1']


# ---------------------------------------------------------------------------
# Step 0: Merge base conservation data with SNP enrichment data
# ---------------------------------------------------------------------------

print("Loading base conservation/annotation data...")
base_data = pd.read_csv(BASE_DATA_PATH)

print("Loading SNP enrichment data...")
snp_data = pd.read_csv(SNP_ENRICHMENT_PATH)

# Only pull in the enrichment-related columns from the SNP table so we don't
# clobber columns that already exist in the base table (e.g. gene_group).
snp_cols_to_merge = ['gene_name'] + [
    c for c in snp_data.columns
    if c.startswith('snp_count_')
    or c.startswith('snp_density_')
    or c.startswith('flank_count_')
    or c.startswith('flank_density_')
    or c.startswith('enrichment_')
]
snp_data_subset = snp_data[snp_cols_to_merge]

# Left join so every gene in the base table is kept, even if it has no
# matching entry in the SNP enrichment table (those rows get NaN enrichment).
data = base_data.merge(
    snp_data_subset,
    left_on='Gene',
    right_on='gene_name',
    how='left'
).drop(columns=['gene_name'])

n_unmatched = data['enrichment_1000genomes'].isna().sum()
if n_unmatched:
    print(f"Warning: {n_unmatched} genes had no matching SNP enrichment data. "
          f"Filling missing enrichment values with 0.")
    for col in snp_cols_to_merge[1:]:
        data[col] = data[col].fillna(0)

data.to_csv(MERGED_DATA_PATH, index=False)
print(f"Saved merged dataset with SNP enrichment features to {MERGED_DATA_PATH}")

# ---------------------------------------------------------------------------
# Step 1: Labels + exclude ambiguous genes
# ---------------------------------------------------------------------------

data['label'] = data['Gene_Type'].apply(lambda x: 1 if x == 'Functional' else 0)

data_non_ambiguous = data[~data['Gene'].isin(ambiguous_genes)]

print("\nTotal counts in full dataset:")
print(f"Functional Genes: {data[data['label'] == 1].shape[0]}")
print(f"Pseudogenes: {data[data['label'] == 0].shape[0]}")

# ---------------------------------------------------------------------------
# Step 2: Train/test split (90/10) using conservation + enrichment features
# ---------------------------------------------------------------------------

X = data_non_ambiguous[FEATURES]
y = data_non_ambiguous['label']

X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.1, random_state=42)

print("\nTraining Set Counts:")
print(f"Functional Genes: {y_train.sum()}")
print(f"Pseudogenes: {len(y_train) - y_train.sum()}")

print("\nTest Set Counts:")
print(f"Functional Genes: {y_test.sum()}")
print(f"Pseudogenes: {len(y_test) - y_test.sum()}")

# ---------------------------------------------------------------------------
# Step 3: Train the Random Forest
# ---------------------------------------------------------------------------

model = RandomForestClassifier(n_estimators=100, random_state=42)
print(f"\nTraining the Random Forest model on features: {FEATURES}")
model.fit(X_train, y_train)

# Feature importances so you can see how much the SNP enrichment features
# are contributing relative to the conservation scores.
importances = pd.Series(model.feature_importances_, index=FEATURES).sort_values(ascending=False)
print("\nFeature Importances:")
print(importances)

# ---------------------------------------------------------------------------
# Step 4: Evaluate on the test set
# ---------------------------------------------------------------------------

print("\nEvaluating the model on the test set...")
test_predictions = model.predict(X_test)
test_probs = model.predict_proba(X_test)[:, 1]

print("\nTest Classification Report:")
print(classification_report(y_test, test_predictions))
print(f"Test ROC AUC Score: {roc_auc_score(y_test, test_probs):.2f}")

# ---------------------------------------------------------------------------
# Step 5: Predict on ambiguous genes
# ---------------------------------------------------------------------------

ambiguous_data = data[data['Gene'].isin(ambiguous_genes)].copy()

if ambiguous_data.empty:
    print("\nNo ambiguous genes found in the dataset.")
else:
    print("\nAmbiguous Genes Counts:")
    print(f"Functional Genes: {ambiguous_data[ambiguous_data['label'] == 1].shape[0]}")
    print(f"Pseudogenes: {ambiguous_data[ambiguous_data['label'] == 0].shape[0]}")

    print("\nPredicting probabilities for ambiguous genes...")
    X_ambiguous = ambiguous_data[FEATURES]
    ambiguous_probs = model.predict_proba(X_ambiguous)[:, 1]
    ambiguous_data['functional_probability'] = ambiguous_probs

    ambiguous_data[['Gene', 'functional_probability']].to_csv(
        AMBIGUOUS_PREDICTIONS_PATH, index=False
    )

    print("\nPredictions for ambiguous genes:")
    print(ambiguous_data[['Gene', 'functional_probability']])

# ---------------------------------------------------------------------------
# Step 6: Score the full dataset
# ---------------------------------------------------------------------------

print("\nAdding predicted probabilities for all genes in the full dataset...")
data['functional_probability'] = model.predict_proba(data[FEATURES])[:, 1]

# ---------------------------------------------------------------------------
# Visualizations
# ---------------------------------------------------------------------------

test_functional_probs = test_probs[y_test == 1]
test_pseudogene_probs = test_probs[y_test == 0]

plt.rcParams['font.family'] = 'DejaVu Serif'
plt.rcParams.update({'font.size': 28,
                     'axes.labelsize': 28,
                     'xtick.labelsize': 28,
                     'ytick.labelsize': 28,
                     'legend.fontsize': 22,
                     'figure.titlesize': 28})

# --- Test set probability distribution ---
plt.figure(figsize=(10, 6))
sns.histplot(test_functional_probs, bins=20, color='firebrick', alpha=0.7, label='Functional Genes')
sns.histplot(test_pseudogene_probs, bins=20, color='cornflowerblue', alpha=0.7, label='Pseudogenes')
plt.xlabel('Probability of Being Functional')
plt.ylabel('Frequency')
plt.legend(title='Gene Type')
plt.grid(axis='y', linestyle='--', alpha=0.7)
plt.tight_layout()
plt.show()

# --- Train vs test overlay ---
train_probs = model.predict_proba(X_train)[:, 1]
train_functional_probs = train_probs[y_train == 1]
train_pseudogene_probs = train_probs[y_train == 0]

plt.figure(figsize=(10, 6))
sns.histplot(test_functional_probs, bins=20, color='firebrick', alpha=0.7, label='Test Functional Genes')
sns.histplot(test_pseudogene_probs, bins=20, color='cornflowerblue', alpha=0.7, label='Test Pseudogenes')
sns.histplot(train_functional_probs, bins=20, color='lightpink', alpha=0.5, label='Train Functional Genes', linestyle='--')
sns.histplot(train_pseudogene_probs, bins=20, color='lightblue', alpha=0.5, label='Train Pseudogenes', linestyle='--')
plt.xlabel('Probability of Being Functional')
plt.ylabel('Frequency')
plt.legend(title='Gene Type')
plt.grid(axis='y', linestyle='--', alpha=0.7)
plt.tight_layout()
plt.show()

# --- Conservation scatter, colored by functional probability (unchanged) ---
plt.figure(figsize=(10, 6))
ambiguous_data_plot = data[data['Gene'].isin(ambiguous_genes)]
plt.scatter(ambiguous_data_plot['PhyloP100_median'], ambiguous_data_plot['ENCODE_max'],
            c=ambiguous_data_plot['functional_probability'], cmap='coolwarm', s=100, edgecolors='black')
plt.xlabel('PhyloP100 Median')
plt.ylabel('ENCODE Max')
plt.colorbar(label='Functional Probability')
plt.grid(False)
plt.tight_layout()
plt.show()

# --- New: SNP enrichment scatter, colored by functional probability ---
# Uses gnomAD and dbSNP enrichment (largest, most robust variant catalogs)
# as the two axes to visualize how enrichment separates predicted classes.
plt.figure(figsize=(10, 6))
plt.scatter(
    ambiguous_data_plot['enrichment_gnomad'],
    ambiguous_data_plot['enrichment_dbsnp'],
    c=ambiguous_data_plot['functional_probability'],
    cmap='coolwarm', s=100, edgecolors='black'
)
plt.xlabel('gnomAD SNP Enrichment')
plt.ylabel('dbSNP SNP Enrichment')
plt.colorbar(label='Functional Probability')
plt.grid(False)
plt.tight_layout()
plt.show()