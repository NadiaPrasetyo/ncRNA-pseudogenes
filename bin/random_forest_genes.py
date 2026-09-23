"""
Random Forest classifier for functional gene vs. pseudogene prediction,
incorporating SNP enrichment features (1000 Genomes, gnomAD, pangenome,
dbSNP) AND promoter motif features (type 1/2/3 scores + motif counts)
alongside conservation/annotation scores.

Pipeline:
  1. Load the base conservation/annotation table (combined_gene_data.csv)
  2. Load the SNP enrichment table (snp_intergenic_genes_pseudogenes.csv)
  3. Load the promoter motif table (promoter_motifs.csv)
  4. Merge all three on gene name -> results/combined_gene_data_with_snp_promoters.csv
  5. Train/evaluate a Random Forest using conservation + enrichment + promoter features
  6. Score ambiguous genes and the full dataset
  7. Plot probability distributions and enrichment/promoter-colored scatters
"""

import pandas as pd
import numpy as np
from sklearn.ensemble import RandomForestClassifier
from sklearn.model_selection import train_test_split
from sklearn.metrics import classification_report, roc_auc_score
import os
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

# Where to save plots. Each plot is saved as its own PNG (individually
# reviewable) and all plots are also collected into one multi-page PDF.
PLOTS_DIR = 'results/plots'
PLOTS_PDF_PATH = 'results/plots/all_plots.pdf'

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
ambiguous_genes = ["RNU1-28P", "RNU1-27P", "RNU1-52P", "RNU2-2", "RNU2-4P", "RNU4-63P", "RNU5F-1", "RNU5B-5P",
 "RNU5E-6P", "RNU5A-8P", "RNU6-1334P", "RNU6-1194P", "RNU6-1059P", "RNU6-202P", "RNU6-693P", "RNU6-751P", "RNU6-42P",
 "RNU6-1078P", "RNU6-12P", "RNU6-231P", "RNU4ATAC15P", "RNU4ATAC5P", "RNU4ATAC13P", "RNU4ATAC4P", "RNU4ATAC9P", "RNU4ATAC18P",
 "RNU4ATAC11P", "RNU4ATAC2P", "RNU4ATAC10P", "RNU4ATAC16P", "RNU4ATAC12P", "RNU6ATAC3P", "RNU6ATAC34P", "RNU6ATAC26P",
 "RNU6ATAC11P", "RNU6ATAC4P", "RNU6ATAC27P", "RNU6ATAC2P", "RNU6ATAC8P", "RNU6ATAC38P", "RNU6ATAC5P", "RNU6ATAC33P",
 "RNU6ATAC30P", "RNU6ATAC29P", "RNU6ATAC6P", "RNU6ATAC7P", "RNU6ATAC37P", "RNU6ATAC20P", "RNU6ATAC41P", "RNU6ATAC10P",
 "RNU6ATAC21P", "RNU6ATAC15P", "RNU6ATAC31P", "RNU6ATAC18P", "RNU6ATAC25P", "RNU6ATAC16P", "RNU6ATAC23P", "RNU11-3P",
 "RNU11-6P", "RNU11-4P", "RNU12-2P", "VTRNA1-1", "VTRNA2-1", "TRR-CCT8-1", "TRK-CTT5-1", "TRK-CTT7-1", "TRY-ATA1-1",
 "TRL-CAA5-1", "TRA-AGC9-1", "TRK-TTT2-1", "TRL-CAG1-3", "TRG-TCC4-1", "TRG-GCC1-4", "TRV-CAC8-1", "TRN-GTT2-8",
 "TRV-AAC1-3", "TRN-GTT2-7", "TRG-CCC4-1", "TRE-TTC4-2", "TRE-CTC9-1", "TRE-CTC11-1", "TRK-CTT6-1", "TRD-GTC10-1",
 "TRU-TCA2-1", "TRE-CTC12-1", "TRE-CTC4-1", "TRE-CTC15-1", "TRK-CTT9-1", "TRUND-NNN1-1", "TRD-GTC2-1", "TRL-CAG1-2",
 "TRE-CTC13-1", "TRL-AAG1-1", "TRC-GCA7-1", "TRE-TTC14-1", "TRE-CTC1-2", "TRS-AGA7-1", "TRE-CTC1-5", "TRV-CAC13-1",
 "TRL-AAG1-2", "TRQ-CTG1-3", "TRL-CAG1-4", "TRL-CAG1-1", "TRA-TGC10-1", "TRL-CAG1-5", "TRE-CTC14-1", "TRR-CCT6-2",
 "TRQ-TTG8-1", "TRUND-NNN5-1", "TRG-TCC2-5", "TRD-GTC2-2", "TRSUP-CTA3-1", "TRV-CAC1-2", "TRE-CTC1-4", "TRL-CAA7-1",
 "TRQ-TTG9-1", "TRG-TCC2-2", "TRE-CTC1-3", "TRSUP-CTA2-1", "TRD-GTC2-5", "TRV-CAC1-3", "TRP-GGG1-1", "TRK-TTT14-1",
 "TRA-AGC15-1", "TRK-TTT8-1", "TRE-TTC15-1", "TRP-AGG5-1", "TRD-GTC2-3", "TRD-GTC2-4", "TRG-TCC2-4", "TRG-TCC2-3",
 "TRG-CCC7-1", "TRL-TAA5-1", "TRQ-TTG10-1", "TRA-TGC9-1", "TRC-GCA25-1", "TRS-ACT1-1", "TRQ-CTG9-1", "TRY-GTA11-1",
 "RN7SL471P", "RNU7-1", "RN7SK", "RN7SKP70", "RN7SKP253", "RN7SKP12", "RN7SKP123", "RN7SKP183", "RN7SKP233",
 "RN7SKP277", "RNU1-125P", "RNU2-63P", "RNU6-1189P", "RNU6-82P", "RNU6ATAC24P", "RN7SL4P", "RN7SL5P", "RN7SL128P",
 "RN7SL81P", "RN7SL359P", "RN7SL657P", "RN7SL731P", "RN7SL698P", "RN7SL674P", "RN7SL648P", "RN7SL396P", "RN7SL431P",
 "RN7SL192P", "RN7SL288P", "RN7SL368P", "RN7SL449P", "RN7SL145P", "RN7SL558P", "RN7SL246P", "RN7SKP255", "RN7SKP95",
 "RNU1-61P", "RNU1-108P", "RNU1-114P", "RNU1-86P", "RNU1-48P", "RNU1-107P", "RNU1-128P", "RNU1-40P", "RNU1-41P",
 "RNU1-95P", "RNU1-97P", "RNU2-57P", "RNU4-86P", "RNU5E-4P", "RNU6-9", "RNU6-7", "RNU6-781P", "RNU6-328P", "RNU6-1218P",
 "RNU6-457P", "RNU6-758P", "RNU6-194P", "RNU6-36P", "RNU6-458P", "RNU6-319P", "RNU6-34P", "RNU6-1077P", "RNU6-397P", 
 "RNU6-789P", "RNU6-437P", "RNU6-417P", "RNU6-842P", "RNU6-513P", "RNU6-1151P", "RNU6-830P", "RNU6-25P", "RNU6-1269P",
 "RNU6-312P", "RNU6-473P", "RNU6-368P", "RNU6-816P", "RNU6-1235P", "RNU6-631P", "RNU6-521P", "RNU6-941P", "RNU6-1314P",
 "RNU6-255P", "RNU6-303P", "RNU6-18P", "RNU6-1268P", "RNU6-749P", "RNU6-1239P", "RNU6-1071P", "RNU6-465P", "RNU6-785P",
 "RNU6-1076P", "RNU6-109P", "RNU6-1171P", "RNU6-184P", "RNY4P16", "RNY3P16", "TRV-CAC11-2", "RN7SL134P", "RN7SL52P",
 "RN7SL96P", "RN7SL440P", "RN7SL850P", "RN7SL787P", "RN7SL774P", "RN7SL429P", "RN7SL381P", "RN7SL759P", "RN7SL318P",
 "RN7SL275P", "RN7SL417P", "RN7SL584P", "RN7SL524P", "RN7SL123P", "RN7SL317P", "RN7SL92P", "RN7SL484P", "RN7SL525P",
 "RN7SL178P", "RN7SL592P", "RN7SL321P", "RN7SL753P", "RN7SL517P", "RN7SL311P", "RN7SL717P", "RN7SL331P", "RN7SL702P",
 "RN7SL722P", "RN7SL640P", "RN7SL656P", "RN7SL544P", "RN7SL199P", "RN7SL763P", "RN7SL725P", "RN7SL818P", "RN7SKP261",
 "RN7SKP159", "RN7SKP80", "RN7SKP176", "RN7SKP64", "RN7SKP111", "RN7SKP83", "RN7SKP283", "RN7SKP88", "RN7SKP126",
 "RN7SKP260", "RN7SKP49", "RN7SKP179", "RN7SKP131", "RN7SKP282"]

# ---------------------------------------------------------------------------
# Step 0a: Merge base conservation data with SNP enrichment data
# ---------------------------------------------------------------------------

print("Loading base conservation/annotation data...")
base_data = pd.read_csv(BASE_DATA_PATH)

print("Loading SNP enrichment data...")
snp_data = pd.read_csv(SNP_ENRICHMENT_PATH)

# Only pull in the enrichment-related columns from the SNP table so we don't
# clobber columns that already exist in the base table (e.g. gene_group).
snp_cols_to_merge = ['gene_name'] + [
    c for c in snp_data.columns
    if c.startswith('snp_density_')
    or c.startswith('enrichment_intergenic_')
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

# ---------------------------------------------------------------------------
# Step 0b: Merge in promoter motif data (type1/2/3 scores + motif counts)
# ---------------------------------------------------------------------------

print("Loading promoter motif data...")
promoter_data = pd.read_csv(PROMOTER_MOTIF_PATH)

# Only pull in the motif score/count columns so we don't duplicate columns
# already present (e.g. enrichment_intergenic_dbsnp is in both this table
# and the SNP enrichment table -- we keep the version already merged above).
promoter_cols_to_merge = ['gene_name'] + PROMOTER_FEATURES
promoter_data_subset = promoter_data[promoter_cols_to_merge]

data = data.merge(
    promoter_data_subset,
    left_on='Gene',
    right_on='gene_name',
    how='left'
).drop(columns=['gene_name'])

# The promoter motif columns come in with blank strings ("") where no motif
# of that type was found, so coerce to numeric first.
for col in PROMOTER_FEATURES:
    data[col] = pd.to_numeric(data[col], errors='coerce')

n_missing_promoter = data[PROMOTER_FEATURES].isna().any(axis=1).sum()
if n_missing_promoter:
    print(f"Note: {n_missing_promoter} genes have missing promoter motif "
          f"values (no motif of that type detected, or no match in the "
          f"promoter table). Filling missing motif scores/counts with 0.")
    data[PROMOTER_FEATURES] = data[PROMOTER_FEATURES].fillna(0)

data.to_csv(MERGED_DATA_PATH, index=False)
print(f"Saved merged dataset with SNP enrichment + promoter features to {MERGED_DATA_PATH}")

# ---------------------------------------------------------------------------
# Step 1: Labels + exclude ambiguous genes
# ---------------------------------------------------------------------------

data['label'] = data['Gene_Type'].apply(lambda x: 1 if x == 'Functional' else 0)

data_non_ambiguous = data[~data['Gene'].isin(ambiguous_genes)]

print("\nTotal counts in full dataset:")
print(f"Functional Genes: {data[data['label'] == 1].shape[0]}")
print(f"Pseudogenes: {data[data['label'] == 0].shape[0]}")

# ---------------------------------------------------------------------------
# Step 2: Train/test split (90/10) using conservation + enrichment + promoter features
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

# Feature importances so you can see how much the SNP enrichment and
# promoter motif features are contributing relative to the conservation scores.
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

    # Build the output table in the requested format
    ambiguous_output = pd.DataFrame({
        'HGNC Gene Symbol': ambiguous_data['Gene'],
        'Current Classification': ambiguous_data['Gene_Type'],
        'Predicted Classification': np.where(
            ambiguous_data['functional_probability'] >= 0.5,
            'Functional',
            'Pseudogene'
        ),
        'Random Forest Predicted Functional Probability': ambiguous_data['functional_probability'],
    })

    ambiguous_output.to_csv(AMBIGUOUS_PREDICTIONS_PATH, index=False)

    print("\nPredictions for ambiguous genes:")
    print(ambiguous_output)


# ---------------------------------------------------------------------------
# Step 6: Score the full dataset
# ---------------------------------------------------------------------------
 
print("\nAdding predicted probabilities for all genes in the full dataset...")
data['functional_probability'] = model.predict_proba(data[FEATURES])[:, 1]
 
# Re-save the combined dataset now that it includes the predicted
# functional probability for every gene.
data.to_csv(MERGED_DATA_PATH, index=False)
print(f"Updated {MERGED_DATA_PATH} with functional_probability column")
 
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
 
os.makedirs(PLOTS_DIR, exist_ok=True)
pdf = PdfPages(PLOTS_PDF_PATH)
 
 
def save_plot(filename):
    """Save the current figure as a PNG in PLOTS_DIR and append it to the
    combined PDF, then close the figure."""
    png_path = os.path.join(PLOTS_DIR, filename)
    plt.savefig(png_path, dpi=300)
    pdf.savefig()
    plt.close()
    print(f"Saved plot: {png_path}")
 
 
# --- Test set probability distribution ---
plt.figure(figsize=(10, 6))
sns.histplot(test_functional_probs, bins=20, color='firebrick', alpha=0.7, label='Functional Genes')
sns.histplot(test_pseudogene_probs, bins=20, color='cornflowerblue', alpha=0.7, label='Pseudogenes')
plt.xlabel('Probability of Being Functional')
plt.ylabel('Frequency')
plt.legend(title='Gene Type')
plt.grid(axis='y', linestyle='--', alpha=0.7)
plt.tight_layout()
save_plot('01_test_probability_distribution.png')
 
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
save_plot('02_train_vs_test_distribution.png')
 
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
save_plot('03_conservation_scatter.png')
 
# --- SNP enrichment scatter, colored by functional probability ---
plt.figure(figsize=(10, 6))
plt.scatter(
    ambiguous_data_plot['enrichment_intergenic_dbsnp'],
    ambiguous_data_plot['PhyloP100_median'],
    c=ambiguous_data_plot['functional_probability'],
    cmap='coolwarm', s=100, edgecolors='black'
)
plt.xlabel('dbSNP SNP Enrichment')
plt.ylabel('PhyloP100 Median')
plt.colorbar(label='Functional Probability')
plt.grid(False)
plt.tight_layout()
save_plot('04_snp_enrichment_scatter.png')
 
# --- Promoter motif scatter, colored by functional probability ---
# Total motif score (sum of type1/2/3 scores) vs. total motif count,
# to visualize how promoter motif strength/abundance separates predicted classes.
ambiguous_data_plot = ambiguous_data_plot.copy()
ambiguous_data_plot['total_motif_score'] = (
    ambiguous_data_plot['type1_score']
    + ambiguous_data_plot['type2_score']
    + ambiguous_data_plot['type3_score']
)
ambiguous_data_plot['total_motif_counts'] = (
    ambiguous_data_plot['type1_motif_counts']
    + ambiguous_data_plot['type2_motif_counts']
    + ambiguous_data_plot['type3_motif_counts']
)
 
plt.figure(figsize=(10, 6))
plt.scatter(
    ambiguous_data_plot['total_motif_score'],
    ambiguous_data_plot['total_motif_counts'],
    c=ambiguous_data_plot['functional_probability'],
    cmap='coolwarm', s=100, edgecolors='black'
)
plt.xlabel('Total Promoter Motif Score (Type 1+2+3)')
plt.ylabel('Total Promoter Motif Count (Type 1+2+3)')
plt.colorbar(label='Functional Probability')
plt.grid(False)
plt.tight_layout()
save_plot('05_promoter_motif_scatter.png')
 
pdf.close()
print(f"\nAll plots saved to {PLOTS_DIR}/ and combined into {PLOTS_PDF_PATH}")
