import pandas as pd
import numpy as np
from sklearn.ensemble import RandomForestClassifier
from sklearn.model_selection import train_test_split
from sklearn.metrics import classification_report, roc_auc_score
import matplotlib.pyplot as plt
import seaborn as sns

# Load the dataset
data = pd.read_csv('results/combined_gene_data.csv')

# Extract relevant features and labels
data['label'] = data['Gene_Type'].apply(lambda x: 1 if x == 'Functional' else 0)

# List of ambiguous genes (assuming you have them defined)
# You can modify this to use any criteria you have for ambiguous genes
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

# Step 1: Exclude ambiguous genes from the dataset
data_non_ambiguous = data[~data['Gene'].isin(ambiguous_genes)]

# Diagnostics: Print number of functional and pseudogenes in the full dataset
print("Total counts in full dataset:")
print(f"Functional Genes: {data[data['label'] == 1].shape[0]}")
print(f"Pseudogenes: {data[data['label'] == 0].shape[0]}")

# Step 2: Split the non-ambiguous data into training and test sets (90% training, 10% testing)
X = data_non_ambiguous[['PhyloP100_median', 'ENCODE_max']]  # Features
y = data_non_ambiguous['label']  # Labels

# Randomly split into 90% training and 10% test
X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.1, random_state=42)

# Diagnostics: Print the number of functional and pseudogenes in the training and test sets
print("\nTraining Set Counts:")
print(f"Functional Genes: {y_train.sum()}")
print(f"Pseudogenes: {len(y_train) - y_train.sum()}")

print("\nTest Set Counts:")
print(f"Functional Genes: {y_test.sum()}")
print(f"Pseudogenes: {len(y_test) - y_test.sum()}")

# Step 3: Initialize the Random Forest Classifier and train the model
model = RandomForestClassifier(n_estimators=100, random_state=42)
print("\nTraining the Random Forest model...")
model.fit(X_train, y_train)

# Step 4: Evaluate the model on the test set
print("\nEvaluating the model on the test set...")
test_predictions = model.predict(X_test)
test_probs = model.predict_proba(X_test)[:, 1]  # Probabilities for functional genes

# Print the classification report and ROC AUC score for the test set
print("\nTest Classification Report:")
print(classification_report(y_test, test_predictions))
print(f"Test ROC AUC Score: {roc_auc_score(y_test, test_probs):.2f}")

# Step 5: Predict on ambiguous genes
# Select ambiguous genes from the original dataset
ambiguous_data = data[data['Gene'].isin(ambiguous_genes)]

# Diagnostics: Check if ambiguous data is empty
if ambiguous_data.empty:
    print("\nNo ambiguous genes found in the dataset.")
else:
    # Print the number of functional and pseudogenes in ambiguous genes
    print("\nAmbiguous Genes Counts:")
    print(f"Functional Genes: {ambiguous_data[ambiguous_data['label'] == 1].shape[0]}")
    print(f"Pseudogenes: {ambiguous_data[ambiguous_data['label'] == 0].shape[0]}")

    # Predict probabilities for ambiguous genes
    print("\nPredicting probabilities for ambiguous genes...")
    X_ambiguous = ambiguous_data[['PhyloP100_median', 'ENCODE_max']]
    ambiguous_probs = model.predict_proba(X_ambiguous)[:, 1]

    # Add the predicted probabilities to the ambiguous data
    ambiguous_data['functional_probability'] = ambiguous_probs

    # Save the ambiguous gene predictions to a CSV file
    ambiguous_data[['Gene', 'functional_probability']].to_csv('results/ambiguous_gene_predictions.csv', index=False)

    # Optionally, print ambiguous gene predictions for review
    print("\nPredictions for ambiguous genes:")
    print(ambiguous_data[['Gene', 'functional_probability']])

# Step 6: Add the predicted probabilities for all genes in the full dataset
print("\nAdding predicted probabilities for all genes in the full dataset...")
data['functional_probability'] = model.predict_proba(data[['PhyloP100_median', 'ENCODE_max']])[:, 1]


# Visualizations for Test Data Distribution

# Predicted probabilities for test data
test_functional_probs = test_probs[y_test == 1]  # Probabilities for true functional genes
test_pseudogene_probs = test_probs[y_test == 0]  # Probabilities for true pseudogenes

# Set font to Times New Roman
plt.rcParams['font.family'] = 'DejaVu Serif'
# Set global font size for everything in the plot
plt.rcParams.update({'font.size': 28,  # Global font size for all text
                     'axes.labelsize': 28,  # Axis labels font size
                     'xtick.labelsize': 28,  # X-axis tick label font size
                     'ytick.labelsize': 28,  # Y-axis tick label font size
                     'legend.fontsize': 22,  # Legend font size
                     'figure.titlesize': 28})  # Global figure title font size

# Create the histogram for test data probabilities
plt.figure(figsize=(10, 6))
sns.histplot(test_functional_probs, bins=20, color='firebrick', alpha=0.7, label='Functional Genes')
sns.histplot(test_pseudogene_probs, bins=20, color='cornflowerblue', alpha=0.7, label='Pseudogenes')

# Add titles and labels
plt.xlabel('Probability of Being Functional')
plt.ylabel('Frequency')
plt.legend(title='Gene Type')
plt.grid(axis='y', linestyle='--', alpha=0.7)

# Show the plot
plt.tight_layout()
plt.show()

# Predicted probabilities for training data
train_probs = model.predict_proba(X_train)[:, 1]
train_functional_probs = train_probs[y_train == 1]
train_pseudogene_probs = train_probs[y_train == 0]

# Overlay training data distributions
plt.figure(figsize=(10, 6))
sns.histplot(test_functional_probs, bins=20, color='firebrick', alpha=0.7, label='Test Functional Genes')
sns.histplot(test_pseudogene_probs, bins=20, color='cornflowerblue', alpha=0.7, label='Test Pseudogenes')
sns.histplot(train_functional_probs, bins=20, color='lightpink', alpha=0.5, label='Train Functional Genes', linestyle='--')
sns.histplot(train_pseudogene_probs, bins=20, color='lightblue', alpha=0.5, label='Train Pseudogenes', linestyle='--')

# Add titles and labels
plt.xlabel('Probability of Being Functional')
plt.ylabel('Frequency')
plt.legend(title='Gene Type')
plt.grid(axis='y', linestyle='--', alpha=0.7)

# Show the plot
plt.tight_layout()
plt.show()


# Scatter plot of ambiguous gene probabilities
plt.figure(figsize=(10, 6))

# Filter out ambiguous data
ambiguous_data = data[data['Gene'].isin(ambiguous_genes)]

# Scatter plot for ambiguous genes
plt.scatter(ambiguous_data['PhyloP100_median'], ambiguous_data['ENCODE_max'], 
            c=ambiguous_data['functional_probability'], cmap='coolwarm', s=100, edgecolors='black')


# Add titles and labels
plt.xlabel('PhyloP100 Median')
plt.ylabel('ENCODE Max')
plt.colorbar(label='Functional Probability')

# Disable grid for better clarity
plt.grid(False)

# Show the plot
plt.tight_layout()
plt.show()