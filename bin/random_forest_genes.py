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

# List of test genes
test_genes_list = ['RNU6-1189P', 'RNU6-82P', 'RNU2-2P', 'RNU1-27P', 'RNU1-28P', 'RNU5B-1']

# Separate training and test data
train_data = data[~data['Gene'].isin(test_genes_list)]  # Exclude test genes from training
test_data = data[data['Gene'].isin(test_genes_list)]  # Include only test genes

# Features for training and testing
X_train = train_data[['PhyloP100_median', 'ENCODE_max']]
y_train = train_data['label']
X_test = test_data[['PhyloP100_median', 'ENCODE_max']]

# Initialize the Random Forest Classifier
model = RandomForestClassifier(n_estimators=100, random_state=42)

# Train the model
model.fit(X_train, y_train)

# Predict probabilities for the test data
test_probabilities = model.predict_proba(X_test)[:, 1]  # Probabilities for functional genes

# Save test predictions
test_data.loc[:, 'functional_probability'] = test_probabilities
test_data[['Gene', 'functional_probability']].to_csv('test_predictions.csv', index=False)

# Evaluate the model on the training data
train_predictions = model.predict(X_train)
train_probs = model.predict_proba(X_train)[:, 1]
print("Training Classification Report:")
print(classification_report(y_train, train_predictions))
print(f"Training ROC AUC Score: {roc_auc_score(y_train, train_probs):.2f}")

print("Test predictions have been saved to 'results/test_predictions.csv'")

# Add the predicted probabilities to the full dataset
data['functional_probability'] = model.predict_proba(data[['PhyloP100_median', 'ENCODE_max']])[:, 1]

# Split the data by gene type for visualization
functional_genes = data[data['label'] == 1]['functional_probability']
pseudogenes = data[data['label'] == 0]['functional_probability']


# Set font to Times New Roman
plt.rcParams['font.family'] = 'DejaVu Serif'
# Set global font size for everything in the plot
plt.rcParams.update({'font.size': 32,  # Global font size for all text
                     'axes.labelsize': 32,  # Axis labels font size
                     'xtick.labelsize': 32,  # X-axis tick label font size
                     'ytick.labelsize': 32,  # Y-axis tick label font size
                     'legend.fontsize': 22,  # Legend font size
                     'figure.titlesize': 32})  # Global figure title font size

# Create the bar plot
plt.figure(figsize=(10, 6))
sns.histplot(functional_genes, bins=20, color='firebrick', alpha=0.7, label='Functional Genes')
sns.histplot(pseudogenes, bins=20, color='cornflowerblue', alpha=0.7, label='Pseudogenes')

# Add titles and labels
plt.xlabel('Probability of Being Functional')
plt.ylabel('Frequency')
plt.legend(title='Gene Type')
plt.grid(axis='y', linestyle='--', alpha=0.7)

# Show the plot
plt.tight_layout()
plt.show()