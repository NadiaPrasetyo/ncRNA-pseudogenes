# ncRNA-Derived Pseudogene Identification and Analysis

## Project Overview

This computational biology project aims to identify non-coding RNA-derived pseudogenes that are under strong negative selection within the human genome. By analyzing the conservation patterns of non-coding spliceosomal RNAs, we hope to reveal pseudogenes that may play important functional roles, despite their pseudogene classification. Identifying these pseudogenes could advance disease research by uncovering genetic elements that warrant further biological and clinical investigation. 

**Project Goals**:
1. Identify non-coding RNA-derived pseudogenes in the human genome from the HGNC database.
2. Analyze these pseudogenes for signs of strong negative selection.
3. Assess the potential biological relevance of these pseudogenes based on their evolutionary conservation.
4. Propose a reevaluation of the pseudogene status for promising candidates, especially in the context of disease-related research.

## Motivation

Pseudogenes derived from non-coding RNAs, particularly those related to the spliceosome, are often disregarded as non-functional. However, strong evolutionary conservation in these regions could indicate important regulatory or structural roles, especially if they exhibit negative selection. Identifying such elements may highlight new targets for research into genetic diseases and reveal overlooked aspects of genomic regulation. Recent studies have highligthed the significance of noncoding RNA RNU4-2 and RNU2-2P, causing neurodevelopmental disorders. The latter is especially interesting as RNU2-2P is classified as a non functional pesudogene of the functional splicesomoal RNA RNU2-1, annotated as so in many popular databases such as HGNC and UCSC. However, the evolutionary conservation tracks from UCSC show high levels of conservation of RNU2-2P, higher even the RNAU2-1. Indicating strong levels of negative selection are maintaining this sequence across the millions of years since it appeared in the genome. 

## Project Structure

-  Documentations
    - [Manuscript Draft](https://docs.google.com/document/d/1e3wGPWg_hJEjmLN69j1Wue_-8xi2F0YqOShb__H43ws/edit?usp=sharing)
    - [Lab Diary/ Notes](https://docs.google.com/document/d/1tBcWpbhzj_dht2hyIGH2wCKs0m5NIL1PiYrqAYb12Cg/edit?usp=sharing)
    - [Background information](https://docs.google.com/document/d/1dFvhfBuW2lBRTbdZu--Ip4hrldj980C-gZIP91_70iw/edit?usp=sharing)

- `data/`: Contains raw and processed datasets, including:
  - **Gene symbols chromosomal locations** (e.g., spliceosomal RNAs).
  - **Comparative genomics data** from other mammals for conservation analysis.

- `bin/`: Contains Python and R scripts for data processing, alignment, evolutionary analysis, and statistical tests:
  - `fetch_ncrna_hgnc.py`: Retrieves and processes HGNC symbols found from the HGNC database, mapped onto the HGNC annotations found in UCSC database.
  - `fetch_ncrna_data.py`: Retrieves and preprocesses HGNC symbols into chromosomal locations by mapping to the Ensembl database.
  - `fetch_conservation_data.py.`: Collects conservation data of each nucleotide position within the location range of each gene symbol (phastCons30way, phyloP100, and phyloP447).
  - `cleanup_txt_data`: Scans through temporary data files to look into no location found or no ensembl transcript IDs found cases and removes incorrect lines. Prints a list of gene symbols where the location was not found in UCSC
  - `fetch_expression_data.py`: Collects maximum expression data from the GTEX tracks downloaded from UCSC
  - `fetch_ENCODE_expr.py`: Collects maximum fpkm expression data from the ENCODE RNA sequence data, downloaded from the ENCODE RNA-Get portal
  - `get_specific_gene_list.py`: Create a list of functional genes with P and pseudogenes without P for exceptions of the rule: every gene that has a P in its gene symbol is a pseudogene
  - `random_forest_genes.py`: Creates a model that is trained on the difference between conservation and max expression data of each gene to calculate the probability of being functional. This is used to test a few ambiguos genes to calculate their functional probability.
  
  - `Boxplot_gene_*.R`: Analyzes conservation of identified pseudogenes across vertebrate species to detect signs of negative selection.
  - `Combined_table_numerical_features.R`: Combines all numerical features into one table
  - `Distplot_*_summary.R`: Creates a distribution histogram for each expression type - labelled for the top 2sd expression
  - `Distplot_gene_*.R`: Creates a distribution histogram for each conservation type - labelled for the top 2sd expression
  - `Jitterplot_*.R`: Calculates robust z-scores of each feature for each gene, normalised to the median of the pseudogenes and create a jitter-plot for the z-score of each feature - separated between functional ncRNA and pseudogene. Outliers are easily spotted with this plot
  - `Kolmogorov-Smirnov_or_wilcoxon_test.R`: Uses Kolmogorov Smirnov test for each gene group and pooled functional ncRNAs and pseudogenes for each conservation type
  - `KS_test_expression.R`: Uses Kolmogorov Smirnov test for each gene group and pooled functional ncRNAs and pseudogenes for each expression type
  - `ScatterPlot_PhyloP100_ENCODE.R`: Creates a scatter-plot of Encode maximum expression against PhyloP100 median. The encode max expression is scaled to improve visualisations.
  - `Violinplot_PhyloP100.R`: Calculates robust z-scores of phyloP100 median for each gene, normalised to the median of the pseudogenes and create a violin-plot for the z-score of each feature - separated between functional ncRNA and pseudogene. Outliers are indicated with jitter-plot overlay.

- `notebooks/`: Jupyter notebooks to explore and visualize the analysis results: not used too much - most documentation in lab diary and manuscript

- `results/`: Final results, including lists of pseudogenes under negative selection and visualizations of conservation patterns.

## Installation and Setup

To set up this project, you will need Python 3.x (Python 3.10.12 was used), R, and the following packages:

### Python Requirements
- `biomart 0.9.2`
- `mysql-connector-python .1.0`
- `numpy 2.1.3`
- `pyBigWig 0.3.23`
- `requests 2.32.3`

### R Requirements
- `ggplot2 3.5.1` 
- `readr 2.1.5` 
- `BiocManager 1.30.25`
- `dplyr 1.1.4`
- `ggrepel 0.9.6`
- `stringr 1.5.1`


Install the required Python packages using:
```bash
pip install biomart mysql-connector-python numpy pyBigWig requests
```

Install the required R packages:
```R
install.packages("ggplot2")
install.packages("BiocManager")
install.packages("readr")
```

### Data Requirements
1. **Comparative Genomics Data**: Available through UCSC Genome Browser’s PhastCons30way scores.

Place all downloaded datasets in the `data/` folder.

## Usage

1. **Identify Pseudogenes and Genes of a Specified Gene Symbol**: Run `fetch_ncrna_data.py` to identify gene symbols and pseudogenes associated with the specified gene, and their chromosomal locations.
   ```bash
   python bin/fetch_ncrna_data.py
   ```

2. **Data Clean up**: Run `cleanup_txt_data.py` to cleanup any invalid data that was received from Ensembl or HGNC. The script will pickup any lines that do not match a specified pattern and delete them automatically.
    ```bash
      python bin/cleanup_txt_data.py
      ```

3. **Analyze Conservation and Selection**: Run `fetch_conservation_data.py` to collect conservation data from the specified chromosomal locations defined in a txt file with a specific format. Run `R-script_plotting.R` to analyze conservation and detect signs of negative selection.
    ```bash
   python bin/fetch_conservation_data.py
   ```
   ```R
   Rscript bin/R/R-script_plotting.R
   ```

4. **Visualize Results**: Open `results_summary.ipynb` in Jupyter Notebook to view and interpret the results. @TODO

## Methods

### Pseudogene Identification
The project computationaly analyses the pseudogene annotations in HGNC along with conservation data to identify likely functional pseudogenes, indicating false annotation and clinical and diagnosis importance of these pseudogenes.

### Cleaning pseudogene symbols found
Crosscheking with the HGNC database to check if symbols found are approved or withdrawn. Withdrawn symbols should be removed from the data in cleanup, while some found approved symbols do not have an associated ensembl transcript ID, must be added in manually by cross checking with UCSC genome browser. Cases like this include: RNU4-3P, RNU5B-5P, RNU5E-2P, RNU5F-5P, RNU6-52P, RNU6-69P, RNU6-1139P.

### Conservation and Selection Analysis
We use comparative genomics data to assess conservation levels across mammalian species, applying tools such as PhastCons scores to detect signs of negative selection. Negatively selected regions are likely to hold functional significance despite being labeled as pseudogenes.

## Results and Interpretation

The final output includes:
- **Lists of conserved pseudogenes** that may be under negative selection.
- **Conservation plots** showing evolutionary patterns across species.
- **Functional hypothesis** for pseudogenes identified as potentially relevant in disease research.

## Future Directions

This pipeline can be expanded in the following ways:
- **Functional Assays**: Experimentally validate the potential functions of identified pseudogenes.
- **Disease Association Analysis**: Investigate links between these pseudogenes and disease phenotypes in genome-wide association studies (GWAS).

## Authors and Acknowledgments

Developed by the Gardner lab - University of Otago, with contributions from Nadia Prasetyo.

## License

This project is licensed under the MIT License.
