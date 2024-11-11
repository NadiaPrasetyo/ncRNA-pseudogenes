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

- `data/`: Contains raw and processed datasets, including:
  - **Gene symbols chromosomal locations** (e.g., spliceosomal RNAs).
  - **Comparative genomics data** from other mammals for conservation analysis.

- `bin/`: Contains Python and R scripts for data processing, alignment, evolutionary analysis, and statistical tests:
  - `fetch_ncrna_data.py`: Retrieves and preprocesses HGNC symbols into chromosomal locations by mapping to the Ensembl database.
  - `fetch_conservation_data.py.`: Collects conservation data of each nucleotide position within the location range of each gene symbol (phastCons30way).
  - `R-script_plotting.R`: Analyzes conservation of identified pseudogenes across vertebrate species to detect signs of negative selection.

- `notebooks/`: Jupyter notebooks to explore and visualize the analysis results: @TODO

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
- `ggplot2` (for visualizations)
- `readr` 
- `BiocManager`


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

2. **Analyze Conservation and Selection**: Run `fetch_conservation_data.py` to collect conservation data from the specified chromosomal locations defined in a txt file with a specific format. Run `R-script_plotting.R` to analyze conservation and detect signs of negative selection.
    ```bash
   python bin/fetch_conservation_data.py
   ```
   ```R
   Rscript bin/R/R-script_plotting.R
   ```

3. **Visualize Results**: Open `results_summary.ipynb` in Jupyter Notebook to view and interpret the results. @TODO

## Methods

### Pseudogene Identification
The project computationaly analyses the pseudogene annotations in HGNC along with conservation data to identify likely functional pseudogenes, indicating false annotation and clinical and diagnosis importance of these pseudogenes.

### Cleaning pseudogene symbols found
Crosscheking with the HGNC database to check if symbols found are approved or withdrawn. Withdrawn symbols should be removed from the data, while some found approved symbols do not have an associated ensembl transcript ID, must be added in manually by cross checking with UCSC genome browser. Cases like this include: RNU4-3P, RNU5B-5P, RNU5E-2P, RNU5F-5P, RNU6-52P, RNU6-69P.

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
