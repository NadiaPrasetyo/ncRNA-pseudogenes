ncRNA-Derived Pseudogene Identification and Analysis
Project Overview

This computational biology project aims to identify non-coding RNA-derived pseudogenes that are under strong negative selection within the human genome. By analyzing the conservation patterns of non-coding spliceosomal RNAs, we hope to reveal pseudogenes that may play important functional roles, despite their pseudogene classification. Identifying these pseudogenes could advance disease research by uncovering genetic elements that warrant further biological and clinical investigation.

Project Goals:

    Identify non-coding RNA-derived pseudogenes in the human genome.
    Analyze these pseudogenes for signs of strong negative selection.
    Assess the potential biological relevance of these pseudogenes based on their evolutionary conservation.
    Propose a reevaluation of the pseudogene status for promising candidates, especially in the context of disease-related research.

Motivation

Pseudogenes derived from non-coding RNAs, particularly those related to the spliceosome, are often disregarded as non-functional. However, strong evolutionary conservation in these regions could indicate important regulatory or structural roles, especially if they exhibit negative selection. Identifying such elements may highlight new targets for research into genetic diseases and reveal overlooked aspects of genomic regulation.
Project Structure

    data/: Contains raw and processed datasets, including:
        Non-coding RNA sequences (e.g., spliceosomal RNAs).
        Genomic annotations for human pseudogenes.
        Comparative genomics data from other vertebrates for conservation analysis.

    scripts/: Contains Python and R scripts for data processing, alignment, evolutionary analysis, and statistical tests:
        fetch_ncrna_data.py: Retrieves and preprocesses non-coding RNA sequences.
        identify_pseudogenes.py: Maps non-coding RNA sequences onto the human genome to identify putative pseudogene regions.
        conservation_analysis.R: Analyzes conservation of identified pseudogenes across vertebrate species to detect signs of negative selection.

    notebooks/: Jupyter notebooks to explore and visualize the analysis results:
        exploratory_analysis.ipynb: Visual exploration of pseudogene features and initial conservation patterns.
        results_summary.ipynb: Summary and visualization of negatively selected pseudogenes with potential functional relevance.

    results/: Final results, including lists of pseudogenes under negative selection and visualizations of conservation patterns.

Installation and Setup

To set up this project, you will need Python 3.x, R, and the following packages:
Python Requirements

    biopython
    pandas
    numpy
    scipy
    pybedtools

R Requirements

    ape (for phylogenetic analysis)
    ggplot2 (for visualizations)

Install the required Python packages using:

pip install biopython pandas numpy scipy pybedtools

Install the required R packages:

install.packages("ape")
install.packages("ggplot2")

Data Requirements

    Genome Annotation Data: Download from Ensembl or UCSC Genome Browser.
    Comparative Genomics Data: Available through UCSC Genome Browser’s PhyloP/PhastCons scores.

Place all downloaded datasets in the data/ folder.
Usage

    Identify Putative Pseudogenes: Run identify_pseudogenes.py to identify potential pseudogenes derived from spliceosomal RNAs.

python scripts/identify_pseudogenes.py

Analyze Conservation and Selection: Run conservation_analysis.R to analyze conservation and detect signs of negative selection.

    Rscript scripts/conservation_analysis.R

    Visualize Results: Open results_summary.ipynb in Jupyter Notebook to view and interpret the results.

Methods
Pseudogene Identification

The project identifies regions in the human genome that resemble known non-coding RNAs (specifically spliceosomal RNAs) but lack protein-coding potential. These regions are then matched to annotated pseudogenes in the human genome.
Conservation and Selection Analysis

We use comparative genomics data to assess conservation levels across vertebrate species, applying tools such as PhyloP and PhastCons scores to detect signs of negative selection. Negatively selected regions are likely to hold functional significance despite being labeled as pseudogenes.
Results and Interpretation

The final output includes:

    Lists of conserved pseudogenes that may be under negative selection.
    Conservation plots showing evolutionary patterns across species.
    Functional hypothesis for pseudogenes identified as potentially relevant in disease research.

Future Directions

This pipeline can be expanded in the following ways:

    Functional Assays: Experimentally validate the potential functions of identified pseudogenes.
    Disease Association Analysis: Investigate links between these pseudogenes and disease phenotypes in genome-wide association studies (GWAS).
    Broader Genomic Scope: Extend this analysis to include other classes of non-coding RNAs (e.g., miRNAs, lncRNAs).

Authors and Acknowledgments

Developed by the [Your Research Group/Institution Name], with contributions from [Team Members’ Names].

For inquiries or collaboration, please contact [Your Email].
License

This project is licensed under the MIT License.
