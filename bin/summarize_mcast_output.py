import argparse
import csv
import pandas as pd
import logging
import os

def compile_mcast_motifs(xml_file):
    """
    Compile MCAST motifs from the given XML file into a pandas DataFrame.
    """

# e.g. line: <motif name="1-WRGGATGACABG" width="12" best-possible-match="AAGGATGACAGG" best-possible-rc-match="CCTGTCATCCTT"/>
    motif_names = []

    for line in xml_file:
        if line.startswith("<motif"):
            # Extract attributes motif name
            name_start = line.find('name="') + len('name="')
            name_end = line.find('"', name_start)
            motif_name = line[name_start:name_end]
            motif_names.append(motif_name)
    return motif_names

def extract_mcast_hits(xml_file, motif_names):
    """
    Extract MCAST hits from the given XML file into a pandas DataFrame.
    """

#     e.g. line:
# <program-name>mcast</program-name>
# <parameters>
# <pattern-file>data/pol3_fasta_by_class/I_streme_out/streme.txt</pattern-file>
# <sequence-file>data/hgnc_gene_seq.fasta</sequence-file>
# <pattern-pvalue-cutoff>0.0005</pattern-pvalue-cutoff>
# <sequence-pvalue-cutoff>1</sequence-pvalue-cutoff>
# </parameters>
# <multi-pattern-scan score="51.8165" pvalue="8.1151e-12">
# <pattern accession="7-CCTGGGAATACCGG" name="7-CCTGGGAATACCGG">
# <scanned-sequence accession="RN7SL690P|pseudogene|chr20" name="RN7SL690P|pseudogene|chr20">
# <matched-element start="20316769" stop="20316756" pvalue="2.2044e-07">
# <sequence>CCTGTAATCCCAGC</sequence>
# </matched-element>
# </scanned-sequence>
# </pattern>
# <pattern accession="11-CTCCTGACCTC" name="11-CTCCTGACCTC">
# <scanned-sequence accession="RN7SL690P|pseudogene|chr20" name="RN7SL690P|pseudogene|chr20">
# <matched-element start="20316811" stop="20316801" pvalue="5.7689e-05">
# <sequence>GAACCCAGGAG</sequence>
# </matched-element>
# </scanned-sequence>
# </pattern>
# <pattern accession="6-CYGAGATCA" name="6-CYGAGATCA">
# <scanned-sequence accession="RN7SL690P|pseudogene|chr20" name="RN7SL690P|pseudogene|chr20">
# <matched-element start="20316829" stop="20316837" pvalue="8.8096e-05">
# <sequence>CCAAGATCG</sequence>
# </matched-element>
# </scanned-sequence>
# </pattern>
# <pattern accession="10-ACTCCAGCCTGG" name="10-ACTCCAGCCTGG">
# <scanned-sequence accession="RN7SL690P|pseudogene|chr20" name="RN7SL690P|pseudogene|chr20">
# <matched-element start="20316847" stop="20316858" pvalue="3.792e-08">
# <sequence>ACTCCAGCCTGG</sequence>
# </matched-element>
# </scanned-sequence>
# </pattern>
# <pattern accession="5-AAGCTAAGCAGGGTC" name="5-AAGCTAAGCAGGGTC">
# <scanned-sequence accession="RN7SL690P|pseudogene|chr20" name="RN7SL690P|pseudogene|chr20">
# <matched-element start="20316891" stop="20316905" pvalue="0.00047381">
# <sequence>AAGGCAGACAGAGTC</sequence>
# </matched-element>
# </scanned-sequence>
# </pattern>
# <pattern accession="10-ACTCCAGCCTGG" name="10-ACTCCAGCCTGG">
# <scanned-sequence accession="RN7SL690P|pseudogene|chr20" name="RN7SL690P|pseudogene|chr20">
# <matched-element start="20316929" stop="20316918" pvalue="3.792e-08">
# <sequence>CCAGGCTGGAGT</sequence>
# </matched-element>
# </scanned-sequence>
# </pattern>
# <pattern accession="6-CYGAGATCA" name="6-CYGAGATCA">
# <scanned-sequence accession="RN7SL690P|pseudogene|chr20" name="RN7SL690P|pseudogene|chr20">
# <matched-element start="20316947" stop="20316939" pvalue="8.8096e-05">
# <sequence>CGATCTTGG</sequence>
# </matched-element>
# </scanned-sequence>
# </pattern>
# <pattern accession="11-CTCCTGACCTC" name="11-CTCCTGACCTC">
# <scanned-sequence accession="RN7SL690P|pseudogene|chr20" name="RN7SL690P|pseudogene|chr20">
# <matched-element start="20316965" stop="20316975" pvalue="5.7689e-05">
# <sequence>CTCCTGGGTTC</sequence>
# </matched-element>
# </scanned-sequence>
# </pattern>
# <pattern accession="7-CCTGGGAATACCGG" name="7-CCTGGGAATACCGG">
# <scanned-sequence accession="RN7SL690P|pseudogene|chr20" name="RN7SL690P|pseudogene|chr20">
# <matched-element start="20317007" stop="20317020" pvalue="2.2044e-07">
# <sequence>GCTGGGATTACAGG</sequence>
# </matched-element>
# </scanned-sequence>
# </pattern>
# <pattern accession="7-CCTGGGAATACCGG" name="7-CCTGGGAATACCGG">
# <scanned-sequence accession="RN7SL690P|pseudogene|chr20" name="RN7SL690P|pseudogene|chr20">
# <matched-element start="20317040" stop="20317027" pvalue="0.00022886">
# <sequence>CCAGCATGCCCGGC</sequence>
# </matched-element>
# </scanned-sequence>
# </pattern>
# <mem:match cluster-id="cluster-1870" seq-name="RN7SL690P|pseudogene|chr20" start="20316756" stop="20317040" evalue="2.6285e-08" qvalue="1.0999e-08">CCTGTAATCCCAGCTACTCGGGAGGCTGAGGCATGATAATCGCTTGAACCCAGGAGGCAGAGTTTGCAGTGAGCCAAGATCGTGCCACTGCACTCCAGCCTGGATGACAGAGTGAGACTCTGTCTGAAAAAAACCAAGGCAGACAGAGTCTCACTCTGTCATCCAGGCTGGAGTGCAGTGGCACGATCTTGGCTCACTGCAACCTCTGCCTCCTGGGTTCAAGCGATTATCATGCCTCAGCTTCCCGAGTAGCTGGGATTACAGGCTTGAGCCAGCATGCCCGGC
# </mem:match>
# </multi-pattern-scan>
# <multi-pattern-scan score="50.6701" pvalue="2.3177e-11">
# <pattern accession="2-GAYGAGATCGGGCRC" name="2-GAYGAGATCGGGCRC">
# <scanned-sequence accession="TSEN15P1|pseudogene|chr17" name="TSEN15P1|pseudogene|chr17">
# <matched-element start="17456573" stop="17456587" pvalue="8.5343e-05">
# <sequence>GATGAAATTAGGCTG</sequence>

    summary_df = pd.DataFrame(columns=["sequence_name", "cluster_id", "score", "p-value"] + motif_names)
    # for each motif, we will create a column in the DataFrame to store the number of hits for that motif in each sequence
    
    with open(xml_file, "r") as f:
        current_seq_name = None
        current_cluster = None
        current_score = None
        current_pvalue = None
        current_hits = {motif: 0 for motif in motif_names}

        for line in f:
            line = line.strip()
            # separate the file into blocks starting with <multi-pattern-scan> and ending with </multi-pattern-scan>
            if line.startswith("<multi-pattern-scan"):
                # Extract score and pvalue
                score_start = line.find('score="') + len('score="')
                score_end = line.find('"', score_start)
                current_score = float(line[score_start:score_end])

                pvalue_start = line.find('pvalue="') + len('pvalue="')
                pvalue_end = line.find('"', pvalue_start)
                current_pvalue = float(line[pvalue_start:pvalue_end])

                # If we were already processing a sequence, save its data before starting a new one
                if current_seq_name is not None:
                    summary_df.loc[len(summary_df)] = {
                        "sequence_name": current_seq_name,
                        "cluster_id": current_cluster,
                        "score": current_score,
                        "p-value": current_pvalue,
                        **current_hits,
                    }
                    # Reset hits for the new sequence
                    current_hits = {motif: 0 for motif in motif_names}
            
            # the seq name is in: # <mem:match cluster-id="cluster-1870" seq-name="RN7SL690P|pseudogene|chr20" start="20316756" stop="20317040" evalue="2.6285e-08" qvalue="1.0999e-08">CCTGTAATCCCAGCTACTCGGGAGGCTGAGGCATGATAATCGCTTGAACCCAGGAGGCAGAGTTTGCAGTGAGCCAAGATCGTGCCACTGCACTCCAGCCTGGATGACAGAGTGAGACTCTGTCTGAAAAAAACCAAGGCAGACAGAGTCTCACTCTGTCATCCAGGCTGGAGTGCAGTGGCACGATCTTGGCTCACTGCAACCTCTGCCTCCTGGGTTCAAGCGATTATCATGCCTCAGCTTCCCGAGTAGCTGGGATTACAGGCTTGAGCCAGCATGCCCGGC
            elif line.startswith("<mem:match"):
                seq_name_start = line.find('seq-name="') + len('seq-name="')
                seq_name_end = line.find('"', seq_name_start)
                current_seq_name = line[seq_name_start:seq_name_end]
                cluster_name_start = line.find('cluster-id="') + len('cluster-id="')
                cluster_name_end = line.find('"', cluster_name_start)
                current_cluster = line[cluster_name_start:cluster_name_end]
            
            # motif matches are in: <pattern accession="7-CCTGGGAATACCGG" name="7-CCTGGGAATACCGG">
            elif line.startswith("<pattern"):
                motif_name_start = line.find('name="') + len('name="')
                motif_name_end = line.find('"', motif_name_start)
                current_motif = line[motif_name_start:motif_name_end]
                if current_motif in current_hits:
                    current_hits[current_motif] += 1

        # After the loop, make sure to save the last sequence's data
        if current_seq_name is not None:
            summary_df.loc[len(summary_df)] = {
                "sequence_name": current_seq_name,
                "cluster_id": current_cluster,
                "score": current_score,
                "p-value": current_pvalue,
                **current_hits,
            }

        return summary_df

def expected_promoter(gene_name):
    """
    Determine the expected promoter type for a given gene name.
    """
    # key: Class of RNA:
    # Type 1: 5s RNA
    # Type 2: TRNA, RN7SL RNA, VTRNA
    # Type 3: RNU6, RNY, RN7SK, RNU6ATAC, tRNASeC
    if gene_name.upper().startswith("5S"):
        return "Type 1"
    elif gene_name.upper().startswith("TRNA") or gene_name.upper().startswith("RN7SL") or gene_name.upper().startswith("VTRNA"):
        return "Type 2"
    elif gene_name.upper().startswith("RNU6") or gene_name.upper().startswith("RNY") or gene_name.upper().startswith("RN7SK") or gene_name.upper().startswith("RNU6ATAC") or gene_name.upper().startswith("TRNASEC"):
        return "Type 3"
    else:
        return "None"

def summarize_types(mcast_input_dir, gene_names_file):
    """
    Summarize MCAST output from the given directory into a pandas DataFrame.
    """
    # separate the three polIII classes based on the subdirectory names, e.g. I_mcast_out, II_mcast_out, III_mcast_out
    classI_dir = os.path.join(mcast_input_dir, "I_mcast_out")
    classII_dir = os.path.join(mcast_input_dir, "II_mcast_out")
    classIII_dir = os.path.join(mcast_input_dir, "III_mcast_out")

    #extract motif names and hits from each class using files: cisml.xml (hits) and mcast.xml (motifs)
    classI_motifs = compile_mcast_motifs(os.path.join(classI_dir, "mcast.xml"))
    classI_hits = extract_mcast_hits(os.path.join(classI_dir, "cisml.xml"), classI_motifs)
    classII_motifs = compile_mcast_motifs(os.path.join(classII_dir, "mcast.xml"))
    classII_hits = extract_mcast_hits(os.path.join(classII_dir, "cisml.xml"), classII_motifs)
    classIII_motifs = compile_mcast_motifs(os.path.join(classIII_dir, "mcast.xml"))
    classIII_hits = extract_mcast_hits(os.path.join(classIII_dir, "cisml.xml"), classIII_motifs)

    # Combine the three classes into a single summary dataframe: 
    # for each gene_name in gene_names_file, we will create a true or false to having any hits in each class
    summary_df = pd.DataFrame(columns=["gene_name", "expected_promoter", "contains_type1_promoter", "contains_type2_promoter", "contains_type3_promoter"])

    # gene_names_file may contain additional comma-separated annotations. Use
    # only the Gene column, and handle exports that quote the complete row.
    with open(gene_names_file, "r", newline="") as f:
        gene_names = []
        for row in csv.reader(f):
            if not row:
                continue
            gene_name = row[0].strip()
            if "," in gene_name:
                gene_name = gene_name.split(",", 1)[0].strip()
            if gene_name.lower() == "gene":
                continue
            gene_names.append(gene_name)

    for gene_name in gene_names:
        # the gene names will only have the first part of the sequence name, e.g. a hit sequence_name of RN7SL690P|pseudogene|chr20 will match RN7SL690P
        contains_type1_promoter = any(gene_name in seq_name for seq_name in classI_hits["sequence_name"]) #if the gene_name is in any of the sequence names in the classI_hits dataframe, then it contains a type1 promoter
        contains_type2_promoter = any(gene_name in seq_name for seq_name in classII_hits["sequence_name"]) #if the gene_name is in any of the sequence names in the classII_hits dataframe, then it contains a type2 promoter
        contains_type3_promoter = any(gene_name in seq_name for seq_name in classIII_hits["sequence_name"]) #if the gene_name is in any of the sequence names in the classIII_hits dataframe, then it contains a type3 promoter

        summary_df.loc[len(summary_df)] = {
            "gene_name": gene_name,
            "expected_promoter": expected_promoter(gene_name),
            "contains_type1_promoter": contains_type1_promoter,
            "contains_type2_promoter": contains_type2_promoter,
            "contains_type3_promoter": contains_type3_promoter
        }

    return summary_df

def main():
    parser = argparse.ArgumentParser(description="Summarize MCAST output")
    parser.add_argument("-d", "--mcast_input_dir", help="Path to the MCAST input directory")
    parser.add_argument("-i", "--gene_names_file", help="Path to the gene names file")
    parser.add_argument("-o", "--output_csv", default="results/mcast_summary.csv", help="Path to the output CSV file")
    args = parser.parse_args()

    logging.basicConfig(level=logging.INFO)

    if not os.path.exists(args.mcast_input_dir):
        logging.error(f"MCAST input directory {args.mcast_input_dir} does not exist.")
        return

    if not os.path.exists(args.gene_names_file):
        logging.error(f"Gene names file {args.gene_names_file} does not exist.")
        return

    if not os.path.exists(os.path.dirname(args.output_csv)):
        logging.error(f"Output directory for CSV file {args.output_csv} does not exist.")
        return

    summary_df = summarize_types(args.mcast_input_dir, args.gene_names_file)
    summary_df.to_csv(args.output_csv, index=False)

if __name__ == "__main__":
    main()