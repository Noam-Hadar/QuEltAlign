#!/usr/bin/env python
# coding: utf-8

# In[ ]:


# Install necessary packages if not already installed
# pip3 install pandas primer3-py pyfaidx biopython networkx matplotlib numpy scipy

# Standard library imports
import os
import random
import pickle
from glob import glob
from time import sleep
from itertools import product
import xml.etree.ElementTree as ET

# External library imports
import pandas as pd
import primer3
from pyfaidx import Fasta
from Bio.Blast import NCBIWWW, NCBIXML
import networkx as nx
import matplotlib.pyplot as plt
import numpy as np

# Define constants
fastq = 'G6PD.fastq'  # Input FASTQ file name
quality = 0.1  # Minimal quality Phred score
min_size = 20  # Minimum sequence length after trimming
min_overlap = 50  # Minimal overlap for de novo assembly
runs = 10

# Phred quality score mapping
phredict = {
    c: 10 ** (-(ord(c) - 33) / 10.0) for c in 
    "!\"#$%&'()*+,-./0123456789:;<=>?@ABCDEFGHIJK"
}

# Function to extract high-quality sequences
def extract_high_quality_sequences(fastq_file, quality_threshold, min_size):
    sequences = []
    with open(fastq_file, 'r', encoding='utf-8-sig') as f:
        lines = f.readlines()
        for i in range(0, len(lines), 4):
            current_sequence = lines[i + 1].strip()
            quality_scores = lines[i + 3].strip()
            #phred_scores = [phredict[char] for char in quality_scores]
            phred_scores = [phredict['F'] for char in quality_scores]

            # Filter low-quality bases
            high_quality_sequence = ''.join(
                [
                    current_sequence[j]
                    for j in range(len(current_sequence))
                    if phred_scores[j] <= quality_threshold
                ]
            )

            # Check minimum size constraint
            if len(high_quality_sequence) >= min_size:
                sequences.append(high_quality_sequence)

    return sequences

# Extract sequences and remove duplicates
high_quality_sequences = extract_high_quality_sequences(fastq, quality, min_size)
unique_sequences = sorted(set(high_quality_sequences))

# Save unique sequences to file
with open('sequences.txt', 'w') as output_file:
    output_file.write('\n'.join(unique_sequences))

# Calculate nucleotide composition
def calculate_nucleotide_composition(sequences):
    total_len = sum(len(seq) for seq in sequences)
    composition = {
        'A': round(sum(seq.count('A') for seq in sequences) / total_len * 100, 2),
        'T': round(sum(seq.count('T') for seq in sequences) / total_len * 100, 2),
        'C': round(sum(seq.count('C') for seq in sequences) / total_len * 100, 2),
        'G': round(sum(seq.count('G') for seq in sequences) / total_len * 100, 2),
    }
    return composition

nucleotide_composition = calculate_nucleotide_composition(unique_sequences)

# Print nucleotide composition
print("Nucleotide Composition:")
for nucleotide, percentage in nucleotide_composition.items():
    print(f"{nucleotide}: {percentage}%")

def getReverseComplement(sequence):
    return sequence.upper().replace('T','a').replace('A', 't').replace('G', 'c').replace('C', 'g').upper()[::-1]

def overlap_length(s1, s2):
    max_len = min(len(s1), len(s2))
    for n in range(max_len, 0, -1):
        if s1.endswith(s2[:n]):
            return n
    return 0

def maxOverlap(s1, s2):
    right = overlap_length(s1, s2)
    left = overlap_length(s2, s1)
    return max([right, left])

sequences = [seq.strip() for seq in open("sequences.txt").readlines() if len(seq.strip()) > min_overlap]
for sequence in sequences:
    reverseComplement = getReverseComplement(sequence)
    if reverseComplement in sequences:
        sequences.remove(sequence)
sequences = sorted(sequences)

# Function to generate an Overlap Graph
def overlap_graph(sequences, min_overlap):
    """Constructs an overlap graph from the given sequences."""
    edges = []
    for s1 in sequences:
        for s2 in sequences:
            if s1 != s2:
                overlap = overlap_length(s1, s2)
                if overlap >= min_overlap:
                    edges.append((s1, s2, overlap))
    return edges

# Generate Overlap Graph
sequences = [seq.strip() for seq in open("sequences.txt").readlines() if len(seq.strip()) > min_overlap]

# Create overlap graph edges
overlap_edges = overlap_graph(sequences, min_overlap)

# Build a directed graph using NetworkX
G = nx.DiGraph()
for s1, s2, overlap in overlap_edges:
    G.add_edge(s1, s2, weight=overlap)

# Draw the overlap graph
def draw_overlap_graph(G, output_filename = 'Results/overlap_graph.svg'):
    #Draws the overlap graph and saves it to a file.
    pos = nx.spring_layout(G, seed=42)
    edge_labels = {
        (u, v): f"{G[u][v]['weight']}"
        for u, v in G.edges
    }
    plt.figure(figsize=(20, 20))
    nx.draw(
        G,
        pos,
        with_labels=True,
        node_size=1,
        node_color='lightblue',
        font_size=1,
        font_weight='bold',
        arrowsize=1,
    )
    nx.draw_networkx_edge_labels(G, pos, edge_labels=edge_labels, font_size=8)
    plt.title('Overlap Graph')
    plt.savefig(output_filename)
    plt.clf()

# Function to assemble genome with "-" in overlaps
def assemble_genome_with_gaps(G):
    """Assembles the genome with overlaps represented by a '-' character."""
    assembled_genomes = []
    utilization_scores = []

    # Iterate through connected components of the graph
    for component in nx.weakly_connected_components(G):
        subgraph = G.subgraph(component)
        genome = ""

        # Find a starting node
        start_node = None
        for node in subgraph.nodes:
            if subgraph.in_degree(node) == 0:
                start_node = node
                break
        if start_node is None:
            # Handle cyclic cases
            start_node = next(iter(subgraph.nodes))

        # Traverse the path to build the genome
        current_sequence = start_node
        genome = current_sequence
        visited = set([current_sequence])

        while True:
            # Find the next node with the highest overlap
            successors = list(subgraph.successors(current_sequence))
            if not successors:
                break

            next_sequence = max(
                successors, key=lambda x: G[current_sequence][x]["weight"]
            )
            overlap = G[current_sequence][next_sequence]["weight"]

            genome += "-" + next_sequence[overlap:]
            visited.add(next_sequence)
            current_sequence = next_sequence

        assembled_genomes.append(genome)
        utilization_scores.append(len(visited) / len(component) * 100)

    return assembled_genomes, utilization_scores

def assemble_genome_from_overlap_graph(G):
    """Assembles the genome by traversing the overlap graph."""
    assembled_genomes = []

    # Iterate through connected components of the graph
    for component in nx.weakly_connected_components(G):
        subgraph = G.subgraph(component)
        genome = ""

        # Find a starting node
        start_node = None
        for node in subgraph.nodes:
            if subgraph.in_degree(node) == 0:
                start_node = node
                break
        if start_node is None:
            # Handle cyclic cases
            start_node = next(iter(subgraph.nodes))

        # Traverse the path to build the genome
        current_sequence = start_node
        genome = current_sequence
        visited = set([current_sequence])

        while True:
            # Find the next node with the highest overlap
            successors = list(subgraph.successors(current_sequence))
            if not successors:
                break

            next_sequence = max(
                successors, key=lambda x: G[current_sequence][x]["weight"]
            )
            overlap = G[current_sequence][next_sequence]["weight"]

            genome += next_sequence[overlap:]
            visited.add(next_sequence)
            current_sequence = next_sequence

        assembled_genomes.append(genome)

    return assembled_genomes


# Assemble genomes
assembled_genomes = assemble_genome_from_overlap_graph(G)
assembled_genomes_with_gaps, utilization_scores = assemble_genome_with_gaps(G)

# Save assembled genomes and utilization scores to files
with open('overlap_assembled_genomes.txt', 'w') as output_file:
    output_file.write("\n".join(assembled_genomes))

with open('overlap_assembled_genomes_with_gaps.txt', 'w') as output_file:
    output_file.write("\n".join(assembled_genomes_with_gaps))

with open('utilization_scores.txt', 'w') as output_file:
    output_file.write("\n".join(f"{score:.2f}%" for score in utilization_scores))

# Draw and save the overlap graph
draw_overlap_graph(G)

df = pd.DataFrame([assembled_genomes, assembled_genomes_with_gaps, utilization_scores]).T
df.columns = ['Assembled Genome', 'Assembled Genome with Junctions', 'Raw Sequences Utilization']
df = df.sort_values('Raw Sequences Utilization', ascending = False, ignore_index = True)

# Read the result file
df = df.head()

def getAlignmentSummary(hsp):
    # Parse XML
    root = ET.fromstring(hsp)
    
    # Extract values
    fields = {
        "Hsp_num": root.find("Hsp_num").text,
        "Bit Score": root.find("Hsp_bit-score").text,
        "Score": root.find("Hsp_score").text,
        "E-value": root.find("Hsp_evalue").text,
        "Query Start": root.find("Hsp_query-from").text,
        "Query End": root.find("Hsp_query-to").text,
        "Hit Start": root.find("Hsp_hit-from").text,
        "Hit End": root.find("Hsp_hit-to").text,
        "Query Frame": root.find("Hsp_query-frame").text,
        "Hit Frame": root.find("Hsp_hit-frame").text,
        "Identity": root.find("Hsp_identity").text,
        "Positive Matches": root.find("Hsp_positive").text,
        "Gaps": root.find("Hsp_gaps").text,
        "Alignment Length": root.find("Hsp_align-len").text,
        "Query Sequence": root.find("Hsp_qseq").text,
        "Hit Sequence": root.find("Hsp_hseq").text,
        "Midline": root.find("Hsp_midline").text,
    }
    # Generate HTML
    alignment_summary = f"""
    <!DOCTYPE html>
    <html lang="en">
    <head>
        <meta charset="UTF-8">
        <meta name="viewport" content="width=device-width, initial-scale=1.0">
        <title>BLASTn Result</title>
        <style>
            body {{ font-family: Arial, sans-serif; line-height: 1.6; margin: 20px; }}
            .table {{ border-collapse: collapse; width: 100%; margin: 20px 0; }}
            .table th, .table td {{ border: 1px solid #ddd; padding: 8px; text-align: left; }}
            .table th {{ background-color: #f4f4f4; }}
            .sequences {{ white-space: pre-wrap; font-family: monospace; background: #f9f9f9; padding: 10px; }}
        </style>
    </head>
    <body>
        <h1>BLASTn Result</h1>
        <table class="table">
            <tr><th>Field</th><th>Value</th></tr>
            {''.join(f"<tr><td>{key}</td><td>{value}</td></tr>" for key, value in fields.items() if key not in ['Query Sequence', 'Hit Sequence', 'Midline'])}
        </table>
        <h2>Alignment</h2>
        <div class="sequences">
            <strong>Query Sequence:</strong><br>{fields['Query Sequence']}<br><br>
            <strong>Midline:</strong><br>{fields['Midline']}<br><br>
            <strong>Hit Sequence:</strong><br>{fields['Hit Sequence']}
        </div>
    </body>
    </html>
    """

    alignment_summary = f"""
    <!DOCTYPE html>
    <html lang="en">
    <head>
        <meta charset="UTF-8">
        <meta name="viewport" content="width=device-width, initial-scale=1.0">
        <title>BLASTn Result</title>
        <style>
            body {{ font-family: Arial, sans-serif; line-height: 1.6; margin: 20px; }}
            .table {{ border-collapse: collapse; width: 100%; margin: 20px 0; }}
            .table th, .table td {{ border: 1px solid #ddd; padding: 8px; text-align: left; }}
            .table th {{ background-color: #f4f4f4; }}
            .sequences {{ white-space: pre-wrap; font-family: monospace; background: #f9f9f9; padding: 10px; }}
        </style>
    </head>
    <body>
        <h1>BLASTn Result</h1>
        <table class="table">
            <tr><th>Field</th><th>Value</th></tr>
            {''.join(f"<tr><td>{key}</td><td>{value}</td></tr>" for key, value in fields.items() if key not in ['Query Sequence', 'Hit Sequence', 'Midline'])}
        </table>
        <h2>Alignment</h2>
        <div class="sequences">
            {fields['Query Sequence']}
            {fields['Midline']}
            {fields['Hit Sequence']}
        </div>
    </body>
    </html>
    """
    return alignment_summary


# BLASTn search function
def blastn_search(seq, query_index):
    print(f"Processing sequence {query_index + 1}...")
    try:
        # Submit BLAST query
        result_handle = NCBIWWW.qblast('blastn', 'nt', seq)
        sleep(10)  # Delay to avoid overloading the server

        # Save result to a temporary file
        result = result_handle.read()
        query_id = f"Query_{query_index + 1}"
        output_file = f"{query_id}_blast_output.xml"
        with open(output_file, 'w') as f:
            f.write(result)

        # Parse BLAST output
        xml = open(output_file, 'r').read()
        hsp = '<Hsp>' + xml.split('<Hsp>')[1].split('</Hsp>')[0] + '</Hsp>'
        return getAlignmentSummary(hsp)
    except Exception as e:
        print(f"Error processing sequence {query_index + 1}: {e}")
        return None


# Apply BLASTn search to all sequences
df['BLASTn'] = df['Assembled Genome'].apply(lambda seq: blastn_search(seq, df.index[df['Assembled Genome'] == seq][0]))

# Clean up temporary files
for output_file in glob('*_blast_output.xml'):
    os.remove(output_file)

print("BLASTn processing completed")

# Function to calculate primers using Primer3
def calculate_primers(sequence):
    try:
        primers = primer3.design_primers(
            {
                'SEQUENCE_ID': 'an_id',
                'SEQUENCE_TEMPLATE': sequence,
                'SEQUENCE_TARGET': [350, 100]
            },
            {
                'PRIMER_TASK': 'generic',
                'PRIMER_PICK_LEFT_PRIMER': 1,
                'PRIMER_PICK_INTERNAL_OLIGO': 0,
                'PRIMER_PICK_RIGHT_PRIMER': 1,
                'PRIMER_NUM_RETURN': 10,
                'PRIMER_OPT_SIZE': 20,
                'PRIMER_MIN_SIZE': 18,
                'PRIMER_MAX_SIZE': 27,
                'PRIMER_OPT_TM': 58.0,
                'PRIMER_MIN_TM': 57.0,
                'PRIMER_MAX_TM': 63.0,
                'PRIMER_MIN_GC': 20.0,
                'PRIMER_MAX_GC': 80.0,
                'PRIMER_MAX_POLY_X': 5,
                'PRIMER_SALT_MONOVALENT': 50.0,
                'PRIMER_DNA_CONC': 50.0,
                'PRIMER_MAX_NS_ACCEPTED': 0,
                'PRIMER_MAX_SELF_ANY': 12,
                'PRIMER_MAX_SELF_END': 8,
                'PRIMER_PAIR_MAX_COMPL_ANY': 12,
                'PRIMER_PAIR_MAX_COMPL_END': 8,
                'PRIMER_PRODUCT_SIZE_RANGE': [[100, 800]]
            }
        )
        return primers
    except Exception as e:
        print(f"Primer design failed for sequence. Error: {e}")
        return None

# Function to generate primers for a genome
def generate_primers(genome):
    header = [
        'Primer left', 'Primer right', 
        'Primer left Tm (°C)', 'Primer right Tm (°C)', 
        'Product size (bp)', 'Product sequence'
    ]
    lines = []
    window_size = 800
    step_size = 400

    for i in range(0, len(genome) - window_size, step_size):
        sequence = genome[i:i + window_size]
        primers = calculate_primers(sequence)

        if primers and 'PRIMER_LEFT_1_SEQUENCE' in primers and 'PRIMER_RIGHT_1_SEQUENCE' in primers:
            try:
                # Extract primer information
                primer_left = primers['PRIMER_LEFT_1_SEQUENCE']
                primer_right = primers['PRIMER_RIGHT_1_SEQUENCE']
                primer_left_Tm = round(primers['PRIMER_LEFT_1_TM'], 2)
                primer_right_Tm = round(primers['PRIMER_RIGHT_1_TM'], 2)
                product_size = primers['PRIMER_PAIR_1_PRODUCT_SIZE']

                # Generate reverse complement of right primer
                primer_right_flipped = primer_right[::-1].translate(
                    str.maketrans("ACGT", "TGCA")
                )

                # Extract product sequence
                product_sequence = (
                    sequence.split(primer_left)[1].split(primer_right_flipped)[0]
                )
                product_sequence = primer_left + product_sequence + primer_right_flipped

                # Store results
                line = [
                    primer_left, primer_right, 
                    primer_left_Tm, primer_right_Tm, 
                    product_size, product_sequence
                ]
                lines.append(line)
            except Exception as e:
                print(f"Error processing primers: {e}")
                continue

    return lines

# Apply primer generation to dataframe
df['PCR_validation_settings'] = df['Assembled Genome'].apply(generate_primers)

# Save the results as a pickle file
df.to_pickle('results.pkl')
print("Primer design completed. Results saved to 'results.pkl'.")

# Constants for nucleotide percentages
minSize = 100  # Adjust as needed
As, Ts, Cs, Gs = 25.0, 25.0, 25.0, 25.0  # Example raw percentages

# Function to calculate nucleotide percentages
def nuc_percentage(genome, As, Ts, Cs, Gs, output_filename):
    gAs = round(genome.count('A') / len(genome) * 100, 2)
    gTs = round(genome.count('T') / len(genome) * 100, 2)
    gCs = round(genome.count('C') / len(genome) * 100, 2)
    gGs = round(genome.count('G') / len(genome) * 100, 2)
    labels = ['A', 'T', 'C', 'G']
    sequence_1 = [As, Ts, Cs, Gs]
    sequence_2 = [gAs, gTs, gCs, gGs]

    bar_width = 0.35
    r1 = np.arange(len(sequence_1))
    r2 = [x + bar_width for x in r1]

    plt.bar(r1, sequence_1, width=bar_width, label='Raw Data', alpha=0.7)
    plt.bar(r2, sequence_2, width=bar_width, label='Assembly', alpha=0.7)

    plt.xticks([r + bar_width / 2 for r in range(len(sequence_1))], labels, fontsize = 40)
    plt.yticks(fontsize = 40)
    plt.ylabel('Percentage', fontsize = 40)
    plt.title('Nucleotide Composition', fontsize = 40)
    plt.legend(fontsize=35, markerscale=5)
    plt.tight_layout()
    plt.savefig(output_filename, format='svg')
    plt.clf()

# Function to generate De Bruijn graph
def de_bruijn(sequence, k, output_filename):
    edges = [(sequence[i:i + k - 1], sequence[i + 1:i + k]) for i in range(len(sequence) - k + 1)]
    G = nx.DiGraph(edges)
    pos = nx.spring_layout(G, seed=42)  # Fixed layout for consistent graph
    nx.draw(G, pos, with_labels=True, node_size=10, font_size=1, arrowsize=5)
    plt.title('De Bruijn Graph')
    plt.savefig(output_filename, format='svg')
    plt.clf()

# HTML formatting functions
def P(text):
    return f"<p>{text}</p>"

def H3(text):
    return f"<h3>{text}</h3>"

# Function to generate reports
def generate_report(data):
    global minSize
    genome = data['Assembled Genome']
    genome_with_junction = data['Assembled Genome with Junctions']
    PCRs = data['PCR_validation_settings']
    
    report_id = f"Prediction_{data.name}"
    path = os.path.join("Results", report_id)
    os.makedirs(os.path.join(path, "assets"), exist_ok=True)

    # Save FASTA files
    fasta_content = f">{report_id}\n" + "\n".join(genome[i:i + 70] for i in range(0, len(genome), 70))
    fasta_with_junctions_content = f">{report_id}_with_junctions\n" + "\n".join(
        genome_with_junction[i:i + 70] for i in range(0, len(genome_with_junction), 70)
    )
    with open(os.path.join(path, f"{report_id}.fasta"), "w") as fasta_file:
        fasta_file.write(fasta_content)
    with open(os.path.join(path, f"{report_id}_with_junctions.fasta"), "w") as fasta_junction_file:
        fasta_junction_file.write(fasta_with_junctions_content)

    # Generate De Bruijn graph
    #de_bruijn_path = os.path.join(path, "assets", f"{report_id}_de_bruijn.svg")
    #de_bruijn(genome, minSize, de_bruijn_path)

    # Generate nucleotide percentage chart
    nuc_percentage_path = os.path.join(path, "assets", f"{report_id}_nucleotide_percentages.svg")
    nuc_percentage(genome, As, Ts, Cs, Gs, nuc_percentage_path)

    # Generate HTML report
    html = f"<html><head><title>{report_id}</title></head><body>"
    html += H3("Sequences")
    html += P(fasta_content.replace("\n", "<br>"))
    html += P(fasta_with_junctions_content.replace("\n", "<br>"))
    if data['BLASTn'] != None:
        html += data['BLASTn']
    html += H3("Recommended PCRs for Validation")
    for PCR in PCRs:
        html += P(f"Forward Primer: {PCR[0]}, Tm: {PCR[2]}°C")
        html += P(f"Reverse Primer: {PCR[1]}, Tm: {PCR[3]}°C")
        html += P(f"Product Size: {PCR[4]} bp")
        html += P(f"Product Sequence: <br>{'<br>'.join(PCR[5][i:i + 70] for i in range(0, len(PCR[5]), 70))}")
    html += f'<img src="{os.path.join("assets", f"{report_id}_nucleotide_percentages.svg")}" alt="Nucleotide Percentages" style="width:50%;">'
    html += "</body></html>"

    with open(os.path.join(path, f"{report_id}_summary.html"), "w") as html_file:
        html_file.write(html)

# Main execution
os.makedirs("Results", exist_ok=True)
df.apply(generate_report, axis=1)
print('De-novo assembly report generated.')

