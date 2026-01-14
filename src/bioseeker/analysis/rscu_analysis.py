import csv
import argparse
from pathlib import Path
from Bio import SeqIO
from Bio.Data import CodonTable
from collections import defaultdict, Counter


def get_standard_genetic_code():
    """
    Returns a dictionary mapping Amino Acids to lists of synonymous Codons.
    Excludes Stop codons.
    """
    table = CodonTable.unambiguous_dna_by_name["Standard"]
    aa_to_codons = defaultdict(list)
    for codon, aa in table.forward_table.items():
        aa_to_codons[aa].append(codon)
    
    return aa_to_codons


def count_codons_for_sequence(sequence):
    """
    Counts codons for a given DNA sequence string.
    Returns a Counter {codon: count}.
    """
    # 1. Clean sequence
    # Remove gaps
    seq_clean = sequence.replace("-", "").upper()
    
    # 2. Count codons
    codon_counts = Counter()
    # We only consider triplets
    for i in range(0, len(seq_clean) - len(seq_clean) % 3, 3):
        codon = seq_clean[i:i+3]
        if "N" in codon: # Skip ambiguous
            continue
        codon_counts[codon] += 1
    
    return codon_counts


def calculate_rscu_from_counts(codon_counts):
    """
    Calculates RSCU from a dictionary/Counter of codon counts.
    Returns a dictionary {codon: rscu_value}.
    """
    aa_to_codons = get_standard_genetic_code()
    rscu_values = {}
    
    # Get all valid codons checks
    
    for aa, synonymous_codons in aa_to_codons.items():
        # n_i: number of synonymous codons for this amino acid
        n_i = len(synonymous_codons)
        
        # Sum of occurrences of all synonymous codons for this amino acid
        total_aa_count = sum(codon_counts[c] for c in synonymous_codons)
        
        if total_aa_count == 0:
            for codon in synonymous_codons:
                rscu_values[codon] = 0.0
            continue
            
        # expected frequency if user equally = total_aa_count / n_i
        expected_count = total_aa_count / n_i
        
        for codon in synonymous_codons:
            observed = codon_counts[codon]
            rscu = observed / expected_count
            rscu_values[codon] = rscu
            
    return rscu_values


def calculate_rscu_for_sequence(sequence):
    """
    Wrapper for backward compatibility or simple use cases.
    """
    counts = count_codons_for_sequence(sequence)
    return calculate_rscu_from_counts(counts)


def process_msa_file(file_path):
    """
    Processes a single MSA file and returns a dict {Species: Counter}.
    """
    file_path = Path(file_path)
    if not file_path.exists():
        print(f"Error: File {file_path} does not exist.")
        return {}

    species_counts = defaultdict(Counter)
    
    try:
        # Assuming format is FASTA
        msa = list(SeqIO.parse(str(file_path), "fasta"))
    except Exception as e:
        print(f"Error parsing {file_path}: {e}")
        return {}

    for record in msa:
        species_name = record.id # e.g., "Dmel", "Dsim"
        sequence = str(record.seq)
        counts = count_codons_for_sequence(sequence)
        species_counts[species_name].update(counts)
            
    return species_counts


def identify_rare_codons(rscu_data, output_path, threshold=1.0):
    """
    Identifies codons where RSCU < threshold for ALL species.
    Writes to output_path.
    rscu_data structure: {Codon: {Species: RSCU_Value}} (values are floats or strings)
    """
    rare_codons = []
    
    # Sort codons for consistent output
    sorted_codons = sorted(rscu_data.keys())
    
    for codon in sorted_codons:
        species_values = rscu_data[codon]
        # Check if all species have RSCU < threshold
        # Values in data might be strings "0.0000". Convert to float.
        if all(float(val) < threshold for val in species_values.values()):
            rare_codons.append(codon)
            
    try:
        with open(output_path, 'w') as f:
            for codon in rare_codons:
                f.write(f"{codon}\n")
        print(f"Generated rare codons file at {output_path} (Threshold < {threshold})")
    except IOError as e:
        print(f"Error writing rare codons: {e}")


def process_directory(input_path, output_csv):
    """
    Iterates over all MSA files, aggregates counts, and writes one CSV.
    output_csv is a FILE path, not directory.
    """
    input_path = Path(input_path)
    output_csv = Path(output_csv)
    
    if not output_csv.parent.exists():
        output_csv.parent.mkdir(parents=True, exist_ok=True)
        
    global_counts = defaultdict(Counter) # {Species: Counter}
    
    # ... (code to read files)
    if input_path.is_file():
        file_counts = process_msa_file(input_path)
        for sp, counts in file_counts.items():
            global_counts[sp].update(counts)
    elif input_path.is_dir():
        files = list(input_path.glob('*.afa'))
        print(f"Found {len(files)} .afa files in {input_path}")
        
        for i, f in enumerate(files):
            if i % 100 == 0:
                print(f"Processing file {i+1}/{len(files)}: {f.name}", end='\r')
            file_counts = process_msa_file(f)
            for sp, counts in file_counts.items():
                global_counts[sp].update(counts)
        print(f"\nFinished processing {len(files)} files.")
    else:
        print("Invalid input path.")
        return

    # Now calculate RSCU for global counts
    code_map = get_standard_genetic_code()
    all_codons = sorted([c for codons in code_map.values() for c in codons])
    
    data = {codon: {} for codon in all_codons}
    species_list = sorted(global_counts.keys())
    
    print("Calculating RSCU...")
    
    # Keep track of numerical values for rare identification
    # Start with explicit loop
    
    for sp in species_list:
        counts = global_counts[sp]
        rscu_map = calculate_rscu_from_counts(counts)
        for codon in all_codons:
            val = rscu_map.get(codon, 0.0)
            data[codon][sp] = val # Store float first
            
    # Write CSV
    try:
        with open(output_csv, 'w', newline='') as csvfile:
            fieldnames = ['Codon'] + species_list
            writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
            
            writer.writeheader()
            for codon in all_codons:
                row = {'Codon': codon}
                for sp in species_list:
                    # Write formatted
                    row[sp] = f"{data[codon][sp]:.4f}"
                writer.writerow(row)
        print(f"Generated {output_csv}") 
    except IOError as e:
        print(f"Error writing to {output_csv}: {e}")
        
    # Generate Rare Codons
    rare_codons_path = output_csv.parent / "rare_codons.txt"
    identify_rare_codons(data, rare_codons_path)


def main():
    parser = argparse.ArgumentParser(description="Calculate Global RSCU from MSA FASTA files.")
    parser.add_argument('input_path', help="Input file or directory containing .afa files")
    parser.add_argument('output_csv', help="Output CSV file path")
    
    args = parser.parse_args()
    
    process_directory(args.input_path, args.output_csv)


if __name__ == "__main__":
    main()
