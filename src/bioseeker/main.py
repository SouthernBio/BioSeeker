import argparse
import os
import glob
import logging
from .io.reader import read_fasta
from .analysis.conservation import calculate_conservation
from .analysis.assembler import ConservationAggregator
from .analysis.statistics import (
    calculate_expected_bicodon_rates,
    perform_regression,
    calculate_z_scores,
    filter_inhibitory_pairs
)
from .visualization.plots import plot_conservation_regression


def setup_logging():
    logging.basicConfig(
        level=logging.INFO,
        format='%(asctime)s - %(levelname)s - %(message)s'
    )


def main():
    setup_logging()
    parser = argparse.ArgumentParser(description="BioSeeker: DNA Sequence Conservation Analyzer")
    parser.add_argument("--input", "-i", type=str, default="FASTA_files", help="Input directory containing FASTA/AFA files")
    parser.add_argument("--output", "-o", type=str, default="results", help="Output directory for results")
    parser.add_argument("--rare-codons", "-r", type=str, help="Path to text file containing rare codons (one per line)")
    args = parser.parse_args()

    input_dir = args.input
    output_dir = args.output
    rare_codons_path = args.rare_codons

    if not os.path.exists(input_dir):
        logging.error(f"Input directory '{input_dir}' does not exist.")
        return

    rare_codons = []
    if rare_codons_path:
        if os.path.exists(rare_codons_path):
            with open(rare_codons_path, 'r') as f:
                rare_codons = [line.strip().upper() for line in f if line.strip()]
            logging.info(f"Loaded {len(rare_codons)} rare codons.")
        else:
            logging.error(f"Rare codons file '{rare_codons_path}' not found.")
            return

    # Find files (supporting .afa, .fasta, .fa)
    files = glob.glob(os.path.join(input_dir, "*.afa")) + \
            glob.glob(os.path.join(input_dir, "*.fasta")) + \
            glob.glob(os.path.join(input_dir, "*.fa"))

    if not files:
        logging.error(f"No FASTA files found in {input_dir}")
        return

    logging.info(f"Found {len(files)} alignment files.")

    # Initialize aggregators for each ORF
    aggregators = {
        0: ConservationAggregator(),
        1: ConservationAggregator(),
        2: ConservationAggregator()
    }

    # Process files
    for file_path in files:
        logging.info(f"Processing {os.path.basename(file_path)}...")
        try:
            sequences_data = read_fasta(file_path)
            raw_sequences = [s[1] for s in sequences_data]
            
            if not raw_sequences:
                logging.warning(f"File {file_path} is empty or unreadable.")
                continue

            for orf in [0, 1, 2]:
                try:
                    c_df, b_df = calculate_conservation(raw_sequences, orf)
                    aggregators[orf].add(c_df, b_df)
                except ValueError as ve:
                    logging.warning(f"Skipping ORF {orf} for {file_path}: {ve}")
                    
        except Exception as e:
            logging.error(f"Failed to process {file_path}: {e}")

    # Save results
    os.makedirs(output_dir, exist_ok=True)
    
    for orf in [0, 1, 2]:
        c_final, b_final = aggregators[orf].get_results()
        
        if not c_final.empty:
            c_out = os.path.join(output_dir, f"codon_conservation_ORF{orf}.csv")
            c_final.to_csv(c_out, index=False)
            logging.info(f"Saved {c_out}")
            
        if not b_final.empty:
            b_out = os.path.join(output_dir, f"bicodon_conservation_ORF{orf}.csv")
            b_final.to_csv(b_out, index=False)
            logging.info(f"Saved {b_out}")
            
            # --- New Statistical Features ---
            # Calculate Expected Rates
            stats_df = calculate_expected_bicodon_rates(c_final, b_final)
            
            # Z-Score Analysis
            stats_df = calculate_z_scores(stats_df)
            
            # Inhibitory Pairs Filter
            if rare_codons:
                inhibitory_df = filter_inhibitory_pairs(stats_df, rare_codons)
                inhib_out = os.path.join(output_dir, f"inhibitory_pairs_ORF{orf}.csv")
                inhibitory_df.to_csv(inhib_out, index=False)
                logging.info(f"Saved {inhib_out}")
            
            # Save full stats
            stats_out = os.path.join(output_dir, f"bicodon_stats_ORF{orf}.csv")
            stats_df.to_csv(stats_out, index=False)
            logging.info(f"Saved {stats_out}")

            # Plotting and Regression
            if orf == 0:
                slope, r2 = perform_regression(stats_df)
                logging.info(f"ORF{orf} Regression: slope={slope:.4f}, R2={r2:.4f}")
                
                plot_out = os.path.join(output_dir, f"conservation_regression_ORF{orf}.png")
                plot_conservation_regression(stats_df, slope, r2, plot_out)
                logging.info(f"Saved {plot_out}")

    logging.info("Analysis complete.")


if __name__ == "__main__":
    main()
