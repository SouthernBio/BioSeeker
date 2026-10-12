"""BioSeeker command-line interface."""

from __future__ import annotations

import argparse
import logging
import sys
import warnings
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path
from typing import List, Optional

from . import __version__
from .analysis.assembler import ConservationAggregator
from .analysis.conservation import conservation_for_records
from .analysis.rscu_analysis import global_rscu, identify_rare_codons
from .analysis.statistics import analyze_pairs, filter_inhibitory_pairs, select_pairs
from .io.reader import discover_alignments, read_fasta
from .visualization.plots import plot_conservation_regression

logger = logging.getLogger("bioseeker")


def _setup_logging(verbose: bool, quiet: bool) -> None:
    level = logging.INFO
    if verbose:
        level = logging.DEBUG
    elif quiet:
        level = logging.WARNING
    logging.basicConfig(level=level, format="%(levelname)s: %(message)s")


def _process_file_conservation(
    file_path: Path,
    orf: int,
    ref: Optional[str],
    thresh: float,
    code: int,
    species: Optional[List[str]] = None,
) -> tuple:
    try:
        records = read_fasta(file_path)
        if not records:
            return None, f"File is empty or unreadable: {file_path}"

        if species:
            record_dict = {r[0]: r[1] for r in records}
            missing = [s for s in species if s not in record_dict]
            if missing:
                return None, f"Skipped {file_path}: missing required species {missing}"
            records = [(s, record_dict[s]) for s in species]

        result = conservation_for_records(
            records, orf, reference_id=ref, threshold=thresh, genetic_code=code
        )
        return result, None
    except Exception as e:
        return None, f"Failed processing {file_path}: {e}"


def cmd_conservation(args: argparse.Namespace) -> int:
    files = discover_alignments(args.input)
    if not files:
        logger.error(f"No alignments found in {args.input}")
        return 1

    logger.info(
        f"Processing {len(files)} alignments for ORFs {args.orfs} (threshold={args.threshold})"
    )

    outdir = Path(args.output)
    outdir.mkdir(parents=True, exist_ok=True)

    has_errors = False
    for orf in args.orfs:
        aggregator = ConservationAggregator()
        with ProcessPoolExecutor(max_workers=args.jobs) as executor:
            futures = [
                executor.submit(
                    _process_file_conservation,
                    f,
                    orf,
                    args.reference,
                    args.threshold,
                    args.genetic_code,
                    getattr(args, "species", None),
                )
                for f in files
            ]
            for f in futures:
                res, err = f.result()
                if err:
                    logger.warning(err)
                    has_errors = True
                elif res:
                    aggregator.add(res)

        if args.strict and has_errors:
            logger.error("Strict mode: aborting due to errors.")
            return 1

        c_df, b_df = aggregator.get_results()
        if not c_df.empty:
            c_out = outdir / f"codon_conservation_ORF{orf}.csv"
            c_df.to_csv(c_out, index=False)
            logger.info(f"Saved {c_out}")
        if not b_df.empty:
            b_out = outdir / f"bicodon_conservation_ORF{orf}.csv"
            b_df.to_csv(b_out, index=False)
            logger.info(f"Saved {b_out}")

    return 0


def cmd_rscu(args: argparse.Namespace) -> int:
    files = discover_alignments(args.input)
    if not files:
        logger.error(f"No alignments found in {args.input}")
        return 1

    outdir = Path(args.output)
    outdir.mkdir(parents=True, exist_ok=True)

    logger.info("Calculating global RSCU...")
    df = global_rscu(files, genetic_code=args.genetic_code)

    rscu_out = outdir / "global_rscu.csv"
    df.to_csv(rscu_out)
    logger.info(f"Saved {rscu_out}")

    rare = identify_rare_codons(df, threshold=args.rare_threshold, genetic_code=args.genetic_code)
    rare_out = outdir / "rare_codons.txt"
    rare_out.write_text("\n".join(rare) + "\n")
    logger.info(f"Saved {rare_out} ({len(rare)} rare codons)")

    return 0


def cmd_analyze(args: argparse.Namespace) -> int:
    outdir = Path(args.output)
    outdir.mkdir(parents=True, exist_ok=True)

    files = discover_alignments(args.input)
    if not files:
        logger.error(f"No alignments found in {args.input}")
        return 1

    rare_codons: List[str] = []
    if args.rare_codons:
        rare_path = Path(args.rare_codons)
        if not rare_path.exists():
            logger.error(f"Rare codons file not found: {rare_path}")
            return 1
        rare_codons = [x.strip().upper() for x in rare_path.read_text().splitlines() if x.strip()]
        logger.info(f"Loaded {len(rare_codons)} rare codons from {rare_path}")

    icp_list: List[str] = []
    if getattr(args, "icp_list", None):
        icp_path = Path(args.icp_list)
        if not icp_path.exists():
            logger.error(f"ICP list file not found: {icp_path}")
            return 1
        icp_list = [x.strip() for x in icp_path.read_text().splitlines() if x.strip()]
        logger.info(f"Loaded {len(icp_list)} known ICPs from {icp_path}")

    for orf in args.orfs:
        logger.info(f"--- ORF {orf} ---")
        aggregator = ConservationAggregator()
        has_errors = False

        with ProcessPoolExecutor(max_workers=args.jobs) as executor:
            futures = [
                executor.submit(
                    _process_file_conservation,
                    f,
                    orf,
                    args.reference,
                    args.threshold,
                    args.genetic_code,
                    getattr(args, "species", None),
                )
                for f in files
            ]
            for f in futures:
                res, err = f.result()
                if err:
                    logger.warning(err)
                    has_errors = True
                elif res:
                    aggregator.add(res)

        if args.strict and has_errors:
            logger.error("Strict mode: aborting due to errors.")
            return 1

        c_df, b_df = aggregator.get_results()
        if c_df.empty or b_df.empty:
            logger.warning(f"No conservation data generated for ORF {orf}.")
            continue

        c_out = outdir / f"codon_conservation_ORF{orf}.csv"
        c_df.to_csv(c_out, index=False)
        b_out = outdir / f"bicodon_conservation_ORF{orf}.csv"
        b_df.to_csv(b_out, index=False)

        stats_df, fit = analyze_pairs(
            c_df,
            b_df,
            aggregator.n_reference_codons,
            genetic_code=args.genetic_code,
            sigma_threshold=args.sigma_threshold,
        )

        stats_out = outdir / f"bicodon_stats_ORF{orf}.csv"
        stats_df.to_csv(stats_out, index=False)
        logger.info(f"ORF {orf} Regression: slope={fit.slope:.4f}, R²={fit.r_squared:.4f}")

        outliers = stats_df[stats_df["Outlier"].isin(["high", "low"])]
        if not outliers.empty:
            outliers_out = outdir / f"outliers_ORF{orf}.csv"
            outliers.to_csv(outliers_out, index=False)
            high_count = (outliers["Outlier"] == "high").sum()
            low_count = (outliers["Outlier"] == "low").sum()
            logger.info(
                f"ORF {orf}: Found {high_count} high outliers (inhibitory pairs) "
                f"and {low_count} low outliers. Saved to {outliers_out.name}"
            )

        if rare_codons:
            inhib = filter_inhibitory_pairs(stats_df, rare_codons)
            if not inhib.empty:
                inhib_out = outdir / f"inhibitory_pairs_ORF{orf}.csv"
                inhib.to_csv(inhib_out, index=False)

        if icp_list:
            known_stats = select_pairs(stats_df, icp_list)
            if not known_stats.empty:
                known_out = outdir / f"known_icp_stats_ORF{orf}.csv"
                known_stats.to_csv(known_out, index=False)

        if not args.no_plot:
            plot_out = outdir / f"conservation_regression_ORF{orf}.png"
            plot_conservation_regression(
                stats_df, fit.slope, fit.r_squared, fit.sigma, str(plot_out)
            )

    return 0


def main(argv: Optional[List[str]] = None) -> int:
    parser = argparse.ArgumentParser(
        description="BioSeeker: DNA Sequence Conservation Analyzer",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument("--version", action="version", version=f"%(prog)s {__version__}")

    # Global options
    group = parser.add_argument_group("Global options")
    group.add_argument("-v", "--verbose", action="store_true", help="Verbose logging")
    group.add_argument("-q", "--quiet", action="store_true", help="Only log warnings and errors")
    group.add_argument("--strict", action="store_true", help="Abort on any file error")
    group.add_argument("-j", "--jobs", type=int, default=None, help="Number of parallel workers")

    # Compatibility arguments for invocation without subcommands
    parser.add_argument("-i", "--input", type=str, help=argparse.SUPPRESS)
    parser.add_argument("-o", "--output", type=str, help=argparse.SUPPRESS)
    parser.add_argument("-r", "--rare-codons", type=str, help=argparse.SUPPRESS)

    subparsers = parser.add_subparsers(title="subcommands", dest="command")

    # Subcommand: conservation
    p_cons = subparsers.add_parser(
        "conservation",
        help="Calculate raw conservation counts",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    p_cons.add_argument("-i", "--input", default="FASTA_files", help="Input file or directory")
    p_cons.add_argument("-o", "--output", default="results", help="Output directory")
    p_cons.add_argument("--reference", help="Reference sequence ID (default: first sequence)")
    p_cons.add_argument(
        "--threshold", type=float, default=1.0, help="Conservation threshold (0, 1]"
    )
    p_cons.add_argument("--genetic-code", type=int, default=1, help="NCBI translation table ID")
    p_cons.add_argument(
        "--orfs", type=int, nargs="+", default=[0, 1, 2], help="Reading frames to analyze"
    )
    p_cons.add_argument(
        "--species",
        nargs="+",
        help="Filter alignments to only include these species",
    )

    # Subcommand: rscu
    p_rscu = subparsers.add_parser(
        "rscu",
        help="Calculate global RSCU and find rare codons",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    p_rscu.add_argument("-i", "--input", default="FASTA_files", help="Input file or directory")
    p_rscu.add_argument("-o", "--output", default="results", help="Output directory")
    p_rscu.add_argument("--genetic-code", type=int, default=1, help="NCBI translation table ID")
    p_rscu.add_argument(
        "--rare-threshold", type=float, default=1.0, help="RSCU threshold for rare codons"
    )

    # Subcommand: analyze
    p_ana = subparsers.add_parser(
        "analyze",
        help="Full analysis pipeline (conservation + stats + plots)",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    p_ana.add_argument("-i", "--input", default="FASTA_files", help="Input file or directory")
    p_ana.add_argument("-o", "--output", default="results", help="Output directory")
    p_ana.add_argument("--reference", help="Reference sequence ID")
    p_ana.add_argument("--threshold", type=float, default=1.0, help="Conservation threshold")
    p_ana.add_argument("--genetic-code", type=int, default=1, help="NCBI translation table ID")
    p_ana.add_argument("--orfs", type=int, nargs="+", default=[0, 1, 2], help="Reading frames")
    p_ana.add_argument(
        "--species",
        nargs="+",
        help="Filter alignments to only include these species",
    )
    p_ana.add_argument(
        "--sigma-threshold", type=float, default=3.0, help="Threshold (in σ) for outliers"
    )
    p_ana.add_argument("--rare-codons", help="File with rare codons to filter inhibitory pairs")
    p_ana.add_argument("--icp-list", help="File with known ICPs to extract their stats")
    p_ana.add_argument("--no-plot", action="store_true", help="Disable plotting")

    args = parser.parse_args(argv)
    _setup_logging(args.verbose, args.quiet)

    # Legacy invocation support
    if args.command is None:
        if args.input:
            warnings.warn(
                "Invoking bioseeker without a subcommand is deprecated. "
                "Use 'bioseeker analyze' instead.",
                DeprecationWarning,
            )
            # Default values for analyze
            args.command = "analyze"
            args.output = args.output or "results"
            args.reference = None
            args.threshold = 1.0
            args.genetic_code = 1
            args.orfs = [0, 1, 2]
            args.species = None
            args.sigma_threshold = 3.0
            args.icp_list = None
            args.no_plot = False
        else:
            parser.print_help()
            return 2

    if args.command == "conservation":
        return cmd_conservation(args)
    elif args.command == "rscu":
        return cmd_rscu(args)
    elif args.command == "analyze":
        return cmd_analyze(args)

    return 0


if __name__ == "__main__":
    sys.exit(main())
