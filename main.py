"""
Main Orchestrator for the Aether - Production Ready.
CLI & Orchestration Engine for High-Throughput QC & Variant Caller.
"""

import os
import sys
import argparse

from rich.console import Console
from rich.panel import Panel
from rich.progress import Progress, SpinnerColumn, TextColumn, TimeElapsedColumn
from rich.table import Table

from src.io_handler import stream_genomic_records
from src.filter_engine import is_high_quality, calculate_aether_purity_score
from src.variant_caller import scan_for_clinical_variants
from src.reporter import generate_clinical_report, compute_execution_summary

def main():
    parser = argparse.ArgumentParser(description="BioProcessor Core - High-Throughput QC & Variant Caller")
    parser.add_argument("-i", "--input", required=True, help="Path to raw or compressed FASTQ file (.fastq, .fastq.gz)")
    parser.add_argument("-o", "--output", default=os.path.join("data", "processed", "clinical_report.csv"), help="Path for final CSV report")
    parser.add_argument("-d", "--db", default=os.path.join("data", "metadata", "clinical_variants.csv"), help="Path to clinical variants CSV database")
    parser.add_argument("--min-phred", type=float, default=30.0, help="Minimum global average Phred quality threshold")
    parser.add_argument("--min-length", type=int, default=35, help="Minimum sequence length after 3' trimming")
    
    args = parser.parse_args()

    console = Console()

    # Branded Header
    console.print(Panel("[bold cyan]AETHER CORE - CASTER BIOTECH[/bold cyan]", expand=False))

    # Exception Handling for Input Files
    if not os.path.exists(args.input):
        console.print(Panel(f"[bold red]Error:[/bold red] Input FASTQ file not found at '{args.input}'", title="File Not Found", border_style="red"))
        sys.exit(1)
        
    if not os.path.exists(args.db):
        console.print(Panel(f"[bold red]Error:[/bold red] Clinical variants database not found at '{args.db}'", title="Database Not Found", border_style="red"))
        sys.exit(1)

    # Ensure output directory exists
    os.makedirs(os.path.dirname(args.output), exist_ok=True)

    detected_variants = []
    purity_scores = []
    total_reads = 0
    passed_reads = 0
    dropped_reads = 0

    try:
        raw_reads = stream_genomic_records(args.input)
        
        with Progress(
            SpinnerColumn(),
            TextColumn("[progress.description]{task.description}"),
            TextColumn("[progress.completed]{task.completed} reads processed"),
            TimeElapsedColumn(),
            console=console
        ) as progress:
            
            task = progress.add_task("[cyan]Streaming and processing FASTQ data...", total=None)
            
            for read in raw_reads:
                total_reads += 1
                
                # Filter Engine
                if is_high_quality(read, min_phred_avg=args.min_phred, min_length=args.min_length):
                    passed_reads += 1
                    
                    # Compute Purity Score
                    purity = calculate_aether_purity_score(read)
                    purity_scores.append(purity)
                    
                    # Variant Caller
                    variant_report = scan_for_clinical_variants(read, db_path=args.db)
                    
                    if variant_report:
                        # Append the purity score to the report dict for the reporter
                        variant_report["Aether_Purity_Score"] = purity
                        detected_variants.append(variant_report)
                else:
                    dropped_reads += 1
                    
                progress.advance(task)
                
    except Exception as e:
        console.print(Panel(f"[bold red]Execution Error:[/bold red] {str(e)}", title="Pipeline Failure", border_style="red"))
        sys.exit(1)

    # Generate Report
    report_generated = generate_clinical_report(detected_variants, args.output)
    
    # Compute Summary
    summary = compute_execution_summary(
        total_reads=total_reads,
        passed_reads=passed_reads,
        dropped_reads=dropped_reads,
        detected_variants_list=detected_variants,
        purity_scores_list=purity_scores
    )

    # Final Stats Table
    table = Table(title="Execution Summary", show_header=True, header_style="bold magenta")
    table.add_column("Metric", style="cyan")
    table.add_column("Value", style="green")

    table.add_row("Total Reads Processed", f"{summary['Total_Reads']:,}")
    table.add_row("Reads Passed QC", f"{summary['Passed_Reads']:,} ({summary['Pass_Rate_Percent']:.2f}%)")
    table.add_row("Reads Dropped", f"{summary['Dropped_Reads']:,} ({summary['Drop_Rate_Percent']:.2f}%)")
    table.add_row("Variants Identified", f"{summary['Variants_Identified']:,}")
    table.add_row("Mean Aether Purity Score", f"{summary['Mean_Aether_Purity_Score']:.2f}")

    console.print(table)

    if report_generated:
        console.print(f"\n[bold green]Success:[/bold green] Clinical report saved at {args.output}")
    else:
        console.print("\n[bold yellow]Info:[/bold yellow] Pipeline completed. No clinical variants detected.")

if __name__ == "__main__":
    main()