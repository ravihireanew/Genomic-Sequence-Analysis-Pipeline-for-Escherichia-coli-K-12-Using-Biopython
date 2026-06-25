# Genomic Sequence Analysis Pipeline for E. coli K-12

A beginner-to-intermediate bioinformatics project using Biopython for genome analysis.

## Overview
This project builds a complete pipeline to:
- Fetch bacterial genome data from NCBI
- Analyze GC content and generate visualizations
- Extract genes and protein sequences
- Perform BLAST homology search
- Analyze codon usage

## Features
- GC content calculation + sliding window plot
- Gene/CDS feature parsing (4,318 CDS found)
- Protein sequence extraction and FASTA export
- Online BLASTp search
- Codon usage statistics
- Summary report generation

## Key Results
- Genome Size: 4,641,652 bp
- GC Content: 50.79%
- Top BLAST hit for thrL: 100% identity to *E. coli*
- Most used codon in thrA: CTG (Leucine)

## Files
- `py.py` – Main analysis script
- `e_coli_k12.gb` – Downloaded genome
- `e_coli_proteins.fasta` – Extracted proteins
- `e_coli_gc_plot.png` – GC visualization

## How to Run
1. Install dependencies:
   ```bash
   pip install biopython matplotlib
