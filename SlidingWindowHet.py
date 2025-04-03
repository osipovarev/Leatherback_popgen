#!/usr/bin/env python3

import sys
import pysam
import os
import gzip


__author__ = "Blair Bentley; edited by chatGPT, 2025"



def get_vcf_samples(vcf_file):
    """Extracts sample names from the VCF file header."""
    with gzip.open(vcf_file, 'rt') as VCF:
        for line in VCF:
            if line.startswith("##"):
                continue
            elif line.startswith("#CHROM"):
                return line.strip().split("\t")[9:]  # Extract sample names
    return []

def get_valid_chromosomes(vcf_file):
    """Extracts chromosomes present in the VCF (with variants)."""
    parsevcf = pysam.Tabixfile(vcf_file)
    return [chrom for chrom in parsevcf.contigs if parsevcf.fetch(chrom)]

def load_chromosome_sizes(chrom_sizes_file):
    """Loads chromosome sizes from a tab-delimited file."""
    with open(chrom_sizes_file, 'r') as file:
        return {line.split()[0]: int(line.split()[1]) for line in file}

def calculate_heterozygosity(parsevcf, chrom, window_start, window_end, samples, output_file):
    """Calculates heterozygosity for a given chromosome window."""
    try:
        rows = tuple(parsevcf.fetch(region=f"{chrom}:{window_start}-{window_end}", parser=pysam.asTuple()))
    except ValueError:
        print(f"Skipping {chrom} (no variants found in this window).")
        return

    sites_total = 0
    calls = [0] * len(samples)
    hets = [0] * len(samples)

    for line in rows:
        if line[6] != "PASS":  # Skip filtered sites
            continue
        sites_total += 1
        for i in range(len(samples)):
            genotype = line[i + 9].split(":")[0]
            if genotype == ".":
                continue  # No call
            calls[i] += 1
            alleles = genotype.replace("|", "/").split("/")
            if len(alleles) == 2 and alleles[0] != alleles[1]:
                hets[i] += 1

    output_file.write(f"{chrom}\t{window_start}\t{sites_total}\t{'\t'.join(map(str, calls))}\t{'\t'.join(map(str, hets))}\n")

def sliding_window_analysis(parsevcf, chrom, window_size, step_size, chrom_sizes, samples, output_file):
    """Performs heterozygosity analysis using a sliding window approach."""
    if chrom not in chrom_sizes:
        print(f"Skipping {chrom}: No chromosome size information available.")
        return
    
    start_pos = 1
    end_pos = chrom_sizes[chrom]
    window_start = start_pos
    window_end = start_pos + window_size - 1

    while window_start < end_pos:
        window_end = min(window_start + window_size - 1, end_pos)
        print(f"Processing: {chrom}:{window_start}-{window_end}")
        calculate_heterozygosity(parsevcf, chrom, window_start, window_end, samples, output_file)
        window_start += step_size

def main():
    """Main function to orchestrate the sliding window heterozygosity analysis."""
    if len(sys.argv) != 5:
        print("Usage: python SlidingWindowHet.py <vcf.gz> <window_size> <step_size> <chrom_sizes.txt>")
        sys.exit(1)

    vcf_file = sys.argv[1]
    window_size = int(sys.argv[2])
    step_size = int(sys.argv[3])
    chrom_sizes_file = sys.argv[4]

    # Load chromosome sizes
    chrom_sizes = load_chromosome_sizes(chrom_sizes_file)

    # Open VCF file and get valid chromosomes
    parsevcf = pysam.Tabixfile(vcf_file)
    valid_chromosomes = get_valid_chromosomes(vcf_file)
    samples = get_vcf_samples(vcf_file)

    if not valid_chromosomes:
        print("No valid chromosomes found in VCF. Exiting.")
        sys.exit(1)

    print(f"Detected {len(valid_chromosomes)} valid chromosomes: {', '.join(valid_chromosomes)}")

    output_filename = f"{vcf_file}_het_{window_size}win_{step_size}step.txt"
    
    with open(output_filename, 'w') as output_file:
        output_file.write("chrom\twindow_start\tsites_total\tcalls\thets\n")

        for chrom in valid_chromosomes:
            print(f"Processing chromosome: {chrom}")
            sliding_window_analysis(parsevcf, chrom, window_size, step_size, chrom_sizes, samples, output_file)

    print(f"Analysis complete. Results saved to {output_filename}")

if __name__ == "__main__":
    main()
