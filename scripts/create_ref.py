################################################################################
# Create reference genome for doing simulations for SCiMS
#
# Author: Kobie Kirven
################################################################################

# Imports
import argparse
import os
import sys
import subprocess
from Bio import SeqIO
import gzip

# Get command line arguments
parser = argparse.ArgumentParser(description='Create reference genome for doing simulations for SCiMS')
parser.add_argument('-r', '--reference', help='Reference genome', required=True)
parser.add_argument('-o', '--output', help='Output file', required=True)

scaffolds = ["NC_000001.11", "NC_000002.12", "NC_000003.12", "NC_000004.12", "NC_000005.10", "NC_000006.12", "NC_000007.14",
             "NC_000008.11", "NC_000009.12", "NC_000010.11", "NC_000011.10", "NC_000012.12", "NC_000013.11", "NC_000014.9",
             "NC_000015.10", "NC_000016.10", "NC_000017.11", "NC_000018.10", "NC_000019.10", "NC_000020.11", "NC_000021.9",
             "NC_000022.11", "NC_000023.11", "NC_000024.10"]
#[x.strip("\n") for x in open("scaffolds.txt", 'r')]
args = parser.parse_args()


# Create female reference genome
# Put two copies of the reference genome in the output file except for 
# the y chromosome
with open(f"{args.output}_female.fasta", 'w') as out:
    with gzip.open(args.reference, 'rt') as fasta_file:
        for record in SeqIO.parse(fasta_file, 'fasta'):
            if record.id != 'NC_000024.10' and record.id in scaffolds:
                out.write(f">{record.id}_1\n{record.seq}\n")
                out.write(f">{record.id}_2\n{record.seq}\n")

# Create the male reference genome
# Put two copies of the reference genome in the output file except for
# the x chromosome and the y chromosome
with open(f"{args.output}_male.fasta", 'w') as out:
    with gzip.open(args.reference, 'rt') as fasta_file:
        for record in SeqIO.parse(fasta_file, 'fasta'):
            if record.id != 'NC_000023.11' and record.id != 'NC_000024.10' and record.id in scaffolds:
                out.write(f">{record.id}_1\n{record.seq}\n")
                out.write(f">{record.id}_2\n{record.seq}\n")
            elif record.id == 'NC_000023.11':
                out.write(f">{record.id}_1\n{record.seq}\n")
            elif record.id == 'NC_000024.10':
                out.write(f">{record.id}_2\n{record.seq}\n")


