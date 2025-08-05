#!/usr/bin/env python3

from Bio import SeqIO
import os
import argparse

def convert_fasta(fasta_path):
    # Check if file exists
    if not os.path.isfile(fasta_path):
        print(f"File not found: {fasta_path}")
        return 1

    out_dir = os.path.dirname(fasta_path)
    out_path = os.path.join(out_dir, "Trinity_1line.fasta")

    print(f"Processing: {fasta_path} -> {out_path}")

    with open(fasta_path, "r") as handle, open(out_path, "w") as output:
        record_count = 0
        for record in SeqIO.parse(handle, "fasta"):
            output.write(f">{record.id}\n{str(record.seq)}\n")
            record_count += 1

        print(f"Wrote {record_count} records to {out_path}")

    return 0

def main():
    parser = argparse.ArgumentParser(description="Convert a multi-line FASTA file to one-line-per-sequence format.")
    parser.add_argument("fasta_file", help="Path to the input FASTA file.")

    args = parser.parse_args()
    convert_fasta(args.fasta_file)

if __name__ == "__main__":
    main()
