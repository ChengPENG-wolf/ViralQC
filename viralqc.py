import os.path

from utils.download_database import download_database
from utils.contamination import contamination
from utils.completeness import completeness
import argparse


def run_download_database(args):
    db_path = args.db
    download_database(db_path)


def run_contamination(args):
    input = args.input
    db = args.db
    output = args.output
    threads = args.threads
    contamination(input, db, output, threads)

def run_completeness(args):
    input = args.input
    db = os.path.join(args.db, 'genome_database')
    output = args.output
    threads = args.threads
    bin = args.bin
    completeness(input, db, output, threads, bin)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(prog='ViralQC', description='ViralQC is a python library for quality assessment of '
                                                 'assembled viral contigs or bins. ViralQC contains two '
                                                 'primary modules. The contamination detection module '
                                                 'identifies and removes non-viral regions within viral contigs. '
                                                 'The completeness module evaluates the expected genome length and '
                                                 'calculates the completeness of the viral contigs.')
    subparsers = parser.add_subparsers()

    download_parser = subparsers.add_parser('download_database', help='Download the database of ViralQC.')
    download_parser.add_argument('--db', type=str, help='The path to save the database.')
    download_parser.set_defaults(func=run_download_database)

    contamination_parser = subparsers.add_parser('contamination', help='Detect and remove contamination in viral contigs.')
    contamination_parser.add_argument('--input', type=str, help='The path to the input fasta file.')
    contamination_parser.add_argument('--db', type=str, help='The path to the database.')
    contamination_parser.add_argument('--output', type=str, help='Output directory.')
    contamination_parser.add_argument('--threads', type=int, help='Number of threads.')
    contamination_parser.set_defaults(func=run_contamination)

    completeness_parser = subparsers.add_parser('completeness', help='Evaluate the completeness of viral contigs.')
    completeness_parser.add_argument('--input', type=str, help='The path to the input fasta file.')
    completeness_parser.add_argument('--db', type=str, help='The path to the database.')
    completeness_parser.add_argument('--output', type=str, help='Output directory.')
    completeness_parser.add_argument('--threads', type=int, help='Number of threads.')
    completeness_parser.add_argument('--bin', action='store_true', help='Evaluate the completeness of bins.')
    completeness_parser.set_defaults(func=run_completeness)



    #end_to_end_parser = subparsers.add_parser('end_to_end', help='Run the end-to-end pipeline of ViralQC.')

    args = parser.parse_args()
    args.func(args)
