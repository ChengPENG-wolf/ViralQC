import os
import argparse


def run_download_database(args):
    from utils.download_database import download_database
    db_path = args.db
    download_database(db_path)


def run_build_custom_database(args):
    from utils.custom_database import build_custom_database
    fasta_file = args.fasta
    database_pth = args.db
    sequence_info_file = args.info
    thread = args.threads
    build_custom_database(fasta_file, database_pth, sequence_info_file, thread)


def run_contamination(args):
    from utils.contamination import contamination
    input = args.input
    db = args.db
    output = args.output
    threads = args.threads
    contamination(input, db, output, threads)


def run_completeness(args):
    from utils.completeness import completeness
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

    custom_parser = subparsers.add_parser('build_custom_database', help='Build a custom database for ViralQC.')
    custom_parser.add_argument('--fasta', '-f', type=str, help='The path to the custom viral genome fasta file.')
    custom_parser.add_argument('--db', '-d', type=str, help='The path to save the custom database.')
    custom_parser.add_argument('--info', '-s', type=str, default=None, help='The path to the sequence info file.')
    custom_parser.add_argument('--threads', '-t', type=int, help='Number of threads.')
    custom_parser.set_defaults(func=run_build_custom_database)

    contamination_parser = subparsers.add_parser('contamination', help='Detect and remove contamination in viral contigs.')
    contamination_parser.add_argument('--input', '-i', type=str, help='The path to the input fasta file.')
    contamination_parser.add_argument('--db', '-d', type=str, help='The path to the database.')
    contamination_parser.add_argument('--output', '-o', type=str, help='Output directory.')
    contamination_parser.add_argument('--threads', '-t', type=int, help='Number of threads.')
    contamination_parser.set_defaults(func=run_contamination)

    completeness_parser = subparsers.add_parser('completeness', help='Evaluate the completeness of viral contigs.')
    completeness_parser.add_argument('--input', '-i', type=str, help='The path to the input fasta file.')
    completeness_parser.add_argument('--db', '-d', type=str, help='The path to the database.')
    completeness_parser.add_argument('--output', '-o', type=str, help='Output directory.')
    completeness_parser.add_argument('--threads', '-t', type=int, help='Number of threads.')
    completeness_parser.add_argument('--bin', '-b', action='store_true', help='Evaluate the completeness of bins.')
    completeness_parser.set_defaults(func=run_completeness)

    args = parser.parse_args()
    args.func(args)
