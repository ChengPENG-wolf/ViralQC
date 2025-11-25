from Bio import SeqIO
import pickle as pkl
import pandas as pd
import subprocess
import shutil
import os


def summarize_sequence_info(fasta_file, database_pth, sequence_info_file):
    reference_genome_length = {}
    for record in SeqIO.parse(fasta_file, 'fasta'):
        reference_genome_length[record.id] = len(record.seq)
    pkl.dump(reference_genome_length, open(f'{database_pth}/reference_genome_length.pkl', 'wb'))

    if os.path.exists(sequence_info_file):
        reference_genome_tax = {}
        df = pd.read_csv(sequence_info_file, sep='\t')
        for i, row in df.iterrows():
            reference_genome_tax[row['seq_id']] = row['Taxonomic Lineage']
        pkl.dump(reference_genome_tax, open(f'{database_pth}/reference_genome_tax.pkl', 'wb'))


def predict_protein(fasta_file, database_pth, thread):
    temp_pth = os.path.join(database_pth, 'temp')

    cmd = f"python parallel-prodigal-gv.py -t {thread} -q -i {fasta_file} -a {temp_pth}/reference_protein.faa"
    _ = subprocess.check_call(cmd, shell=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)

    cmd = f"diamond makedb --in {database_pth}/reference_protein.faa -d {database_pth}/reference_protein -p {thread}"
    _ = subprocess.check_call(cmd, shell=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)


def cluster(database_pth, thread):
    temp_pth = os.path.join(database_pth, 'temp')

    cmd = f'mmseqs easy-cluster {database_pth}/reference_protein.faa {temp_pth}/reference_protein {temp_pth} --threads {thread} --min-seq-id 0.5 -c 0.8'
    _ = subprocess.check_call(cmd, shell=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)


def construct_protein_cluster(database_pth):
    temp_pth = os.path.join(database_pth, 'temp')

    cluster = {}
    cluster_member = {}
    for line in open(os.path.join(temp_pth, 'reference_protein_cluster.tsv'), 'r'):
        line = line.strip().split()
        cluster[line[1]] = line[0]
        if line[0] not in cluster_member:
            cluster_member[line[0]] = set()
        cluster_member[line[0]].add(line[1])

    # sort protein by cluster size
    cluster_size = {cluster_id: len(cluster_member[cluster_id]) for cluster_id in cluster_member}
    cluster_size = {k: v for k, v in sorted(cluster_size.items(), key=lambda item: item[1], reverse=True)}
    cluster2id = {cluster_id: i for i, cluster_id in enumerate(cluster_size)}
    pkl.dump(cluster2id, open(os.path.join(database_pth, 'reference_protein_cluster2id.pkl'), 'wb'))

    for protein in cluster:
        cluster[protein] = cluster2id[cluster[protein]]
    pkl.dump(cluster, open(os.path.join(database_pth, 'reference_protein_cluster.pkl'), 'wb'))

    genome_info = {}
    protein_info = {}
    f = open(os.path.join(database_pth, 'reference_protein_info.tsv'), 'w')
    f.write('Genome Accession\tProtein ID\tStart\tEnd\tStrand\tCluster\n')
    for record in SeqIO.parse(f'{database_pth}/reference_protein.faa', 'fasta'):
        protein_id = record.id
        genome_id = protein_id.rsplit('_', 1)[0]
        representative = cluster[protein_id]
        description = record.description.split(' # ')
        start = int(description[1])
        end = int(description[2])
        strand = int(description[3])
        if genome_id not in genome_info:
            genome_info[genome_id] = []
            protein_info[genome_id] = {}

        genome_info[genome_id].append([protein_id, representative])
        protein_info[genome_id][protein_id] = {'start': start, 'end': end, 'strand': strand, 'cluster': representative}
        f.write(f'{genome_id}\t{protein_id}\t{start}\t{end}\t{strand}\t{representative}\n')

    pkl.dump(genome_info, open(f'{database_pth}/temp/reference_genome_info.pkl', 'wb'))
    f.close()


def build_custom_database(fasta_file, database_pth, sequence_info_file=None, thread=1):
    genome_database_pth = os.path.join(database_pth, 'genome_database')
    temp_pth = os.path.join(database_pth, 'temp')

    if not os.path.exists(database_pth):
        os.makedirs(database_pth)
    if not os.path.exists(genome_database_pth):
        os.makedirs(genome_database_pth)
    if not os.path.exists(temp_pth):
        os.makedirs(temp_pth)

    summarize_sequence_info(fasta_file, genome_database_pth, sequence_info_file)
    predict_protein(fasta_file, database_pth, thread)
    cluster(database_pth, thread)
    construct_protein_cluster(database_pth)

    if os.path.exists(temp_pth):
        shutil.rmtree(temp_pth)


