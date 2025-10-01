import numpy as np
from Bio import SeqIO
import pickle as pkl
import subprocess
from utils.genome import GenomeStructure
import scipy.stats as stats
import resource
import time
import os


def get_query_length_dtr(fasta_pth, out_pth, min_length=21):
    temp_pth = os.path.join(out_pth, 'midfolder')
    query_length = {}
    query_dtr = {}
    for record in SeqIO.parse(fasta_pth, 'fasta'):
        query_length[record.id] = len(record.seq)

        seq = str(record.seq).casefold()
        substring = seq.casefold()[0:min_length]
        pos = seq.casefold().rfind(substring)
        if pos < len(seq) / 2:
            query_dtr[record.id] = False
        else:
            substring = seq.casefold()[pos:]
            query_dtr[record.id] = (seq.casefold()[:len(substring)] == substring)
    pkl.dump(query_length, open(os.path.join(temp_pth, 'query_seq_length.pkl'), 'wb'))
    pkl.dump(query_dtr, open(os.path.join(temp_pth, 'query_dtr.pkl'), 'wb'))


def get_query_length_bin(fasta_pth, out_pth):
    temp_pth = os.path.join(out_pth, 'midfolder')
    query_length = {'bin': 0}
    for record in SeqIO.parse(fasta_pth, 'fasta'):
        query_length['bin'] += len(record.seq)
    pkl.dump(query_length, open(f'{temp_pth}/query_seq_length.pkl', 'wb'))


def protein_prediction_and_alignment(database_pth, fasta_pth, out_pth, thread):
    print("[1/4] Calling genes with prodigal...")
    temp_pth = os.path.join(out_pth, 'midfolder')
    if not os.path.exists(os.path.join(temp_pth, 'query_protein.faa')):
        cmd = f"python utils/parallel-prodigal-gv.py -t {thread} -q -i {fasta_pth} -a {temp_pth}/query_protein.faa"
        _ = subprocess.check_call(cmd, shell=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    print("[2/4] Running all against all alignment...")
    # diamond
    cmd = f"diamond blastp -d {database_pth}/reference_protein.dmnd -q {temp_pth}/query_protein.faa -p {thread} -o {temp_pth}/query_protein.blastp --query-cover 50 --subject-cover 50 -k 1000 --ultra-sensitive"
    _ = subprocess.check_call(cmd, shell=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)


def parse_alignment_result(database_pth, out_pth):
    temp_pth = os.path.join(out_pth, 'midfolder')
    reference_protein_cluster = pkl.load(open(f'{database_pth}/reference_protein_cluster.pkl', 'rb'))

    f = open(f'{temp_pth}/query_protein_blastp.tsv', 'w')
    f.write('Contig_id\tQuery\tSubject\tCluster\tIdentity\tE-value\tBitscore\n')

    # parse diamond result
    query_protein_alignment = {}
    for line in open(f'{temp_pth}/query_protein.blastp', 'r'):
        line = line.strip().split('\t')
        query = line[0]
        query_contig_id = query.rsplit('_', 1)[0]
        subject = line[1]
        cluster = reference_protein_cluster[subject]
        identity = float(line[2])
        evalue = float(line[10])
        bitscore = float(line[11])

        if query_contig_id not in query_protein_alignment:
            query_protein_alignment[query_contig_id] = {}
        if query not in query_protein_alignment[query_contig_id]:
            query_protein_alignment[query_contig_id][query] = []
        query_protein_alignment[query_contig_id][query].append([subject, cluster, identity, evalue, bitscore])
        f.write(f'{query_contig_id}\t{query}\t{subject}\t{reference_protein_cluster[subject]}\t{identity}\t{evalue}\t{bitscore}\n')

    pkl.dump(query_protein_alignment, open(f'{temp_pth}/query_protein_alignment.pkl', 'wb'))
    f.close()

    genome_info = {}
    protein_info = {}
    for record in SeqIO.parse(f'{temp_pth}/query_protein.faa', 'fasta'):
        protein_id = record.id
        genome_id = protein_id.rsplit('_', 1)[0]
        description = record.description.split(' # ')
        start = int(description[1])
        end = int(description[2])
        strand = int(description[3])

        if genome_id not in genome_info:
            genome_info[genome_id] = []
            protein_info[genome_id] = {}

        if genome_id in query_protein_alignment and protein_id in query_protein_alignment[genome_id]:
            representative = query_protein_alignment[genome_id][protein_id][0][1]
            genome_info[genome_id].append([protein_id, representative])
            protein_info[genome_id][protein_id] = {'length': len(record.seq), 'start': start, 'end': end, 'strand': strand,
                                                   'cluster': representative}
        else:
            genome_info[genome_id].append([protein_id, f'UNK_{protein_id.rsplit("_", 1)[-1]}'])
            protein_info[genome_id][protein_id] = {'length': len(record.seq), 'start': start, 'end': end,
                                                   'strand': strand,
                                                   'cluster': 'NA'}

    pkl.dump(genome_info, open(f'{temp_pth}/query_seq_info.pkl', 'wb'))
    pkl.dump(protein_info, open(f'{temp_pth}/query_protein_info.pkl', 'wb'))


def parse_alignment_result_bin(database_pth, out_pth):
    temp_pth = os.path.join(out_pth, 'midfolder')
    reference_protein_cluster = pkl.load(open(f'{database_pth}/reference_protein_cluster.pkl', 'rb'))

    f = open(f'{temp_pth}/query_protein_blastp.tsv', 'w')
    f.write('Contig_id\tQuery\tSubject\tCluster\tIdentity\tE-value\tBitscore\n')

    # parse diamond result
    query_protein_alignment = {'bin': {}}
    for line in open(f'{temp_pth}/query_protein.blastp', 'r'):
        line = line.strip().split('\t')
        query = line[0]
        query_contig_id = query.rsplit('_', 1)[0]
        subject = line[1]
        cluster = reference_protein_cluster[subject]
        identity = float(line[2])
        evalue = float(line[10])
        bitscore = float(line[11])

        if query not in query_protein_alignment['bin']:
            query_protein_alignment['bin'][query] = []
        query_protein_alignment['bin'][query].append([subject, cluster, identity, evalue, bitscore])
        f.write(f'{query_contig_id}\t{query}\t{subject}\t{reference_protein_cluster[subject]}\t{identity}\t{evalue}\t{bitscore}\n')

    pkl.dump(query_protein_alignment, open(f'{temp_pth}/query_protein_alignment.pkl', 'wb'))
    f.close()

    genome_info = {}
    protein_info = {'bin': {}}
    for record in SeqIO.parse(f'{temp_pth}/query_protein.faa', 'fasta'):
        protein_id = record.id
        genome_id = protein_id.rsplit('_', 1)[0]
        description = record.description.split(' # ')
        start = int(description[1])
        end = int(description[2])
        strand = int(description[3])

        if genome_id not in genome_info:
            genome_info[genome_id] = []

        if protein_id in query_protein_alignment['bin']:
            representative = query_protein_alignment['bin'][protein_id][0][1]
            genome_info[genome_id].append([protein_id, representative])
            protein_info['bin'][protein_id] = {'length': len(record.seq), 'start': start, 'end': end, 'strand': strand,
                                                   'cluster': representative}
        else:
            genome_info[genome_id].append([protein_id, f'UNK_{protein_id.rsplit("_", 1)[-1]}'])
            protein_info['bin'][protein_id] = {'length': len(record.seq), 'start': start, 'end': end,
                                                   'strand': strand,
                                                   'cluster': 'NA'}

    assembly_seq_info = {'bin': []}
    for genome_id in genome_info:
        assembly_seq_info['bin'].extend(genome_info[genome_id])

    pkl.dump(assembly_seq_info, open(f'{temp_pth}/query_seq_info.pkl', 'wb'))
    pkl.dump(protein_info, open(f'{temp_pth}/query_protein_info.pkl', 'wb'))


def compute_bit_aai_score(out_pth):
    temp_pth = os.path.join(out_pth, 'midfolder')
    protein_info = pkl.load(open(f'{temp_pth}/query_protein_info.pkl', 'rb'))

    alignment_scores = {}
    query_protein_alignment = pkl.load(open(f'{temp_pth}/query_protein_alignment.pkl', 'rb'))

    for contig in query_protein_alignment:

        contig_alignments = query_protein_alignment[contig]
        contig_amino_acid_len = sum([protein_info[contig][protein]['length'] for protein in protein_info[contig]])
        alignment_scores[contig] = {}

        for protein in contig_alignments:

            protein_len = protein_info[contig][protein]['length']
            protein_alignments = contig_alignments[protein]
            bit_scores = {}
            aai_scores = {}

            for record in protein_alignments:
                subject = record[0]
                genome_id = subject.rsplit('_', 1)[0]
                identity = record[2]
                bitscore = record[4]

                bit_scores[genome_id] = max(bit_scores.get(genome_id, 0), bitscore)
                aai_scores[genome_id] = max(aai_scores.get(genome_id, 0), identity)

            for genome_id in bit_scores:
                if genome_id not in alignment_scores[contig]:
                    alignment_scores[contig][genome_id] = {'bit_score': [], 'aai_score': [], 'amino_acid_fraction': []}

                alignment_scores[contig][genome_id]['bit_score'].append(bit_scores[genome_id])
                alignment_scores[contig][genome_id]['aai_score'].append(aai_scores[genome_id] * protein_len)
                alignment_scores[contig][genome_id]['amino_acid_fraction'].append(protein_len)

        for genome_id in alignment_scores[contig]:
            alignment_scores[contig][genome_id]['bit_score'] = np.sum(alignment_scores[contig][genome_id]['bit_score'])
            alignment_scores[contig][genome_id]['aai_score'] = np.sum(alignment_scores[contig][genome_id]['aai_score']) / np.sum(alignment_scores[contig][genome_id]['amino_acid_fraction'])
            alignment_scores[contig][genome_id]['amino_acid_fraction'] = np.sum(alignment_scores[contig][genome_id]['amino_acid_fraction']) * 100 / contig_amino_acid_len

    pkl.dump(alignment_scores, open(f'{temp_pth}/query_alignment_score.pkl', 'wb'))


def compute_bit_aai_score_bin(database_pth, out_pth):
    temp_pth = os.path.join(out_pth, 'midfolder')
    protein_info = pkl.load(open(f'{temp_pth}/query_protein_info.pkl', 'rb'))
    genome2assembly = pkl.load(open(f'{database_pth}/reference_genome2assembly.pkl', 'rb'))

    alignment_scores = {'bin': {}}
    query_protein_alignment = pkl.load(open(f'{temp_pth}/query_protein_alignment.pkl', 'rb'))

    for contig in query_protein_alignment:

        contig_alignments = query_protein_alignment[contig]
        contig_amino_acid_len = sum([protein_info[contig][protein]['length'] for protein in protein_info[contig]])

        for protein in contig_alignments:

            protein_len = protein_info[contig][protein]['length']
            protein_alignments = contig_alignments[protein]
            bit_scores = {}
            aai_scores = {}

            for record in protein_alignments:
                subject = record[0]
                genome_id = genome2assembly[subject.rsplit('_', 1)[0]]
                identity = record[2]
                bitscore = record[4]

                bit_scores[genome_id] = max(bit_scores.get(genome_id, 0), bitscore)
                aai_scores[genome_id] = max(aai_scores.get(genome_id, 0), identity)

            for genome_id in bit_scores:
                if genome_id not in alignment_scores[contig]:
                    alignment_scores[contig][genome_id] = {'bit_score': [], 'aai_score': [], 'amino_acid_fraction': []}

                alignment_scores[contig][genome_id]['bit_score'].append(bit_scores[genome_id])
                alignment_scores[contig][genome_id]['aai_score'].append(aai_scores[genome_id] * protein_len)
                alignment_scores[contig][genome_id]['amino_acid_fraction'].append(protein_len)

        for genome_id in alignment_scores[contig]:
            alignment_scores[contig][genome_id]['bit_score'] = np.sum(alignment_scores[contig][genome_id]['bit_score'])
            alignment_scores[contig][genome_id]['aai_score'] = np.sum(alignment_scores[contig][genome_id]['aai_score']) / np.sum(alignment_scores[contig][genome_id]['amino_acid_fraction'])
            alignment_scores[contig][genome_id]['amino_acid_fraction'] = np.sum(alignment_scores[contig][genome_id]['amino_acid_fraction']) * 100 / contig_amino_acid_len

    pkl.dump(alignment_scores, open(f'{temp_pth}/query_alignment_score.pkl', 'wb'))


def compute_confidence(num_ref, num_cluster, num_singleton, num_common, num_a, num_b):
    max_sig = 1000

    total_pcs = num_cluster + num_singleton
    num_a, num_b = sorted([num_a, num_b])
    pval = stats.hypergeom.sf(num_common - 1, total_pcs, num_a, num_b)
    sig = min(max_sig, np.nan_to_num(-np.log10(pval) - np.log10(num_ref)))
    if sig >= 5:
        confidence = 'high'
    elif sig >= 1:
        confidence = 'medium'
    else:
        confidence = 'low'
    return sig, confidence

def get_protein_set(seq_info):
    protein_set = set()
    num_singleton = 0
    for protein in seq_info:
        protein_set.add(protein[1])
        if isinstance(protein[1], str):
            num_singleton += 1
    return protein_set, num_singleton


def compute_sv(database_pth, out_pth):
    temp_pth = os.path.join(out_pth, 'midfolder')
    query_seq_info = pkl.load(open(f'{temp_pth}/query_seq_info.pkl', 'rb'))
    query_seq_list = list(query_seq_info.keys())
    query_alignment_score = pkl.load(open(f'{temp_pth}/query_alignment_score.pkl', 'rb'))

    reference_genome_info = pkl.load(open(f'{database_pth}/reference_genome_info.pkl', 'rb'))
    reference_genome_list = list(reference_genome_info.keys())

    run_time = []
    result = {}
    for i, query in enumerate(query_seq_list):
        start_t = time.perf_counter()
        result[query] = []
        q_proteins, num_singleton = get_protein_set(query_seq_info[query])

        max_shared_protein = 0
        for j, ref in enumerate(reference_genome_list):

            r_proteins, _ = get_protein_set(reference_genome_info[ref])
            common_proteins = q_proteins & r_proteins

            #if len(common_proteins) >= max(1, int(max_shared_protein * 0.5)):
                #max_shared_protein = max(max_shared_protein, len(common_proteins))
            if len(common_proteins) > 0:

                if len(common_proteins) == 1:
                    events = {
                        'match': 1,
                        'mutation': 0,
                        'insertion': 0,
                        'deletion': 0,
                        'translocation': 0,
                        'duplication': 0,
                        'outlier': len(q_proteins) - 1,
                        'gapopen': 0
                    }
                else:
                    g = GenomeStructure(query, query_seq_info[query], ref, reference_genome_info[ref], common_proteins)
                    events = g.get_events()

                #sig, confidence = compute_confidence(num_singleton, len(common_proteins), len(q_proteins), len(r_proteins))
                #sig, confidence = 0, 'high'
                alignment_scores = query_alignment_score[query].get(ref, {'bit_score': 0, 'aai_score': 0, 'amino_acid_fraction': 0})
                bit_score = alignment_scores['bit_score']
                aai_score = alignment_scores['aai_score']
                amino_acid_fraction = alignment_scores['amino_acid_fraction']
                protein_faction = len(common_proteins) * 100 / len(q_proteins)
                result[query].append([query, len(q_proteins), ref, len(r_proteins), len(common_proteins),
                                      events, bit_score, aai_score, protein_faction, amino_acid_fraction])


        end_t = time.perf_counter()
        run_time.append(end_t - start_t)

    pkl.dump(result, open(f'{temp_pth}/result_info.pkl', 'wb'))


def compute_sv_bin(database_pth, out_pth):
    temp_pth = os.path.join(out_pth, 'midfolder')
    query_seq_info = pkl.load(open(f'{temp_pth}/query_seq_info.pkl', 'rb'))
    query_seq_list = list(query_seq_info.keys())
    query_alignment_score = pkl.load(open(f'{temp_pth}/query_alignment_score.pkl', 'rb'))

    reference_genome_info = pkl.load(open(f'{database_pth}/reference_assembly_info.pkl', 'rb'))
    reference_genome_list = list(reference_genome_info.keys())

    run_time = []
    result = {}
    for i, query in enumerate(query_seq_list):
        start_t = time.perf_counter()
        result[query] = []
        q_proteins, num_singleton = get_protein_set(query_seq_info[query])

        max_shared_protein = 0
        for j, ref in enumerate(reference_genome_list):

            r_proteins, _ = get_protein_set(reference_genome_info[ref])
            common_proteins = q_proteins & r_proteins

            #if len(common_proteins) >= max(1, int(max_shared_protein * 0.5)):
                #max_shared_protein = max(max_shared_protein, len(common_proteins))
            if len(common_proteins) > 0:

                if len(common_proteins) == 1:
                    events = {
                        'match': 1,
                        'mutation': 0,
                        'insertion': 0,
                        'deletion': 0,
                        'translocation': 0,
                        'duplication': 0,
                        'outlier': len(q_proteins) - 1,
                        'gapopen': 0
                    }
                else:
                    g = GenomeStructure(query, query_seq_info[query], ref, reference_genome_info[ref], common_proteins)
                    events = g.get_events()

                #sig, confidence = compute_confidence(num_singleton, len(common_proteins), len(q_proteins), len(r_proteins))
                #sig, confidence = 0, 'high'
                alignment_scores = query_alignment_score[query].get(ref, {'bit_score': 0, 'aai_score': 0, 'amino_acid_fraction': 0})
                bit_score = alignment_scores['bit_score']
                aai_score = alignment_scores['aai_score']
                amino_acid_fraction = alignment_scores['amino_acid_fraction']
                protein_faction = len(common_proteins) * 100 / len(q_proteins)
                result[query].append([query, len(q_proteins), ref, len(r_proteins), len(common_proteins),
                                      events, bit_score, aai_score, protein_faction, amino_acid_fraction])


        end_t = time.perf_counter()
        run_time.append(end_t - start_t)

    pkl.dump(result, open(f'{temp_pth}/result_info.pkl', 'wb'))


def compute_structural_score(events):
    a0 = 1
    a1 = -0.05
    a2 = -0.1
    a3 = -0.05
    a4 = -0.05
    a5 = -0.05
    a6 = -0.05

    return (a0 * events['match'] + a1 * events['mutation'] + a2 * events['gapopen'] + a3 * (events['insertion'] + events['deletion'] - events['gapopen']) + a4 * events['translocation'] + a5 * events['duplication'] + a6 * events['outlier'])


def compute_structural_score_bin(events):
    a0 = 1
    a1 = -0.05
    a2 = -0.1
    a3 = -0.05
    a4 = 0
    a5 = -0.05
    a6 = -0.05

    return (a0 * events['match'] + a1 * events['mutation'] + a2 * events['gapopen'] + a3 * (events['insertion'] + events['deletion'] - events['gapopen']) + a4 * events['translocation'] + a5 * events['duplication'] + a6 * events['outlier'])


def write_result(database_pth, out_pth):
    temp_pth = os.path.join(out_pth, 'midfolder')
    result_info = pkl.load(open(f'{temp_pth}/result_info.pkl', 'rb'))
    query_seq_len = pkl.load(open(f'{temp_pth}/query_seq_length.pkl', 'rb'))
    query_dtr = pkl.load(open(f'{temp_pth}/query_dtr.pkl', 'rb'))
    reference_genome_length = pkl.load(open(f'{database_pth}/reference_genome_length.pkl', 'rb'))
    reference_genome_tax = pkl.load(open(f'{database_pth}/reference_genome_tax.pkl', 'rb'))
    query_seq_info = pkl.load(open(f'{temp_pth}/query_seq_info.pkl', 'rb'))

    num_ref = len(reference_genome_tax)
    num_cluster = len(pkl.load(open(f'{database_pth}/reference_protein_cluster2id.pkl', 'rb')))

    f = open(f'{out_pth}/completeness_result.csv', 'w')
    f.write('seq_name,seq_length,num_protein,dtr,expected_length,completeness,confidence,significance_score,ref_id,ref_taxonomy,match,mutation,insertion,'
            'deletion,translocation,duplication,outlier\n')

    for query in result_info:
        records = result_info[query]
        if len(records) == 0:
            continue

        for record in records:
            events = record[5]
            struc_score = compute_structural_score(events)
            record.append(struc_score)

        # sort by structural score and bit score
        records = sorted(records, key=lambda x: (x[-1], x[6]), reverse=True)
        if records[0][-1] < 1:
            records = sorted(records, key=lambda x: x[6], reverse=True)

        # select the best alignment
        best_alignment = records[0][2]
        expect_genome_length = reference_genome_length[best_alignment]
        completeness = query_seq_len[query] * 100 / expect_genome_length
        completeness = min(100, completeness)

        q_proteins, num_singleton = get_protein_set(query_seq_info[query])
        sig, confidence = compute_confidence(num_ref, num_cluster, num_singleton, records[0][4], records[0][1],
                                             records[0][3])

        f.write(f'{query},{query_seq_len[query]},{records[0][1]},{query_dtr[query]},{expect_genome_length},{completeness},{confidence},{sig},{best_alignment},'
            f'{reference_genome_tax[best_alignment]},{records[0][4]},{records[0][5]["mutation"]},'
            f'{records[0][5]["insertion"]},{records[0][5]["deletion"]},{records[0][5]["translocation"]},'
            f'{records[0][5]["duplication"]},{records[0][5]["outlier"]}\n')

    f.close()


def write_result_bin(database_pth, out_pth):
    temp_pth = os.path.join(out_pth, 'midfolder')
    result_info = pkl.load(open(f'{temp_pth}/result_info.pkl', 'rb'))
    query_seq_len = pkl.load(open(f'{temp_pth}/query_seq_length.pkl', 'rb'))
    reference_genome_length = pkl.load(open(f'{database_pth}/reference_assembly_length.pkl', 'rb'))
    reference_genome_tax = pkl.load(open(f'{database_pth}/reference_assembly_tax.pkl', 'rb'))
    query_seq_info = pkl.load(open(f'{temp_pth}/query_seq_info.pkl', 'rb'))

    num_ref = len(reference_genome_tax)
    num_cluster = len(pkl.load(open(f'{database_pth}/reference_protein_cluster2id.pkl', 'rb')))

    f = open(f'{out_pth}/completeness_result.csv', 'w')
    f.write('seq_name,seq_length,num_protein,expected_length,completeness,confidence,significance_score,ref_id,ref_taxonomy,match,mutation,insertion,'
            'deletion,duplication,outlier\n')

    for query in result_info:
        records = result_info[query]
        if len(records) == 0:
            continue

        for record in records:
            events = record[5]
            struc_score = compute_structural_score_bin(events)
            record.append(struc_score)

        # sort by structural score and bit score
        records = sorted(records, key=lambda x: (x[-1], x[6]), reverse=True)
        if records[0][-1] < 1:
            records = sorted(records, key=lambda x: x[6], reverse=True)

        # select the best alignment
        best_alignment = records[0][2]
        expect_genome_length = reference_genome_length[best_alignment]
        completeness = query_seq_len[query] * 100 / expect_genome_length
        completeness = min(100, completeness)

        q_proteins, num_singleton = get_protein_set(query_seq_info[query])
        sig, confidence = compute_confidence(num_ref, num_cluster, num_singleton, records[0][4], records[0][1],
                                             records[0][3])

        f.write(f'{query},{query_seq_len[query]},{records[0][1]},{expect_genome_length},{completeness},{confidence},{sig},{best_alignment},'
            f'{reference_genome_tax[best_alignment]},{records[0][4]},{records[0][5]["mutation"]},'
            f'{records[0][5]["insertion"]},{records[0][5]["deletion"]},'
            f'{records[0][5]["duplication"]},{records[0][5]["outlier"]}\n')

    f.close()


def completeness(input, db, output, threads, bin):
    if not os.path.exists(output):
        os.makedirs(output)

    if not os.path.exists(os.path.join(output, 'midfolder')):
        os.makedirs(os.path.join(output, 'midfolder'))

    print("Running ViralQC completeness estimation...")
    if bin:
        get_query_length_bin(input, output)
        protein_prediction_and_alignment(db, input, output, threads)
        parse_alignment_result_bin(db, output)
        compute_bit_aai_score_bin(db, output)
        print("[3/4] Computing structural variants...")
        compute_sv_bin(db, output)
        print("[4/4] Writing results...")
        write_result_bin(db, output)
    else:
        get_query_length_dtr(input, output)
        protein_prediction_and_alignment(db, input, output, threads)
        parse_alignment_result(db, output)
        compute_bit_aai_score(output)
        print("[3/4] Computing structural variants...")
        compute_sv(db, output)
        print("[4/4] Writing results...")
        write_result(db, output)

    max_mem_self = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    max_mem_child = resource.getrusage(resource.RUSAGE_CHILDREN).ru_maxrss
    if platform.system() == "Linux":
        peak_mem =  (max_mem_self + max_mem_child) / float(1e6)
    else:
        peak_mem =  (max_mem_self + max_mem_child) / float(1e9)
    print(f"Peak mem: {round(peak_mem, 2)} GB")






