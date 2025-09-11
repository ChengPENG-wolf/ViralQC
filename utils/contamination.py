import os
import subprocess
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
import pickle as pkl
import pandas as pd
import numpy as np
import re


def special_match(strg, search=re.compile(r'[^ACGT]').search):
    return not bool(search(strg))

def get_input_length(input_pth, output_pth):
    temp_pth = os.path.join(output_pth, 'midfolder')
    contig_len = {}
    for record in SeqIO.parse(input_pth, 'fasta'):
        contig_len[record.id] = len(record.seq)
    pkl.dump(contig_len, open(f"{temp_pth}/input_len.pkl", 'wb'))


def predict_protein(input_pth, output_pth, threads):
    temp_pth = os.path.join(output_pth, 'midfolder')

    cmd = f"python utils/parallel-prodigal-gv.py -t {threads} -q -i {input_pth} -a {temp_pth}/query_protein.faa"
    _ = subprocess.check_call(cmd, shell=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)

    r = []
    for record in SeqIO.parse(f"{temp_pth}/query_protein.faa", 'fasta'):
        if str(record.seq)[-1] == '*':
            r.append(SeqRecord(Seq(str(record.seq)[:-1]), id=record.id, description=''))
        else:
            r.append(SeqRecord(record.seq, id=record.id, description=''))
    SeqIO.write(r, f"{temp_pth}/query_protein_temp.faa", 'fasta')


def protein_branch(db, output_pth):
    temp_pth = os.path.join(output_pth, 'midfolder')
    cmd = f"python utils/extract.py esm2_t33_650M_UR50D {temp_pth}/query_protein_temp.faa {temp_pth} --include mean"
    _ = subprocess.check_call(cmd, shell=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)

    cmd = f"python utils/protein_branch.py --input {temp_pth}/embed.pt --db {db}/protein_model --output {temp_pth}"
    _ = subprocess.check_call(cmd, shell=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)

    os.remove(f"{temp_pth}/query_protein_temp.faa")


def parse_protein_result(output_pth):
    temp_pth = os.path.join(output_pth, 'midfolder')
    protein_classifcy_result = pkl.load(open(f"{temp_pth}/protein_result.pkl", 'rb'))

    f = open(f"{temp_pth}/protein_result.csv", 'w')
    f.write('protein_id,start,end,score,predict\n')

    for record in SeqIO.parse(f'{temp_pth}/query_protein.faa', 'fasta'):
        protein_name = record.id
        description = record.description.split(' # ')
        start = int(description[1])
        end = int(description[2])
        contig_name = protein_name.rsplit('_', 1)[0]
        if protein_name in protein_classifcy_result:
            if protein_classifcy_result[protein_name] < 0.5:
                f.write(f'{protein_name},{start},{end},{protein_classifcy_result[protein_name]},non-virus\n')
            else:
                f.write(f'{protein_name},{start},{end},{protein_classifcy_result[protein_name]},virus\n')
    f.close()


def dna_branch(input_pth, db, output_pth):
    temp_pth = os.path.join(output_pth, 'midfolder')

    records = []
    frag_len = 2000
    for record in SeqIO.parse(input_pth, 'fasta'):
        sequence = str(record.seq).upper()
        if len(sequence) >= frag_len:
            for i in range(0, len(sequence) - frag_len + 1, 500):
                sequence1 = sequence[i:i + frag_len]
                if special_match(sequence1):
                    records.append(SeqRecord(Seq(sequence1), id=f'{record.id}_{i+1}_{i+frag_len}', description=''))
    SeqIO.write(records, f"{temp_pth}/query_seq_{frag_len}.fasta", 'fasta')

    cmd = f"python utils/dna_branch.py --input {temp_pth}/query_seq_{frag_len}.fasta --db {db}/dna_model --output {temp_pth}/dna_branch"
    _ = subprocess.check_call(cmd, shell=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)


def branch_aggregattion(output_pth, virus_threshold_p, host_threshold_p, virus_threshold_d):
    temp_pth = os.path.join(output_pth, 'midfolder')
    protein_classifcy_result = pkl.load(open(f"{temp_pth}/protein_result.pkl", 'rb'))
    contig_len = pkl.load(open(f"{temp_pth}/input_len.pkl", 'rb'))

    # get protein result
    protein_result = {}
    for record in SeqIO.parse(f'{temp_pth}/query_protein.faa', 'fasta'):
        protein_name = record.id
        description = record.description.split(' # ')
        start = int(description[1])
        end = int(description[2])
        contig_name = protein_name.rsplit('_', 1)[0]
        if protein_name in protein_classifcy_result:
            if protein_classifcy_result[protein_name] < host_threshold_p:
                if contig_name not in protein_result:
                    protein_result[contig_name] = []
                protein_result[contig_name].append({'start': start, 'end': end, 'predict': 'non-virus'})
            elif protein_classifcy_result[protein_name] > virus_threshold_p:
                if contig_name not in protein_result:
                    protein_result[contig_name] = []
                protein_result[contig_name].append({'start': start, 'end': end, 'predict': 'virus'})
            else:
                if contig_name not in protein_result:
                    protein_result[contig_name] = []
                protein_result[contig_name].append({'start': start, 'end': end, 'predict': 'unknown'})


    # get dna result
    dna_result = {}
    df = pd.read_csv(f"{temp_pth}/dna_branch/dna_result.csv")
    for i, row in df.iterrows():
        segment_name = row['seq_name']
        virus_score = float(row['virus_score'])
        contig_name = segment_name.rsplit('_', 2)[0]
        start = int(segment_name.rsplit('_', 2)[1])
        end = int(segment_name.rsplit('_', 2)[2])

        if virus_score > virus_threshold_d:
            if contig_name not in dna_result:
                dna_result[contig_name] = []
            dna_result[contig_name].append({'start': start, 'end': end, 'predict': 'virus'})


    for contig_name in dna_result:

        y = np.zeros(contig_len[contig_name] + 1)
        for i in range(len(dna_result[contig_name])):
            start = dna_result[contig_name][i]['start']
            end = dna_result[contig_name][i]['end']
            predict = 1 if dna_result[contig_name][i]['predict'] == 'virus' else -1
            y[start:end + 1] += predict

        new_result = []
        new_result.append({'start': 1, 'end': 500, 'predict': y[1]})
        for i in range(501, len(y), 500):
            if y[i] == new_result[-1]['predict']:
                new_result[-1]['end'] = min(i+499, contig_len[contig_name])
            else:
                new_result.append({'start': i, 'end': min(i+499, contig_len[contig_name]), 'predict': y[i]})

        dna_result[contig_name] = new_result

    for contig_name in dna_result:
        new_result = []
        for i in range(0, len(dna_result[contig_name])):

            if dna_result[contig_name][i]['predict'] >= 3:
                if len(new_result) > 0 and new_result[-1]['end'] == dna_result[contig_name][i]['start'] - 1:
                    new_result[-1]['end'] = dna_result[contig_name][i]['end']
                else:
                    new_result.append(dna_result[contig_name][i])
                    new_result[-1]['predict'] = 'virus'
        dna_result[contig_name] = new_result

    for contig_name in dna_result:
        for i in range(0, len(dna_result[contig_name])):
            if dna_result[contig_name][i]['predict'] == 'virus':
                start = dna_result[contig_name][i]['start']
                end = dna_result[contig_name][i]['end']
                if start == 1001:
                    dna_result[contig_name][i]['start'] = 1
                else:
                    dna_result[contig_name][i]['start'] = start - 500
                if end == contig_len[contig_name] - 1000:
                    dna_result[contig_name][i]['end'] = contig_len[contig_name]


    # merge protein and dna result
    result = {}
    for contig_name in protein_result:
        seq = np.zeros(contig_len[contig_name] + 1)
        result[contig_name] = []

        if contig_name in dna_result:
            for segment in dna_result[contig_name]:
                seq[segment['start']:segment['end']+1] += 1

        for i in range(len(protein_result[contig_name])):
            if protein_result[contig_name][i]['predict'] == 'non-virus':
                start = protein_result[contig_name][i]['start']
                end = protein_result[contig_name][i]['end']
                length = end - start + 1
                if sum(seq[start:end+1]) >= length * 0.5:
                    result[contig_name].append({'start': start, 'end': end, 'predict': 'virus'})
                else:
                    result[contig_name].append({'start': start, 'end': end, 'predict': 'non-virus'})
            else:
                result[contig_name].append(protein_result[contig_name][i])


    for contig_name in result:

        seq = np.zeros(len(result[contig_name]))
        scores = np.zeros(len(result[contig_name]) - 1)
        for i in range(len(result[contig_name])):
            if result[contig_name][i]['predict'] == 'virus':
                seq[i] = 1
            elif result[contig_name][i]['predict'] == 'non-virus':
                seq[i] = -1

        window_size = min(50, max(10, int(len(result[contig_name]) * 0.1)))
        flag = True
        breakpoint = []
        while flag:
            if len(breakpoint) == 0:
                start_pos = 0
            else:
                start_pos = breakpoint[-1] + 1

            if start_pos >= len(result[contig_name]) - 1:
                break

            for i in range(start_pos, len(result[contig_name]) - 1):
                s_left = sum(seq[start_pos:i+1]) / (i - start_pos + 1)
                s_right = sum(seq[i+1:i+1+window_size]) / len(seq[i+1:i+1+window_size])
                scores[i] = abs(s_left - s_right)

            max_score = np.max(scores[start_pos:])
            max_pos = np.argmax(scores[start_pos:]) + start_pos

            if max_score > 1:
                breakpoint.append(max_pos)
            else:
                flag = False

        new_result = []
        for i in range(len(breakpoint)):
            if i == 0:
                predict = 'virus' if sum(seq[0:breakpoint[i]+1]) > sum(seq[breakpoint[i]+1:breakpoint[i]+1+window_size]) else 'non-virus'
                new_result.append({'start': result[contig_name][0]['start'], 'end': result[contig_name][breakpoint[i]]['end'], 'predict': predict})
            else:
                predict = 'virus' if sum(seq[breakpoint[i-1]+1:breakpoint[i]+1]) > sum(seq[breakpoint[i]+1:breakpoint[i]+1+window_size]) else 'non-virus'
                new_result.append({'start': result[contig_name][breakpoint[i-1]+1]['start'], 'end': result[contig_name][breakpoint[i]]['end'], 'predict': predict})

        if len(breakpoint) > 0:
            if breakpoint[-1] < len(result[contig_name]) - 1:
                if len(breakpoint) == 1:
                    predict = 'non-virus' if sum(seq[0:breakpoint[-1]+1]) > sum(seq[breakpoint[-1]+1:]) else 'virus'
                    new_result.append({'start': result[contig_name][breakpoint[-1] + 1]['start'],
                                       'end': result[contig_name][-1]['end'], 'predict': predict})
                else:
                    predict = 'non-virus' if sum(seq[breakpoint[-2]+1:breakpoint[-1]+1]) > sum(seq[breakpoint[-1]+1:]) else 'virus'
                    new_result.append({'start': result[contig_name][breakpoint[-1]+1]['start'], 'end': result[contig_name][-1]['end'], 'predict': predict})
        else:
            if np.mean(seq) <= -0.5:
                new_result.append({'start': result[contig_name][0]['start'], 'end': result[contig_name][-1]['end'], 'predict': 'non-virus'})
        result[contig_name] = new_result


    # extend boundary to contig boundary
    for contig_name in result:
        new_result = []
        for i in range(len(result[contig_name])):
            if i == 0:
                new_result.append({'start': 1, 'end': result[contig_name][i]['end'],'predict': result[contig_name][i]['predict']})
            else:
                if result[contig_name][i]['predict'] == new_result[-1]['predict']:
                    new_result[-1]['end'] = result[contig_name][i]['end']

                else:
                    if new_result[-1]['predict'] == 'virus':
                        new_result[-1]['end'] = max(new_result[-1]['end'], result[contig_name][i]['start'] - 1)
                        new_result.append({'start': new_result[-1]['end'] + 1,
                                           'end': result[contig_name][i]['end'],
                                           'predict': result[contig_name][i]['predict']})
                    else:
                        new_result.append({'start': new_result[-1]['end'] + 1,
                                           'end': result[contig_name][i]['end'],
                                           'predict': result[contig_name][i]['predict']})

        if len(new_result) != 0:
            new_result[-1]['end'] = contig_len[contig_name]

        result[contig_name] = new_result


    # record the number of proteins and length information of each region
    for contig_name in result:
        new_result = []
        for i in range(len(result[contig_name])):
            left = result[contig_name][i]['start']
            right = result[contig_name][i]['end']
            region_length = right - left + 1
            predict = result[contig_name][i]['predict']
            total_protein = 0
            viral_protein = 0
            non_viral_protein = 0
            viral_region = 0

            if contig_name in protein_result:
                for segment in protein_result[contig_name]:
                    if (left <= segment['start'] <= right) or (left <= segment['end'] <= right):
                        if min(segment['end'], right) - max(segment['start'], left) + 1 >= (segment['end'] - segment['start'] + 1) * 0.5:
                            total_protein += 1
                            if segment['predict'] == 'virus':
                                viral_protein += 1
                            elif segment['predict'] == 'non-virus':
                                non_viral_protein += 1

            if contig_name in dna_result:
                for segment in dna_result[contig_name]:
                    if (left <= segment['start'] <= right) or (left <= segment['end'] <= right):
                        viral_region += min(segment['end'], right) - max(segment['start'], left) + 1

            if len(new_result) != 0 and new_result[-1]['predict'] == predict:
                new_result[-1]['end'] = right
                new_result[-1]['length'] += region_length
                new_result[-1]['total_protein'] += total_protein
                new_result[-1]['viral_protein'] += viral_protein
                new_result[-1]['non_viral_protein'] += non_viral_protein
                new_result[-1]['viral_region'] += viral_region
            else:
                new_result.append({'start': left, 'end': right, 'length': region_length, 'predict': predict,
                                   'total_protein': total_protein, 'viral_protein': viral_protein,
                                   'non_viral_protein': non_viral_protein, 'viral_region': viral_region})
        result[contig_name] = new_result

    new_result = {}
    for contig_name in result:
        if len(result[contig_name]) == 0:
            continue
        records = []
        for i in range(len(result[contig_name])):

            records.append(result[contig_name][i])
        new_result[contig_name] = records
    result = new_result

    pkl.dump(result, open(f"{temp_pth}/contamination_result.pkl", 'wb'))
    f = open(f"{temp_pth}/provirus_result.csv", 'w')
    f.write('accession,region_length,predict,predict left,predict right,#total protein,'
            '#viral protein,#non-viral protein,viral region\n')
    for contig_name in result:
        for segment in result[contig_name]:
            region_length = segment['length']
            predict = segment['predict']
            left = segment['start']
            right = segment['end']
            total_protein = segment['total_protein']
            viral_protein = segment['viral_protein']
            non_viral_protein = segment['non_viral_protein']
            viral_region = segment['viral_region']

            f.write(f'{contig_name},{region_length},{predict},{left},{right},{total_protein},{viral_protein},'
                    f'{non_viral_protein},{viral_region}\n')

    f.close()


def write_result(input_pth, output_pth):
    temp_pth = os.path.join(output_pth, 'midfolder')
    contig_len = pkl.load(open(f"{temp_pth}/input_len.pkl", 'rb'))
    result = pkl.load(open(f"{temp_pth}/contamination_result.pkl", 'rb'))

    f = open(f"{output_pth}/contamination_result.csv", 'w')
    f.write(f'seq_name,seq_length,total_genes,virus_genes,non-virus_genes,virus_length,non-virus_length,region_types,region_coords_bp\n')

    for contig_id in result:
        if len(result[contig_id]) <= 1 and result[contig_id][0]['predict'] == 'virus':
            continue
        contig_length = contig_len[contig_id]
        total_genes = 0
        virus_genes = 0
        non_virus_genes = 0
        virus_length = 0
        non_virus_length = 0
        region_types = ''
        region_coords_bp = ''
        for region in result[contig_id]:
            total_genes += region['total_protein']
            virus_genes += region['viral_protein']
            non_virus_genes += region['non_viral_protein']
            if region['predict'] == 'virus':
                virus_length += region['length']
            else:
                non_virus_length += region['length']
            region_types += f'{region["predict"]};'
            region_coords_bp += f'{region["start"]}-{region["end"]};'
        f.write(f'{contig_id},{contig_length},{total_genes},{virus_genes},{non_virus_genes},{virus_length},{non_virus_length},{region_types},{region_coords_bp}\n')
    f.close()

    decontaminated_records = []
    for record in SeqIO.parse(input_pth, 'fasta'):
        if record.id in result:
            if len(result[record.id]) <= 1:
                decontaminated_records.append(record)
                continue
            for region in result[record.id]:
                if region['predict'] == 'virus':
                    decontaminated_records.append(SeqRecord(record.seq[region['start']-1:region['end']], id=f'{record.id}_{region["start"]}_{region["end"]}', description=''))
        else:
            decontaminated_records.append(record)
    SeqIO.write(decontaminated_records, f"{output_pth}/extracted_virus.fasta", 'fasta')


def contamination(input, db, output, threads):
    if not os.path.exists(output):
        os.makedirs(output)

    if not os.path.exists(os.path.join(output, 'midfolder')):
        os.makedirs(os.path.join(output, 'midfolder'))

    get_input_length(input, output)
    predict_protein(input, output, threads)
    protein_branch(db, output)
    parse_protein_result(output)
    dna_branch(input, db, output)
    branch_aggregattion(output, virus_threshold_p=0.5, host_threshold_p=0.1, virus_threshold_d=0.5)
    write_result(input, output)




