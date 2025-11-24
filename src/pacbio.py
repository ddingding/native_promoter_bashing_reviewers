import os
import pandas as pd
from Bio import SeqIO
from Bio import Align
from tqdm import tqdm
import matplotlib.pyplot as plt
import numpy as np

import pysam
from collections import Counter
import re

from itertools import groupby
from operator import itemgetter
import nucleotides as nt
from nucleotides import get_ind_muts
from constants import sample_to_positions, sample_to_length, sample_to_5prime_bases

# create fastq file quality mapping
ordered_fastq_vals = '''!"#$%&'()*+,-./0123456789:;<=>?@ABCDEFGHI'''
FASTQSCORE_TO_NUM = dict(zip(ordered_fastq_vals, list(range(len(ordered_fastq_vals)))))

# todo: move to config file
# 2 means the plasmid only has the promoter and 5'utr region, but is driving the expression of a GFP protein. 
# 3 stands for 'full length' meaning the whole gene is found in that plasmid, ie. the promoter+5'utr is driving the native gene (wih introns and 2kb after the gene)
DIC_GENE_NAMES = {  
    '2A': 'SBPase',
    '2B': 'Raf1',
    '2C': 'Psbs',
    '3A': 'SBPase_fl',
    '3B': 'Raf1_fl',
    '3C': 'Psbs_fl',
}

# list_vals indicate: where does promoter start, where does promoter end, where does 5'UTR end
DIC_GENE_STR_SYN = {
    '2A': [2211,4211,4498],     # 2A:SBPase starts at 2211 until 4211, 5'UTR until 4498 until GFP starts
    '2B': [2196,4196,4280],     # 2B: Raf1 starts at 2196 until 4196, 5"UTR until 4280 when GFP starts
    '2C': [2196,4196,4396],     # 2C: Psbs: 2196 until 4196, 5" until 4396
}

def get_muts(r):
    # expects a pysam.AlignedSegment object
    # returns a mut_str with insertions, deletions and SNPs
    # TODO this misses short reads that don't span the full length, 
    # ie. anything that goes after the short read is not reported

    query_n = r.query_name # read sequence
    query_seq = r.query_sequence # reference sequence to align against

    try:
        aligned_tuples = r.get_aligned_pairs(with_seq=True)
    except ValueError:
        return 'noMDtag'
    
    list_muts = []

    last_query_pos = -1
    for base_tuple in aligned_tuples:
        query_pos, ref_pos, ref_base = base_tuple

        if ref_pos != None:
            # SNPs: lower case, and fetch from query sequence
            # works for barcode too
            if ref_base.islower():
                mut = ref_base.upper() + str(ref_pos) + query_seq[query_pos]
                list_muts.append(mut)
            
            # deletions in query show up as None in the read, but ref position stays.
            # works for the one deletion
            # note that each individual character will be marked
            if query_pos == None:
                mut = ref_base + str(ref_pos) + '_'
                list_muts.append(mut)
            
            # insertions are not flagged, show up as a skip in the sequence of tuples, 
            # and need to be reconstructed from the query_sequence, 
            if last_query_pos != None and query_pos != None:
                if int(query_pos) != (int(last_query_pos)+1):
                    #print('insertion')
                    mut = 'I'+ str(ref_pos) + query_seq[last_query_pos+1:query_pos]
                    list_muts.append(mut)

            last_query_pos = query_pos

    return ':'.join(list_muts)

def get_bc(str_muts):
    # expects a string of mutations separated by ':'
    # extracts a particular position
    
    return

def summarize_sam(samfile, file_out_name):
    # expects an aligned sam file, and outputs a mutation summary of mapped mutations
    # for the reads that were mapped fine, extract mutations and get the barcodes
    samfile = pysam.AlignmentFile(samfile, 'r')
    with open(file_out_name, 'w') as fout:
        for r in tqdm(samfile.fetch()):
            str_muts = get_muts(r)
            fout.write('\t'.join([r.query_name,str_muts])+ '\n')

def tuple_muts_to_mutstr(list_tuples):
    # colelct a tuple of mtuations into a mutstr
    list_muts = [''.join(map(str,list(t))) for t in list_tuples]
    return ':'.join(list_muts)

def get_bc_from_muts(muts, pos_lim_min=140, pos_lim_max=180):
    # get bc from mutations, some are coded as insertions, other as SNPs, see examples, based on being in a defined window.

    # examples:
    # get_bc_from_muts('I153AAAAAAAAAAAACCGT:T2810_:G2811_:T4313_') returns ('AAAAAAAAAAAACCGT', 'T2810_:G2811_:T4313_')
    # get_bc_from_muts('I0AATT:N153T:N154A:N155C:N156T:N157A:N158T:N159T:N160A:N161G:N162G:N163G:N164T:N165T:N166T:N167A:C264A:G2483_:A2484_:G2485_:A3494G') returns ('TACTATTAGGGTTTA', 'I0AATT:C264A:G2483_:A2484_:G2485_:A3494G')
    # get_bc_from_muts2('I147TCTGTAGAGTACATT', pos_lim_min=140, pos_lim_max=180) returns ('TCTGTAGAGTACATT', 'wt')
    # get_bc_from_muts('N153T:N154A:N155C:N156T:N157A:N158T:N159T:N160A:N161G:N162G:N163G:N164T:N165T:N166T:N167A', pos_lim_min=140, pos_lim_max=180) returns ('TACTATTAGGGTTTA', 'wt')

    list_muts_sep = nt.get_ind_muts(muts)

    barcode_muts = []
    for m_tuple in list_muts_sep:
        # #catching those barcodes that are flagged as insertions in the right position
        # if m_tuple[1] > pos_lim_min and m_tuple[1] < pos_lim_max and m_tuple[0] =='I':

        #     bc = m_tuple[2]
        #     muts_wo_bc = ':'.join([m[0]+str(m[1])+m[2] for m in list_muts_sep if m != m_tuple])
        #     # catch wt cases where this muts_wo_bc is empty
        #     if len(muts_wo_bc) == 0:
        #         muts_wo_bc = ['wt']
        #     return bc, ':'.join(muts_wo_bc)

        # otherwise coded as N SNPs, so add them to the barcode muts list.
        if m_tuple[0] == 'N':
            barcode_muts.append(m_tuple)

    # filter for the correct position
    bc_right_pos = []
    for m_tuple in barcode_muts:
        pos = int(m_tuple[1])
        if pos >= pos_lim_min and pos <= pos_lim_max:
            bc_right_pos.append(m_tuple)

    # sort the barcodes by position to create barcode
    pos_to_base = dict(zip([t[1] for t in bc_right_pos], [t[2] for t in bc_right_pos]))
    bc = ''
    for p in sorted(pos_to_base.keys()):
        bc+=pos_to_base[p]

    # find the remaining mutations that are not classified as barcodes.
    muts_wo_bc = [m[0]+str(m[1])+m[2] for m in list_muts_sep if m not in bc_right_pos]
    if len(muts_wo_bc) == 0:
        muts_wo_bc = ['wt']
    return bc, ':'.join(muts_wo_bc)


def plot_bar_dic_tuples(s_to_num, label1, label2, label3=None, color2='maroon'):
    plt.figure()
    plt.bar(s_to_num.keys(), [v[0] for v in s_to_num.values()], label = label1)
    plt.bar(s_to_num.keys(), [v[1] for v in s_to_num.values()], color=color2, alpha=0.1, label=label2)

    if len(list(s_to_num.values())[0]) ==3:
        plt.bar(s_to_num.keys(), [v[2] for v in s_to_num.values()], color='orange', alpha=0.3, label=label3)

    plt.title(f'{label1} and {label2} reads')
    plt.ylabel('# reads')
    plt.xlabel('sample')
    plt.legend()
    plt.show()

def plot_along_gene(dic_vals, plot_n='counts along gene'):
    plt.figure()
    sorted_vals = sorted(dic_vals.items(), key= lambda x: x[0])
    sorted_x = [v[0] for v in sorted_vals]
    sorted_y = [v[1] for v in sorted_vals]
    plt.plot(sorted_x, sorted_y)
    plt.xlim([0,6000])
    plt.title(plot_n)
    plt.show()

    # zoomed in version
    plt.figure()
    sorted_vals = sorted(dic_vals.items(), key= lambda x: x[0])
    sorted_x = [v[0] for v in sorted_vals]
    sorted_y = [v[1] for v in sorted_vals]
    plt.plot(sorted_x, sorted_y)
    plt.xlim([0,6000])
    plt.ylim([0,500])
    plt.title(plot_n)
    plt.show()

def filter_muts_by_pos(list_mut_tuples, min_p = 2000, max_p=4500, complement = False):
    # filter out mutations that are withing a certain region only
    if complement == False:
        return tuple_muts_to_mutstr([t for t in list_mut_tuples if (t[1]>=min_p) and (t[1]<=max_p)])
    else:
        return tuple_muts_to_mutstr([t for t in list_mut_tuples if (t[1]<min_p) or (t[1]>max_p)])


def read_fastq(filename):
    """Reads FASTQ file and remove the special characters!"""
    n_to_sequences = {}
    n_to_qualities = {}
    with open(filename) as fh:
        while True:
            n= fh.readline().rstrip()
            seq = fh.readline().rstrip() # read base sequence
            fh.readline() # skip placeholder line
            qual = fh.readline().rstrip() #base quality line
            if len(seq) == 0:
                break
            n_to_sequences[n] = seq
            n_to_qualities[n]= qual
    return n_to_sequences, n_to_qualities

def get_min_fastq_scores_reads(
        fin,
        fastq_score_cutoff = 30,
        dout = None
    ):
    #expects a fastq file and returns dictionaries of sequence name to sequence, the list of scores, and to the min score value

    f_name = fin.split('''/''')[-1]
    # read fastq file
    n_to_s, n_to_q = read_fastq(fin)

    # plot some example q value scores
    plt.figure(figsize=(15,3))
    c=0
    for n, qr in n_to_q.items():
        if c>4 and c< 9:
            plt.plot(list(range(len(qr))), list(map(int,[FASTQSCORE_TO_NUM[qb] for qb in qr])))
        c+=1
    plt.title(f'example scores for {f_name}')
    plt.show()
    if dout != None:
        plt.savefig(dout + f'{f_name[:2]}_individual_example_scores.svg', format='SVG')

    # get all the fastq vals
    n_to_num_qs = {}
    for n, qseq in n_to_q.items():
        n_to_num_qs[n] = [FASTQSCORE_TO_NUM[qb] for qb in qseq]

    # find the minimum fastq val across read
    n_to_min_q = {}
    for n,qs in n_to_num_qs.items():
        n_to_min_q[n] = min(qs)

    # plot the distribution of minima
    plt.figure()
    plt.hist(n_to_min_q.values(), bins=100)
    plt.axvline(fastq_score_cutoff, c='red')
    plt.xlabel('qscore')
    plt.ylabel('frequency')
    plt.title('min qscore over read')
    plt.show()
    if dout != None:
        plt.savefig(dout + f'{f_name[:2]}_dist_min_scores.svg', format='SVG')

    return n_to_s, n_to_num_qs, n_to_min_q

def fetch_ee_from_fastq(fastq_in):
    # for a fastq file that was processed by vsearch with the fastq_eeout option, extract the expected error values from the header
    ee_dict = {}
    with open(fastq_in, 'r') as f:
        for line in f:
            if line.startswith('@'):
                read_id,ee = line.strip().split(';ee=')
                ee_dict[read_id[1:]] = float(ee)
    return ee_dict

def get_r_to_minq(fastq):
    # for a fastq file, find the min Q score for each read.
    r_to_minq = {}
    for r in SeqIO.parse(fastq, "fastq"):
        read_id = r.id
        read_minq = min(r.letter_annotations["phred_quality"])
        r_to_minq[read_id] = float(read_minq)
    return r_to_minq

def grep_reads_with_bc(
    bcoi,
    df,
    fastq_original = '/data/szchen/dms/snakemake/consensus-calling/fqs/2A_all.fastq',
    fastq_out_prefix = '/data/davidding/dms/dms_plants/02_pacbio/misc/tmp_',
    flank5prime = 'TGAACTATACAAATAA',
    flank3prime = 'CTCGAGGTTCGAGT'
    ):
    # tool for debugging:
    # given a particular barcode, grep it in the original fastq file, and write those reads to a new temporary file.
    # these fastq reads can then be inspected by manual alignment against the reference.

    # df is for example: '/data/szchen/dms/snakemake/consensus-calling/dfs/mutsBCDFs/df_2A_AC.csv'
    # flanking reads are the bases immediately next to these sequences.

    # print the expected mutation 
    print(df.loc[df['bc'] == bcoi].muts_clean.values)

    search_str = flank5prime + bcoi + flank3prime
    print(search_str)

    cmd = f"""grep -B 1 -A 2 --no-group-separator {search_str} {fastq_original} > {fastq_out_prefix}{bcoi}.fastq"""
    print(cmd)
    os.system(cmd)



###############################  shift mutant positions
# check the base of the offset wrt to the fasta file
def check_offset(seq_in, offset, expected_base):
    return expected_base == seq_in[offset]

def shift_filter_mutants(mut, offset, min_orig_pos, max_orig_pos):
    # filter a list of muts being within the promoter/5'UTR, ie. for the ones within min and max orig pos
    # offset the mutants so that the 5'UTR is defined as 0
    list_muts = get_ind_muts(mut) # returns list of tuples of (wt, pos, mut)

    # filter for mut tuples that are within range
    list_muts_OI = [m for m in list_muts if m[1]>=min_orig_pos and m[1]<=max_orig_pos]

    # shift the positions of the mutants
    list_muts_pos_updated = [wtMut +str(int(pos) - offset)+ mutMut for wtMut, pos, mutMut in list_muts_OI]
    print(list_muts_pos_updated)
    # if there are no mutations found in the promoter region, then it's a wild type read
    if len(list_muts_pos_updated) == 0:
        return 'wt'
    else:
        return ':'.join(list_muts_pos_updated)

############################### get the forbidden sequences.



'''
3B: 
    A23:
    Present in many barcodes
    T28:
    present in 32866 out of 34282 barcodes
    30, 31
    Probably questionable, close to 23 and 28
    G52
    Too many A’s
    A136_
    Repeating A’s
    583
    Repeating A’s
    634
    Repeating A’s
    A2027G:T2028A:C2029T':
    repeating C's, in almost all barcodes
    2030
    Repeating C’s

3C:
    614
    Repeating gt’s
    639, 638, 637
    Flanked by repeating gt’s and repeating t’s
    842
    13 repeating As
    1097
    Many t’s nearby
'''
dic_muts_ignore ={
    # '3B': [23, 28, 2027, 2028, 2029], # we see A23G, T28A in most samples, and A2027G:T2028A:C2029T in almost all
    '3B': list(range(32)) + [52, 136, 583, 634, 2027, 2028, 2029, 2030],
    '3C': [614, 637, 638, 639, 842, 1097], # homo region or just too common
}


def filter_muts_to_ignore(mut, list_pos_ignore):
    # filter a mutstring like 'A23G:T28A:G581A:G1429_________________________________________________________________________________________________________________________:A2027G:T2028A:C2029T'
    # to remove the positions in the list_pos_ignore to: G581A:G1429_________________________________________________________________________________________________________________________
    # input 'A23G:T28A:A2027G:T2028A:C2029T' -> 'wt'
    if mut == 'wt':
        return 'wt'
    tuple_muts = get_ind_muts(mut)
    muts_keep = [m for m in tuple_muts if int(m[1]) not in list_pos_ignore]
    if len(muts_keep) == 0:
        return 'wt'
    else:
        return ':'.join([''.join(map(str,m)) for m in muts_keep])
