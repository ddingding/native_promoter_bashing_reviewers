import re
import numpy as np
from itertools import groupby
from collections import Counter
import constants as co

def hamming(str1, str2):
    assert len(str1) == len(str2)
    return sum(c1 != c2 for c1, c2 in zip(str1, str2))

def complement(inbase):
    cDic = {"A": "T", "T": "A", "C": "G", "G": "C"}
    return cDic[inbase]

def rev_complement(instr):
    compl = [complement(c) for c in instr]
    return "".join(compl[::-1])


def fasta_iter(fasta_name):
    """
	given a fasta file,  yield tuples of header, sequence
	https://github.com/GuyAllard/fasta_iterator/blob/master/fasta_iterator/__init__.py
    """
    rec = None
    for line in open(fasta_name, "r"):
        if line[0] == ">":
            if rec:
                yield rec
            rec = FastaRecord(line.strip()[1:])
        else:
            rec.sequence += line.strip()

    if rec:
        yield rec

class FastaRecord:
    """
    Fasta Record, contains a header and a sequence
	# from https://github.com/GuyAllard/fasta_iterator/blob/master/fasta_iterator/__init__.py
    """

    def __init__(self, header="", sequence=""):
        self.header = header
        self.sequence = sequence

    def __str__(self):
        return ">{}\n{}".format(self.header, self.sequence)

def mutstr_to_dic(mutStr):
    # take mutstr like '4:T,5:A,8:A' and make a dictionary of position to mutant nt
    return dict(
        zip(
            [int(x.split(":")[0]) for x in mutStr.split(",")],
            [x.split(":")[1] for x in mutStr.split(",")],
        )
    )

def mutdic_to_mutcodon(mutpos_to_base, wtCodon):
    # take  {'100': 'A', '101': 'G'} and convert to mutcodon, which is a 3 letter base string.
    mutCodon = ""
    codonPos = min(mutpos_to_base.keys()) // 3 * 3

    for i in range(3):
        if i in [int(x) % 3 for x in mutpos_to_base]:
            mutCodon += mutpos_to_base[codonPos + i]
        else:
            mutCodon += wtCodon[i]
    return mutCodon

def break_mut(m):
    # breaks a mutation into wt, pos, mut, handling minus signs.
    # nt.break_mut('T-3575AAA') -> ('T', -3575, 'AAA')

    try:
        wt = m[0]
    except IndexError:
        print('indexerror taking the first element of m for: ',m, ind_muts)
        return []
    try:
        pos = int(re.findall(r'-?\d+', m)[0]) # can handle minus signs.
    except IndexError:
        print('indexerror in finding position for: ',m)
        return []
    mut = m.split(str(pos))[1]

    return (wt, pos, mut)

def get_ind_muts(muts):
    # pass it a string so it can handle nan
    # splits a mutstr like C247_:T324_:T384_:T409A
    # into a list of tuples indicating wt, pos and mut

    if muts.lower() == 'nan':
        return []
    
    if muts == 'wt':
        return ['wt']

    try:
        ind_muts = muts.split(':')
    except AttributeError:
        #print('attributeError:',muts)
        return []

    list_muts_sep = []
    for m in ind_muts:
        (wt, pos, mut) = break_mut(m)

        list_muts_sep.append((wt, pos, mut))
    return list_muts_sep

def shift_muts(muts, shift_pos):
    # shift all the positions in a mutstr by a given amount
    # nt.shift_mut('C247_:T324_:T384_:T409A', 100) -> 'C347_:T424_:T484_:T509A'
    if muts == 'wt':
        return 'wt'
    else:
        return ':'.join([f'{wt}{pos+shift_pos}{mut}' for wt, pos, mut in get_ind_muts(muts)])

def is_in_promoter_utr(muts, start_promoter, end_utr):
    # check whether all mutations are in promoter region defined by position fo start_promoter and end_utr
    if str(muts).lower() == 'nan':
        return 'nan'
    if muts == 'wt':
        return True
    for wt, pos, mut in get_ind_muts(muts):
        if pos < start_promoter or pos > end_utr:
            return False
    return True

def add_tss0_mut_col(df,gene_n, mut_col ='muts_oprc'):
    # shifting a column to be 0 wrt to the TSS

    # subtract one to get the 0-based position of the TSS
    tss_pos = co.sample_to_positions[gene_n][1] -1
    len_5utr = co.sample_to_positions[gene_n][-1] - tss_pos +1
    #print('length of 5UTR:', len_5utr)

    df['mut_tss0'] = df[mut_col].apply(lambda m: shift_muts(m, -tss_pos))
    return df

def count_ind_muts(df, mut_col = 'muts_wo_bc'):
    # for a dataframe, get a counter of all mutations in a given column
    list_ind_muts = []
    for mutstr in df[mut_col].dropna().values:
        for m in mutstr.split(':'):
            list_ind_muts.append(m)

    counter_ind_muts = Counter(list_ind_muts)
    print('found ',len(counter_ind_muts), ' individual mutations.')
    return counter_ind_muts

def count_ind_mut_pos(df, mut_col = 'muts_wo_bc'):
    # for a dataframe, get a counter of all mutation positions in a given column
    list_ind_mut_pos = []
    for mutstr in df[mut_col].dropna().values:
        for wt, m_pos, mut in get_ind_muts(mutstr):
            list_ind_mut_pos.append(m_pos)

    counter_ind_mut_pos = Counter(list_ind_mut_pos)
    print(len(counter_ind_mut_pos))
    return counter_ind_mut_pos

def clean_muts(str_muts):
    # get rid of any insertions at position 0 (wrong reference or pacbio remnants)
    # collate all the consecutive deletions (like A2917_, A2918_) into a single deletion
    # TODO: inspect for repetitive base insertion elements (like an extra A in an AAAAA)
    if str(str_muts).lower() == 'nan':
        return np.nan
    if str(str_muts).lower() == 'wt':
        return 'wt'


    agg_pos_to_mut_tuple = {}

    list_mut_tuples = get_ind_muts(str_muts)
    pos_to_mut_tuple = dict(zip([p for _,p,_ in list_mut_tuples],list_mut_tuples))

    ####### aggregate consecutive deletions into a single deletion

    # make a dictionary of position to mut_tuple for deletions
    pos_to_mut_tuple_del = dict(zip([p for p,v in pos_to_mut_tuple.items() if v[2].endswith('_')],
                                    [v for p,v in pos_to_mut_tuple.items() if v[2].endswith('_')]
                                    )
                                )
    del_pos = list(pos_to_mut_tuple_del.keys())

    # Group by differences (consecutive integers have equal differences)  
    gb = groupby(enumerate(del_pos), key=lambda x: x[0] - x[1]) # x[0] is the index in list, and x[1] is the value

    # Repack elements from each group into list
    all_groups = ([i[1] for i in g] for _, g in gb) # ie. [[247], [324], [384], [580], [602], [927], [1009], [1141], [1230], [1369], [1692, 1693]]

    # write to aggregation dictionary
    for l_pos_del in all_groups:
        p = min(l_pos_del)
        mut_tuple_del = pos_to_mut_tuple_del[p]
        agg_pos_to_mut_tuple[p] = (mut_tuple_del[0],str(p),'_'*len(l_pos_del))

    ###### add remaining mutations without insertions at the beginning    

    # add all other mutations
    for p, t in pos_to_mut_tuple.items():
        if p not in del_pos:
            if p != 0: # get rid of pos 0 insertions that might be 
                agg_pos_to_mut_tuple[p] = t

    ####### write out

    #sort it by position and return
    new_list_muts = []
    for p in sorted(agg_pos_to_mut_tuple.keys()):
        new_list_muts.append(''.join(map(str,list(agg_pos_to_mut_tuple[p]))))

    if new_list_muts != []:
        return ':'.join(new_list_muts)
    else:
        return 'wt'

def fetch_high_reads(df, min_reads = 20, min_reads_rna = 0):
    return df.loc[
        (df.bc_count_dna_rep1 > min_reads) &((df.bc_count_rna_rep1 > min_reads_rna) & (df.bc_count_rna_rep2 > min_reads_rna))]

def get_df_with_muts(df, mut, mut_col = 'muts_wo_bc', strict=True):
    # select a subset of the rows that contain a particular mutation exactly 
    ##get_df_with_muts(df, 'T5442G') should return 148 rows for 3B
    # if using 'strict', then only report rows that contain that exact mutations
    # if not using 'strict', then report rows that contain that mutation anywhere in the string, 
    # ie. looking for deletions of a minimal size.

    df_nona = df.dropna(subset=[mut_col])

    if strict:
        df_nona['contains_mut'] = df_nona[mut_col].apply(lambda x: mut in x.split(':'))
        df_mut = df_nona.loc[df_nona.contains_mut]
    else:
        df_mut = df_nona.loc[df_nona[mut_col].str.contains(mut)]
    return df_mut

def get_df_wt_mut_del(
    df, 
    mut_col = 'muts_clean',
    mut_str = 'AGCTGGCTTGTGGGGACCAGACAAAAAAGGAATG',
    mut_strict = True,
    n_del = 300 ):
    # prepare data
    df_wt = df.loc[df[mut_col] == 'wt'] # fetching the wildtype reads
    print('found wt reads:', len(df_wt))
    df_mut = get_df_with_muts(df,mut_str, mut_col=mut_col, strict=mut_strict) # assuming this is the beginning of the enhancer.
    print(f'found mut reads:', len(df_mut))
    df_del = get_df_with_muts(df,'_'*n_del, mut_col=mut_col, strict=False) # need to use column where consecutive deletions are merged into a string
    print('found deletion reads:', len(df_del))
    return df_wt, df_mut, df_del