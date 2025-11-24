# functions for generating the promoter libraries

from constants import bases, re_to_site
from nucleotides import hamming, rev_complement
from collections import Counter
import numpy as np
import primer3 as p3

def gen_singles(seq, seq_index):
    # receive a segment and return a dictionary of name to all single bp segments
    len_wt_seq = len(seq)
    dict_singles = {}
    for i in range(len(seq)):
        curr_base = seq[i]
        for b in [b for b in bases if b!= curr_base]:
            seq_single = seq[:i] + b + seq[i+1:]
            n= 'single_'+curr_base + str(seq_index+i) + b
            #make sure all singles have the right length
            assert len(seq_single) == len_wt_seq
            dict_singles[n] = seq_single
    return dict_singles

def gen_dels(seq, seq_index, del_len):
    # receive a segment (+ del_len at the 5prime end), and return dictionary of all del
    # name of segment corresponds to end of deletion
    dict_del = {}
    for i in range(len(seq)-del_len+1):
        n= str(del_len)+'bp_del_'+ str(seq_index + i + del_len)
        seq_del = seq[:i] + seq[i+del_len:]
        dict_del[n] = seq_del
    return dict_del


def gen_abe_mimics(seq, seq_index, win_len):
    # takes a sequence and converts all A-->G within a certain window
    # as well as T--> C within a certain window (opposite strand)
    #name corresponds to editing window end.
    
    abe_dict = {}
    for i in range(win_len, len(seq)):
        # get editing window.
        seq_win = seq[i:i+win_len] 
        
        # if only one A in it, skip
        if Counter(seq_win)['A'] >1:
            # replace all A with Gs.
            edited_seq_win = seq_win.replace('A', 'G')
            
            n= 'ABE_AtoG_win'+str(win_len)+'_p'+str(seq_index + i +win_len)
            seq_abe_ag = seq[:i] + edited_seq_win + seq[i+win_len:]
            abe_dict[n] = seq_abe_ag
            
        if Counter(seq_win)['T'] > 1:
            # replace all T with Cs:
            edited_seq_win = seq_win.replace('T', 'C')
            
            n= 'ABE_TtoC_win'+str(win_len)+'_p'+str(seq_index + i+ win_len)
            seq_abe_ag = seq[:i] +  edited_seq_win + seq[i+win_len:]
            abe_dict[n] = seq_abe_ag
            
    return abe_dict

def gen_random_inserts(seq, start_pos, end_pos, len_rand, n_seqs):
    #takes a sequence, chooses a random position between start and stop of that sequence
    # and inserts a random sequence of length len_rand
    dict_rand={}
    for c in range(n_seqs):
        pos_insert = start_pos + n_seqs% (end_pos - start_pos)
        # choose random positions
        #pos_insert = np.random.choice(list(range(start_pos, end_pos)))
        # generate a random seq
        rand_seq = ''.join([np.random.choice(bases) for _ in range(len_rand)])
        
        n = 'rand_p'+str(pos_insert)+'_'+rand_seq
        
        full_rand_seq = seq[:pos_insert] + rand_seq + seq[pos_insert:]
        dict_rand[n] = full_rand_seq
    return dict_rand


    
def gen_motif_inserts(seq, start_pos, end_pos, dic_motifs, stagger = 5):
    # for a sequence isnerts all motifs  in dic_motifs (name to motifsequence) every 'stagger' position between start and stop pos
    dic_n_to_seq = {}
    for i in range(start_pos, end_pos, stagger):
        pre_seq = seq[:i]
        suff_seq = seq[i:]
        for im, m in dic_motifs.items():
            n= 'ins_motif_m{}_p{}_fwd'.format(im, i)
            ins_seq = pre_seq + m + suff_seq
            dic_n_to_seq[n] = ins_seq
            
            n= 'ins_motif_m{}_p{}_rev'.format(im, i)
            ins_seq = pre_seq + rev_complement(m) + suff_seq
            dic_n_to_seq[n] = ins_seq
    
    return dic_n_to_seq

# write small list
def write_fasta_df_muts(df_muts, fout, c_write = 'seq_order', fasta_name = 'Unnamed: 0' ):
    with open(fout, 'w') as f:
        sorted_idxs = sorted(df_muts[fasta_name])
        c=0
        for idx in sorted_idxs:
            #if c<100:
            f.write('>'+ idx + '\n')
            f.write(df_muts.loc[df_muts[fasta_name]==idx][c_write]+'\n')   
            c+=1

'''
# testing
d = gen_random_inserts('AAAAAAAAAAAAAAAAAAAAA', -5, 0, 10, 10)
d

df_rand = pd.DataFrame.from_dict(d, orient='index', columns= [ 'seq'])
df_rand

dic_abe = gen_abe_mimics('AATAAAATCATT', -1000, 5)
df_abe = pd.DataFrame.from_dict(dic_abe, orient='index', columns= [ 'seq'])
pd.concat([df_rand, df_abe], axis=0)

#test
print(len('mydogisacat'))
gen_dels('mydogisacat', -1000, 2)

# tst motif insertions
gen_motif_inserts('AAAAAAAAAAAAAAAAAAAAAAAAA', 5, 11, n_to_m_ins)
'''

def get_gene_library(prime_5, utr_5, exon_1, c, gene_name, n_to_ms= None):
    '''
    INPUT:
    - take a 5' untranscribed region, 5'UTR and exon_1, 
    - count of segment (specifies primer pair to choose)
    - name of gene
    
    and make 
    - all singles in 1 kb within 1kb
    - all 2 bp deletions
    - all 10 bp deletions
    - all A's in a particular window
    
    - prime edits.

    RETURN: 
    - df with all mutations
    - count for primer pairs
    '''
    
    prime_utr_5 = prime_5 + utr_5

    full_5_prime = ''

    df_1kb_all = pd.DataFrame()

    # iterate through all 200 bp segments in one promoter including 5'UTR:
    for i in range((-1000-len(utr_5)), 0, 200):
        start = i

        ############## get the segment to mutate,seg_curr,and the bsmbI restriction overhangs
        
        # adding bmsbI restriction sites
        # generate Bsa-HF-v2 restriction sites.
        # https://www.neb.com/products/r3733-bsai-hf-v2#Product%20Information
        # pre RE is invariant to where 3' end falls
        bp6_pre_seg = prime_utr_5[start-6: start]
        
        
        # variable to add sequence at the 5'end to prevent mutated region being too small, might fall through PCR purification after RE digest ()
        prime_5_add = ''
        # to catch the last segment that goes to translational start and is shorter than 200 bp
        # 0 denotes the end of the gene, or translational start
        if start+200 >= 0:
            print('doing last segment')
            # get the whole end of the 5' utr
            seg_curr = prime_utr_5[start:]
            
            # using the first 6 bases of translational start to cut and ligate if segment is shorter
            prime_5_add =prime_utr_5[start-200+len(seg_curr):start]
            bp6_pre_seg =  prime_utr_5[start-200+len(seg_curr)-6:start-200+len(seg_curr)]
            bp6_post_seg = exon_1[:6]
        
        # if segment to mutate is in 5' region and RE overhang is also in 5' region
        elif (start+200 < 0) and (start+206 < 0):
            seg_curr = prime_utr_5[start:start+200]
            bp6_post_seg = prime_utr_5[start+200:start+206]
        
        # if segment to mutate is within 5' region, but the 6bp RE overhang overlaps the start codon
        elif (start+200 <0) and (start+206 >= 0):
            bp6_post_seg_5prime = prime_utr_5[start+200:]
            coding_segment =  exon_1[:(6-len(bp6_post_seg_5prime))]
            bp6_post_seg =  bp6_post_seg_5prime+ coding_segment
            print('created a mixed overhang for RE enzyme for gene {} at segment pos{},from 5:{} and exon seq:{} '.format(gene_name, start, bp6_post_seg_5prime,coding_segment))
            

            
        print('current segment to mutate:', seg_curr)
        print('len of current segment to mtuate', len(seg_curr))

        # reconstruct full 5 prime and utr and compare for sanity at end
        full_5_prime += seg_curr


        ###### generate mutations
        # generate all singles
        singles_dic = gen_singles(seg_curr, i)
        df_singles = pd.DataFrame.from_dict(singles_dic, orient='index', columns= [ 'min_chunk'])

        print('gene {}, segment {}, singles #: {}'.format(gene_name,start, len(df_singles)))
        
        # generate all 2 bp deletions
        del_2bp_dic = gen_dels(seg_curr,i, del_len=2)
        df_2bp_del = pd.DataFrame.from_dict(del_2bp_dic, orient='index', columns= [ 'min_chunk'])
        print('gene {}, segment{}, 2bp dels #: {}'.format(gene_name, start,len(df_2bp_del)))

        del_12bp_dic = gen_dels(seg_curr,i, del_len=12)
        df_12bp_del = pd.DataFrame.from_dict(del_12bp_dic, orient='index', columns= [ 'min_chunk'])
        print('gene {}, segment {}, 12bp dels #: {}'.format(gene_name,start, len(df_12bp_del)))

        # generate abe mimics
        win_len = 5
        abe_dic = gen_abe_mimics(seg_curr, i, win_len)
        df_abe = pd.DataFrame.from_dict(abe_dic, orient='index', columns= [ 'min_chunk'])

        print('gene {}, segment {}, abe #: {}'.format(gene_name, start, len(df_abe)))


        ######## consolidating mutations
        # creating a dataframe with all mutants
        df_segment = pd.concat([df_singles, df_2bp_del, df_12bp_del, df_abe], axis=0)

        # adding restriction enzymes
        bsa1_pre = 'GGTCTC' + bp6_pre_seg # make sense 5' overhang
        bsa1_post = bp6_post_seg + 'GAGACC' #  makes antisense 3' overhang

        df_segment['bsa1_pre'] = bsa1_pre
        df_segment['bsa1_post'] = bsa1_post

        # add the bases at 3' end to make sure no sequences are too small
        
        df_segment['prime_5_add'] = prime_5_add
        ####### adding orthogonal primers
        p1 = df_primers_choose.iloc[c*2].Sequence
        p2 = rev_complement(df_primers_choose.iloc[c*2+1].Sequence)
        df_segment['amp_fwd'] = p1
        df_segment['amp_rev'] = p2
        df_segment['primer_ids'] =  str(df_primers_choose.iloc[c*2]['Primer id']) +'_' + str(df_primers_choose.iloc[c*2+1]['Primer id'])
        df_segment['segment_loc'] = start
        # add to overall dataframe
        df_1kb_all = pd.concat([df_1kb_all, df_segment], axis=0)


        c+=1
        print(c)

    print('number of singles, abes, 2+10bp deletions:',len(df_1kb_all))

    '''
    Todos: insert transcriptional activator motifs.
    '''
    
    '''
    # generate random isnertions in 5' promoter region only (no 5' UTR inserts)
    random_insert_dic = gen_random_inserts(prime_5[-200:], -200, 0, 10, 10000- len(df_1kb_all))
    
    df_random_insert = pd.DataFrame.from_dict(random_insert_dic, orient='index', columns= [ 'min_chunk'])
    
    # generate BsmBi-v2 restriction sites.
    bp6_pre_seg = prime_5[-206: -200]
    bp6_post_seg = utr_5[:6]

    bsa1_pre = 'GGTCTC' + bp6_pre_seg # make sense 5' overhang
    bsa1_post = bp6_post_seg + 'GAGACC' # 

    df_random_insert['bsa1_pre'] = bsa1_pre
    df_random_insert['bsa1_post'] = bsa1_post

    # add orthogonal primers
    p1 = df_primers_choose.iloc[c*2].Sequence
    p2 = rev_complement(df_primers_choose.iloc[c*2+1].Sequence)
    df_random_insert['amp_fwd'] = p1
    df_random_insert['amp_rev'] = p2
    df_random_insert['primer_ids'] =  str(df_primers_choose.iloc[c*2]['Primer id']) +'_' + str(df_primers_choose.iloc[c*2+1]['Primer id'])
    df_random_insert['segment_loc'] = -200
    '''
    if n_to_ms != None:
        
        tss_minus_start = 150
        tss_plus_end = 49
        
        motif_chunk = (prime_5+utr_5)[-len(utr_5)-tss_minus_start:-len(utr_5)+tss_plus_end]
        print('motif_chunk length for motif inserts',len(motif_chunk))
        d_m_ins = gen_motif_inserts(motif_chunk,0, len(motif_chunk) , n_to_m_ins, stagger=5)

        df_m_insert = pd.DataFrame.from_dict(d_m_ins, orient='index', columns= [ 'min_chunk'])

         # generate BsmBi-v2 restriction sites.
        bp6_pre_seg = prime_5[-tss_minus_start-6:-tss_minus_start]
        bp6_post_seg = utr_5[tss_plus_end: tss_plus_end+6]

        bsa1_pre = 'GGTCTC' + bp6_pre_seg # make sense 5' overhang
        bsa1_post = bp6_post_seg + 'GAGACC' # 

        df_m_insert['bsa1_pre'] = bsa1_pre
        df_m_insert['bsa1_post'] = bsa1_post

        # add orthogonal primers
        p1 = df_primers_choose.iloc[c*2].Sequence
        p2 = rev_complement(df_primers_choose.iloc[c*2+1].Sequence)
        df_m_insert['amp_fwd'] = p1
        df_m_insert['amp_rev'] = p2
        df_m_insert['primer_ids'] =  str(df_primers_choose.iloc[c*2]['Primer id']) +'_' + str(df_primers_choose.iloc[c*2+1]['Primer id'])
        df_m_insert['segment_loc'] = -tss_minus_start
        df_m_insert['prime_5_add'] = ''
    c+=1
    print(c)

    df_2kb = pd.concat([df_1kb_all,df_m_insert ], axis=0)
    df_2kb['seq_order'] = df_2kb.apply(lambda r: r.amp_fwd + r.bsa1_pre + r.prime_5_add + r.min_chunk +r.bsa1_post + r.amp_rev, axis=1)
    df_2kb['len_seq_order'] = df_2kb.apply(lambda r: len(r.seq_order), axis=1)
    
    # check the length of the segments adds up to the original fragment
    assert hamming(full_5_prime, prime_utr_5[-1000-len(utr_5):]) == 0

    df_2kb['gene_name'] = gene_name
    return df_2kb, c

def write_fasta_df_muts(df_muts, fout, c_write = 'seq_order' ):
    with open(fout, 'w') as f:
        sorted_idxs = sorted(df_muts.index)
        
        for idx in sorted_idxs:
            f.write('>'+ idx + '\n')
            f.write(df_muts.loc[idx][c_write]+'\n')     

def print_min_max_lens(df):
    print(min(df.len_seq_order))
    print(max(df.len_seq_order))

def filter_df_muts(df, n_motifs_keep):
    m_names_to_keep = ['m'+ str(n) for n in n_motifs_keep] # add an m for each number
    #print(len(m_names_to_keep))
    non_ins_rows = [k for k in df['Unnamed: 0'] if not (k.startswith('ins'))] 
    #print(len(non_ins_rows))
    ins_rows = [k for k in df['Unnamed: 0'] if (k.startswith('ins'))]
    #print(len(ins_rows))
    ins_rows_keep = [k for k in ins_rows if (k.split('_')[2] in m_names_to_keep)]
    #print(len(ins_rows_keep))
    rows_keep = non_ins_rows + ins_rows_keep
    #print(len(rows_keep))

    return df.loc[df['Unnamed: 0'].isin(rows_keep)]

def write_fasta_df_muts(df_muts, fout, c_write = 'seq_order', fasta_name = 'Unnamed: 0' ):
    with open(fout, 'w') as f:
        sorted_idxs = sorted(df_muts[fasta_name])
        c=0
        for idx in sorted_idxs:
            #if c<100:
            f.write('>'+ idx + '\n')
            f.write(df_muts.loc[df_muts[fasta_name]==idx][c_write].values[0]+'\n')   
            c+=1

# for each segment and gene, get the matching backbone sequences
def get_matching_backbone_seqs_library(prime_5, utr_5, first_10_aa, c, gene_name, seq_len_show = 70):
    '''
    INPUT:
    - take a 5' untranscribed region, 5'UTR and exon_1, 
    - count of segment (specifies primer pair to choose)
    - name of gene
    
    and make 
    - all singles in 1 kb within 1kb
    - all 2 bp deletions
    - all 10 bp deletions
    - all A's in a particular window
    
    - prime edits.

    RETURN: 
    - the locations where the primers should be designed
    '''
    
    dict_pnum_to_bb_seq = {}
    
    #prime_utr_5 = prime_5 + utr_5
    prime_utr_5_gene10aa =  prime_5 + utr_5 + first_10_aa


    # iterate through all 200 bp segments in one promoter including 5'UTR:
    for i in range((-1000-len(utr_5) - len(first_10_aa)), -len(first_10_aa), 200):
        start = i

        ############## get the segment to mutate,seg_curr,and the bsmbI restriction overhangs
        
        # adding bmsbI restriction sites
        # generate Bsa-HF-v2 restriction sites.
        # https://www.neb.com/products/r3733-bsai-hf-v2#Product%20Information
        # pre RE is invariant to where 3' end falls
        #bp6_pre_seg = prime_utr_5_gene10aa[start-6: start]
        backbone_seq_rev = rev_complement(prime_utr_5_gene10aa[start-seq_len_show: start]) # unless changed downstream

        
        # variable to add sequence at the 5'end to prevent mutated region being too small, 
        #might fall through PCR purification after RE digest ()

        # to catch the last segment that goes to translational start and is shorter than 200 bp
        # 0 denotes the end of the gene, or translational start
        prime_5_add = '' # to add if sequences are too small
        if start+200 >= -len(first_10_aa):
            print('doing last segment')
            # get the whole end of the 5' utr
            seg_curr = prime_utr_5_gene10aa[start:-len(first_10_aa)]
            
            
            #prime_5_add =prime_utr_5[start-200+len(seg_curr):start]

            # using the first 6 bases of translational start to cut and ligate if segment is shorter
            #prime_5_add =prime_utr_5_gene10aa[start-200+len(seg_curr):start]
            #bp6_pre_seg =  prime_utr_5_gene10aa[start-200+len(seg_curr)-6:start-200+len(seg_curr)]
            rev_seq_start = start-200+len(seg_curr) # this is to account for the fact that the 5'utr chunks were made a bit longer
            #rev_seq_start = start
            rev_seq_end = rev_seq_start-seq_len_show
            print('rev_seq_indices', rev_seq_start, rev_seq_end)
            backbone_seq_rev = rev_complement(prime_utr_5_gene10aa[rev_seq_end:rev_seq_start])
            print(backbone_seq_rev)
            #bp6_post_seg = first_10_aa[:6]
            backbone_seq_fwd = first_10_aa[:seq_len_show]
        # if segment to mutate is in 5' region and RE overhang is also in 5' region
        elif (start+200 <  -len(first_10_aa)) and (start+206 <  -len(first_10_aa)):
            #seg_curr = prime_utr_5_gene10aa[start:start+200]
            #bp6_post_seg = prime_utr_5_gene10aa[start+200:start+206]
            backbone_seq_fwd = prime_utr_5_gene10aa[start+200:start+200 + seq_len_show]
            
        # if segment to mutate is within 5' region, but the 6bp RE overhang overlaps the start codon
        elif (start+200 < -len(first_10_aa)) and (start+206 >=  -len(first_10_aa)):
            #bp6_post_seg_5prime = prime_utr_5_gene10aa[start+200:]
            #coding_segment =  first_10_aa[:(6-len(bp6_post_seg_5prime))]
            #bp6_post_seg =  bp6_post_seg_5prime+ coding_segment
            backbone_seq_fwd = prime_utr_5_gene10aa[start+200:] + first_10_aa[:seq_len_show]
            print('created a mixed overhang for RE enzyme for gene {} at segment pos{},from 5:{} and exon seq:{} '.format(gene_name, start, bp6_post_seg_5prime,coding_segment))
        

        # get the primer sequence associated with fragment    
        p1_num = df_primers_choose.iloc[c*2]['Primer id']
        p2_num = df_primers_choose.iloc[c*2+1]['Primer id']
        #
        print(gene_name + ' curr_primers '+ str(p1_num) + ' '+ str(p2_num))
        print(start, prime_utr_5_gene10aa[start:start+10])
        print(start, prime_utr_5_gene10aa[start+200: start+210])

        # make a dictionary to record backbone primer sequences
        dict_pnum_to_bb_seq[str(p1_num)+'_'+str(p2_num)+'_rev'] = backbone_seq_rev
        dict_pnum_to_bb_seq[str(p1_num)+'_'+str(p2_num)+'_fwd'] = backbone_seq_fwd

        #todo
            

        c+=1
        print(c)

    #generate backbones for the TF insertion chunks
    
    tss_minus_start = 150 + len(first_10_aa) + len(utr_5)
    tss_plus_end = 49 - len(first_10_aa) - len(utr_5)

     # generate BsmBi-v2 restriction sites.
    bp6_pre_seg = prime_utr_5_gene10aa[-tss_minus_start-6:-tss_minus_start]
    bp6_post_seg = prime_utr_5_gene10aa[tss_plus_end: tss_plus_end+6]

    backbone_seq_rev = rev_complement(prime_utr_5_gene10aa[-tss_minus_start-seq_len_show:-tss_minus_start])
    backbone_seq_fwd= prime_utr_5_gene10aa[tss_plus_end: tss_plus_end+seq_len_show]
    
    # get the primer numbers
    # get the primer sequence associated with fragment
    p1_num = df_primers_choose.iloc[c*2]['Primer id']
    p2_num = df_primers_choose.iloc[c*2+1]['Primer id']

    # make a dictionary to record backbone primer sequences
    dict_pnum_to_bb_seq[str(p1_num)+'_'+str(p2_num)+'_rev'] = backbone_seq_rev
    dict_pnum_to_bb_seq[str(p1_num)+'_'+str(p2_num)+'_fwd'] = backbone_seq_fwd
    #to do: add gatacaGGTCTC to backbone_seq_rev
    # todo: add gatacaGGTCTC to backbone_seq_fwd
    c+=1
    

    return dict_pnum_to_bb_seq, c



def get_deletions_primers(n_gene, fl_seq, prime_5, utr_5, exon_1, chunk_size=200, max_len_primer = 90, min_tm = 60):

    prime5_exon1_full = prime_5 + utr_5 + exon_1
    
    # check prime5_full and first 30 bases for any restriction sites --> choose the resitrction enzyme based on that.
    
    # TODO: need to check the entire gene
    bsa1_site = 'GGTCTC' # cut N1 and 3'overhang N5
    bsmb1_site = 'CGTCTC' # cut N1 and 3'overhang N5
    aar1_site = 'CACCTGC' # cut N4 and 3'overhang N8 https://www.neb.com/products/r0745-paqci#Product%20Information
    re_enzyme_poss =[]
    for re_name, re_site in [('bsa1',bsa1_site), ('bsmb1', bsmb1_site), ('aar1', aar1_site)]:
        if (re_site not in fl_seq) and (rev_complement(re_site) not in fl_seq):
            re_enzyme_poss.append((re_name, re_site))
    print(n_gene, re_enzyme_poss)

    re_enzyme_to_use = re_enzyme_poss[0]
    re_name_to_use = re_enzyme_to_use[0]
    re_site_to_use = re_enzyme_to_use[1]
    print(len(exon_1))
    dic_primers = {}
    for del_start in range(-2000-len(exon_1),-len(exon_1),chunk_size):
        rev_prim = ''
        for l_rev in range(30,max_len_primer-12-6):
            curr_rev_seq = prime5_exon1_full[del_start-l_rev:del_start]
            if p3.calcTm(curr_rev_seq)>min_tm:
                rev_prim= rev_complement(curr_rev_seq)
                #print('found rev backbone sequece')
                break

        # add the first 6bp to cut
        if re_name_to_use.startswith('b'):
            full_rev_prim = 'GATACA' + re_site_to_use.lower() +rev_complement(prime5_exon1_full[del_start+chunk_size:del_start+chunk_size+6]) + rev_prim.lower()

        # if using aarI, need to use 8 bp from the deletion site at the 3' end (of the fwd primer) to make the right cut size
        elif re_name_to_use == 'aar1':
            full_rev_prim = 'GATACA'+ re_site_to_use.lower() + rev_complement(prime5_exon1_full[del_start+chunk_size:del_start+chunk_size+8]) + rev_prim.lower()
        if rev_prim =='':
            print('no rev primer found')
        
        
        if re_name_to_use.startswith('b'):
            fwd_del_start = del_start
        
        # if using aarI, need to go back 4 bases, to make the right 4bp 3' overhang for the rev primer.
        elif re_name_to_use == 'aar1':
            fwd_del_start = del_start -4
                
        fwd_prim = ''
        
        for l_fwd in range(30,max_len_primer-12):
            curr_fwd_seq = prime5_exon1_full[fwd_del_start+chunk_size:fwd_del_start+chunk_size+l_fwd]
            if p3.calcTm(curr_fwd_seq)>min_tm:
                fwd_prim= curr_fwd_seq
                #print('found fwd backbone sequece')
                break    
        # this should cover both aarI and bsmbI/bsaI sites.
        full_fwd_prim = 'gataca' + re_site_to_use + fwd_prim.lower()
        
        if fwd_prim =='':
            print('no rev primer found')

        dic_primers['g_{}_p_{}_dsize_{}_re_{}_{}'.format(n_gene, del_start, chunk_size, re_name_to_use, 'f')] = full_fwd_prim
        dic_primers['g_{}_p_{}_dsize_{}_re_{}_{}'.format(n_gene, del_start, chunk_size, re_name_to_use, 'r')] = full_rev_prim

        #print(fwd_prim)
        #print(rev_prim)

        #print(full_fwd_prim, len(full_fwd_prim))
        #print(full_rev_prim, len(full_rev_prim))

        #print(del_start)
        
    return dic_primers

################################### transcription factor motif functions ####################


def get_motif_from_meme_file(f):
    # takes as inptu a filepath to a .meme file downloaded form plantreg db
    # returns the highest motif

    bases = 'ACGT'

    bases_to_get=False
    motif = ''
    for l in open(f, 'r'):
        # check that alphabet hasn't changed
        if l.startswith('ALPHABET'):
            assert l == 'ALPHABET= ACGT\n'

        # check whether base frequency lines have already ended
        if l.startswith('\n'):
            bases_to_get=False    
        # if the base frequencly line is true: add the base
        if bases_to_get == True:
            split_freqs = l.rstrip('\n').rstrip('\t').split('\t')

            # remove whitespace from left
            split_freqs_no_space = [k.lstrip(' ') for k in split_freqs]
            #print(split_freqs_no_space)
            freq_list = list(map(float, split_freqs_no_space))

            # to check we got all the 4 frequences
            assert len(freq_list) == 4

            index_top_base = freq_list.index(max(freq_list))
            top_base = bases[index_top_base]
            motif += top_base

        # if we reach lines to start reading the base:
        if l.startswith('letter-probability matrix:'):
            bases_to_get=True

    return(motif)



def reduce_motifs(db_ms):
    # take a list of motifs and return a list of motifs that 
    # do not have any smaller motifs that are found in other motifs in the list
    # have to convert to set so that duplicates don't interfere with each other
    
    ms_use = set(db_ms)
    small_list = []
    for k in db_ms:
        found_in_other = sum([(k in m) for m in ms_use])
        num_dups = sum([(k == m) for m in ms_use])
        if found_in_other >1 and found_in_other > num_dups:
            continue
        else:
            small_list.append(k)
            
    print(len(ms_use))
    print(len(set(small_list)))
    
    return set(small_list)



def get_primers_df(df_filtered, gene_n):
    primer_pairs = sorted(
        [(int(k.split('_')[0]), int(k.split('_')[1])) for k in set(df_filtered.primer_ids)]
    )

    df_primers = pd.DataFrame(primer_pairs, columns=['p1', 'p2'])
    df_primers['p1_seq'] = df_primers.apply(
        lambda r: dic_pc_num_to_seq[r.p1], axis=1
                                        )
    df_primers['p2_seq'] = df_primers.apply(
        lambda r: dic_pc_num_to_seq[r.p2], axis=1
                                        )

    df_primers['gene'] = gene_n
    return df_primers


def add_t2s_sites(d_bb_seq, re_to_use):
    d_keep = {}
    for n, s in d_bb_seq.items():
        d_keep[n] = 'GATACA'+re_to_site[re_to_use] + s
    return d_keep

def get_tm_bb_primers(d_bb_seqs, min_tm = 63):
    d_to_keep = {}
    
    for n, seq in d_bb_seqs.items():
        for i in range(len(seq)+1):
            curr_seq = seq[:i]
            if p3.calcTm(curr_seq)>min_tm:
                d_to_keep[n] = curr_seq
                print(n, p3.calcTm(curr_seq))
                break
            if i == len(seq):
                print('hit total length')
                d_to_keep[n] = curr_seq
                print(n, p3.calcTm(curr_seq))
    return d_to_keep