from phylogeny_utilities.utilities import *
import numpy as np
import tipp_dev_utils.ext.nwops as nwo
import datetime, sys, multiprocessing, subprocess as sp
from collections import Counter
import os

nproc = 20
work = '/projects/tallis/nute/data/tippref_2018'
tippref_fold='/projects/tallis/nute/data/tippref_2018'

cogs = ['COG0012','COG0016','COG0018','COG0048','COG0049','COG0052',
        'COG0080','COG0081','COG0087','COG0088','COG0090','COG0085', #<--removed for now
        'COG0091','COG0092','COG0093','COG0094','COG0096','COG0097',
        'COG0098','COG0099','COG0100','COG0102','COG0103','COG0124',
        'COG0172','COG0184','COG0185','COG0186','COG0197','COG0200',
        'COG0201','COG0202','COG0215','COG0256','COG0495','COG0522',
        'COG0525','COG0533','COG0541','COG0552']

nucleo_alnmt_lens={'COG0012':5997,'COG0016':18762,'COG0018':33552,'COG0048':6264,
                    'COG0049':5598,'COG0052':29463,'COG0080':6372,'COG0081':6333,
                    'COG0085':101391,'COG0087':10275,'COG0088':10410,'COG0090':3816,
                    'COG0091':12645,'COG0092':17901,'COG0093':2868,'COG0094':7263,
                    'COG0096':2232,'COG0097':4476,'COG0098':13278,'COG0099':2247,
                    'COG0100':4497,'COG0102':5511,'COG0103':8007,'COG0124':33450,
                    'COG0172':18024,'COG0184':5961,'COG0185':3414,'COG0186':7644,
                    'COG0197':4113,'COG0200':10707,'COG0201':13467,'COG0202':10086,
                    'COG0215':41127,'COG0256':4356,'COG0495':68034,'COG0522':4950,
                    'COG0525':75000,'COG0533':40638,'COG0541':21147,'COG0552':115590}

# cog='COG0552'
# print("Running %s" % cog)

def td2flt(td):
    '''
    python time delta object converted to seconds as a floating point.
    :param td:
    :return:
    '''
    return td.seconds+td.microseconds/1000000

def nuc2aa(s):
    '''
    Converts a string of DNA into a string of proteins by looking the amll up.
    '''
    return ''.join(map(lambda x: cv.codons_nostop.get(x,'X'), map(''.join, zip(*[iter(s.upper())]*3))))

def nuc2codons(s):
    return list(map(''.join, zip(*[iter(s.upper())]*3)))

def str_to_position_map(long_string):
    '''
    This takes a string from a gappy alignment and returns a list of tuples where
    each each non-gapped character is represented by a tuple of the form
    (position, character).
    '''
    inds = np.where(np.frombuffer(long_string.encode('utf-8'), dtype=np.uint8) != 45)[0]
    return list(zip(inds, [long_string[j] for j in inds]))

# aa_raw, nu_raw, aa_aln, nu_aln, nm_old2new, nm_new2old = ts.get_raw_alnd_nuke_aa_fastadicts(cog)
# ks = list(aa_aln.keys())
# ks_old = list(map(lambda x: nm_new2old[ks[x]], range(len(ks))))


#
#   Ran this on 3/8 without a problem
#
# def gather_all_leftover_seqs():
#     left_f = open(os.path.join(tippref_fold, 'dna_to_aa_matching', 'leftover_seqs.txt'),'w')
#     sys.stdout = left_f
#     cogs = list(filter(lambda x: x[:3] == 'COG', map(lambda x: x[:7], os.listdir('refpkg'))))
#     for cog in cogs:
#         aa_raw, nu_raw, aa_aln, nm_old2new, nm_new2old = ts.get_raw_alnd_nuke_aa_fastadicts(cog)
#         ks = list(aa_aln.keys())
#         ks_old = list(map(lambda x: nm_new2old[ks[x]], range(len(ks))))
#         seq_scores = np.loadtxt(os.path.join('COGS_aligned_nuc','seq_status_%s.txt' % cog))
#         zero_inds = np.where(seq_scores==0)[0]
#
#         # printing this CSV thing to stdout
#         for i in zero_inds:
#             res1 = run_nw_scoresonly(aa_raw[ks_old[i]], nuc2aa(nu_raw[ks_old[i]]))
#             res2 = run_nw_scoresonly(aa_raw[ks_old[i]], nuc2aa(nu_raw[ks_old[i]][1:]))
#             res3 = run_nw_scoresonly(aa_raw[ks_old[i]], nuc2aa(nu_raw[ks_old[i]][2:]))
#             print("%s, %s, %s, " % (i, ks[i], ks_old[i]) , end = '')
#             print(''.join([" %s,", ] * 18) % (res1 + res2 + res3))
#     sys.stdout= sys.__stdout__
#     left_f.close()
#             # print("fields: (aln-len, [seq len], M_ct, D_ct, G_ct)")
#             # print("nw1:  %d,  %d,  %d,  %d,  %d,  %d," % res1, end = '')
#             # print("  %d,  %d,  %d,  %d,  %d,  %d," % res2, end = '')
#             # print("  %d,  %d,  %d,  %d,  %d,  %d," % res3)

def get_raw_alnd_nuke_aa_fastadicts(cog):
    aa_raw = read_from_fasta(os.path.join(tippref_fold,'COGS_verybest_scrap',cog + '.faa'))
    nu_raw = read_from_fasta(os.path.join(tippref_fold,'COGS_verybest_scrap', cog + '.fna'))
    aa_aln = read_from_fasta(os.path.join(tippref_fold,'COGS_aligned_aa','aln_' + cog + '.faa'))
    # nu_aln = read_from_fasta(os.path.join('COGS_aligned_aa', 'aln_' + cog + '.fna'))
    nm_old2new = get_dict_from_file(os.path.join(tippref_fold,'COGS_name_map',cog + '.txt'))
    nm_new2old = get_dict_from_file(os.path.join(tippref_fold,'COGS_name_map', cog + '.txt'),keysfirst=False)
    return aa_raw, nu_raw, aa_aln, nm_old2new, nm_new2old

def get_aligned_nuc_sequence(aarw, nurw, aaln, i, aln1, aln2, nuc_tgt_alignment_len,frameshift):
    '''This function does the heavy lifting for converting an individual
    Codon sequence to an AA sequence.
    :param aarw:
    :param nurw: must be passed in with the proper frame-shifting already performed.
    :param aaln:
    :param i:
    :param aln1:
    :param aln2:
    :param nuc_tgt_alignment_len:
    :return: (string) nuc_seq_str
    '''
    p2 = part2_mapping_from_alignment(aln1, aln2)
    p3 = get_part3_mapping(aaln)
    p23 = combine_p2_p3_mapping(p2, p3)
    num_aa = len(aarw)
    nuc_seq = np.ones(nuc_tgt_alignment_len, dtype=np.uint8) * 45
    codons = nuc2codons(nurw[frameshift:])
    if len(codons) == num_aa + 1:
        nuc_seq[-3:] = str2nparr(codons[-1])
    for (k, v) in p23.items():
        if v >= 0:
            nuc_seq[v * 3:(v * 3 + 3)] = str2nparr(codons[k])
    nuc_seq_str = nparr2str(nuc_seq)
    return nuc_seq_str

def combine_p2_p3_mapping(p2map, p3map, quiet=True):
    '''
    Combines the mappings for part2 and part2 into a single dictionary.
    :param p2map:
    :param p3map:
    :param quiet:
    :return:
    '''
    p3ks = len(list(p3map.keys()))
    p23map = {}
    p2ct=0
    for k in p2map.keys():
        if p2map[k]!= -1:
            p23map[k] = p3map[p2map[k]]
            p2ct+=1
        else:
            p23map[k] = -1
    if not quiet:
        print("p23 has %s outputs." % len(list(filter(lambda x: x>=0, p23map.values()))))
        print("p2map had %s inputs and p3map had %s outputs" % (p2ct, p3ks))
    return p23map

def get_part3_mapping(aa_aligned):
    '''
    Part 3 mapping maps from the raw AA sequence positions to the aligned AA
    sequence positions. Returns in the form of a dictionary in the form of
    (raw position, aligned position). This is a little basic for the dictionary
    because the raw positions are 0-indexed and sequential, but not a problem.
    :param aa_aligned:
    :return:
    '''
    inds=np.where(str2nparr(aa_aligned)!=45)[0]
    return dict(zip(range(inds.shape[0]),inds))
    # return np.vstack((np.arange(inds.shape[0], dtype=np.int32), inds)).transpose()

def part2_mapping_from_alignment(putative_seq, actual_seq):
    '''
    Part 2 mapping goes from putative AA position to actual AA position (unaligned).

    Takes a set of two aligned sequences and returns a dict that maps the (raw) position
    of the first sequence to the (raw) position of the second. Characters in the second
    that are not mapped onto are not included.

    Returned order starts at zero for both sequences.
    '''
    put_np = np.frombuffer(putative_seq.encode('utf-8'), dtype=np.uint8)
    act_np = np.frombuffer(actual_seq.encode('utf-8'), dtype=np.uint8)
    put_cum = np.cumsum(np.ones(len(putative_seq), dtype=np.int32) * (put_np != 45)) * (put_np != 45)
    act_cum = np.cumsum(np.ones(len(actual_seq), dtype=np.int32) * (act_np != 45)) * (act_np != 45)
    return dict(zip((put_cum[put_cum != 0] - 1), (act_cum[put_cum != 0] - 1)))
    # return np.vstack(((put_cum[put_cum != 0] - 1), (act_cum[put_cum != 0] - 1))).transpose()

def run_nw(seq1,seq2):
    '''
    Run the needleman wunsch on this in the terminal. Aligned sequences are entries
    7 and 8 (0-indexed) in the output. Rest are various alignment stats from the C
    program.
    '''
    nw_path = '/projects/tallis/nute/code/misc_c_utilities/build/nw'
    nw_input=seq1 + '\n' + seq2 + '\n'
    c = sp.run(nw_path, input=nw_input, encoding='utf-8', stdout=sp.PIPE)
    alnmt = c.stdout.strip().split('\t')
    return alnmt
    # return (int(alnmt[0]), int(alnmt[1]), int(alnmt[2]), alnmt[3], alnmt[4])

def run_nw_frametest_pyext(seq1nuc, seq2):
    return nwo.find_best_start_pos(nuc2aa(seq1nuc), nuc2aa(seq1nuc[1:]), nuc2aa(seq1nuc[2:]), seq2)
    # return (-1, -1, res[1], res[3], len(res[5]), res[2], res[4], res[5], res[6])


def run_nw_frametest(seq1nuc,seq2,pretty=False):
    '''
    Run the needleman wunsch on this in the terminal. Aligned sequences are entries
    7 and 8 (0-indexed) in the output. Rest are various alignment stats from the C
    program.
    '''
    nw_path = '/projects/tallis/nute/code/misc_c_utilities/build/nw'
    nw_input=nuc2aa(seq1nuc)+ '\n' + nuc2aa(seq1nuc[1:]) + '\n' + nuc2aa(seq1nuc[2:]) + '\n' + seq2 + '\n'
    if pretty:
        sp.run([nw_path, '-3seq', '-p'], input=nw_input, encoding='utf-8')

        return
    else:
        proc = sp.Popen([nw_path,'-3seq'], stdin=sp.PIPE, stdout=sp.PIPE, encoding='utf-8' )
        c, e = proc.communicate(input=nw_input)
        # c = sp.run([nw_path,'-3seq'], input=nw_input, encoding='utf-8', stdout=sp.PIPE)
        # (s1len, s4len, match_ct, diff_ct, aln_length, n_gaps, score, aln1, aln4, fr
        alnmt = c.strip().split('\t')
        return alnmt

def run_nw_scoresonly(seq1, seq2):
    nw_path = '/projects/tallis/nute/code/misc_c_utilities/build/nw'
    nw_input = seq1 + '\n' + seq2 + '\n'
    c = sp.run(nw_path, input=nw_input, encoding='utf-8', stdout=sp.PIPE)
    a = c.stdout.strip().split('\t')
    # scores_line="(aln-len, seq1, seq2, match_ct, diff_ct, gap_ct)"
    return (int(a[4]), int(a[0]), int(a[1]), int(a[2]), int(a[3]), int(a[5]) )


def get_seq_status_and_alignment(aa_rw, nu_rw):
    '''
    NB: This function is old and was replaced by the one called
    "..._scorebased" below.

    Takes the raw nucleotides and the raw AA sequences and determines their
    compatability. Starts by computing the putative AA sequence from the raw
    nucleotide sequence, then comparing it to the raw AA sequence and grouping
    the result into a few categories, which is returned. Possible
    outcomes are:
        0 - (default)
        1 - raw/putative match exactly.
        2 - raw/putative match by and NW alignment with 5 or fewer gaps.
        3 - if not 1 or 2, frame-shift the nucleotide sequence by 1 character
            and recompute the putative AA sequence. If the result matches with
            5 or fewer gaps, status=3
        4 - if not 1,2 or 3, frame-shirt by 2 characters and recompute putative.
            if result matches with <=5 gaps, status=4.

    The returned value is a tuple (seq1, seq2, status) where seq1 and seq2 are
    the putative and raw AA sequences as aligned by the NW aligner. If the status
    value is 0, that is essentially a flag that the sequences could not be neartly
    matched and require manual attention.
    :param aa_rw:
    :param nu_rw:
    :return:
    '''
    seq_status=0
    seq1=''
    seq2=''
    aa_putative = nuc2aa(nu_rw)
    if aa_rw==aa_putative:
        seq_status=1
        seq1 = aa_putative
        seq2 = aa_rw
    else:
        res = run_nw(aa_putative, aa_rw)
        if int(res[5])<=5:
            seq_status=2
            seq1=res[7]
            seq2=res[8]
        else:
            # res1 = nwo.get_pairwise_alignment(nuc2aa(nu_rw[1:]),aa_rw, Fmat, Dmat)
            res1 = run_nw(nuc2aa(nu_rw[1:]),aa_rw)
            if int(res1[5])<=5:
                seq_status=3
                seq1=res1[7]
                seq2=res1[8]
                return (seq1, seq2, seq_status)
            # res2 = nwo.get_pairwise_alignment(nuc2aa(nu_rw[2:]),aa_rw,Fmat, Dmat)
            res2 = run_nw(nuc2aa(nu_rw[2:]),aa_rw)
            if int(res2[5]) <=5:
                seq_status=4
                seq1=res2[7]
                seq2=res2[8]
    return (seq1, seq2, seq_status)

def get_seq_status_and_alignment_scorebased(aa_rw, nu_rw):
    '''
    Just check all three alignments and take the one with the best score.
    :param aa_rw:
    :param nu_rw:
    :return:
    '''
    seq_score = None
    seq1 = ''
    seq2 = ''
    aa_putative = nuc2aa(nu_rw)
    if abs(len(aa_rw)-len(aa_putative))/max(len(aa_putative),len(aa_rw))>0.5:
        return
    if aa_rw == aa_putative:
        seq_score = 99999
        seq1 = aa_putative
        seq2 = aa_rw
        fs = 0  # frameshift
        return (seq1, seq2, seq_score, fs, -1, -1, -1, -1)
    else:
        # return (res[7], res[8], res[6], int(res[9]))
        res = run_nw_frametest_pyext(nu_rw, aa_rw)
        # (seq1, seq2, score, frameshift, match_ct, gap_ct, diff_ct, aln_length)
        return (res[5], res[6], res[4], res[0], res[1], res[2], res[3], len(res[5]))


#****************************************
# END UTILITY FUNCTIONS
#****************************************

# def read_final_seqs_for_codon_addback():
#     seq_fix_fp = os.path.join('/projects/tallis/nute/data/tippref_2018', 'dna_to_aa_matching',
#                               'final_seqs_for_codon_alignment_repair.txt')
#     seqfile = open(seq_fix_fp,'r')
#
#     temp_partial_fasta_loc='/projects/tallis/nute/data/tippref_2018/COGS_aligned_nuc/partial'
#     headers=seqfile.readline().strip().split('\t')
#     headers_small = headers[0:5]
#     data={}
#     for ln in seqfile:
#         flds=ln.strip().split('\t')[0:5]
#         if flds[0] not in data:
#             data[flds[0]] = [tuple(flds),]
#         else:
#             data[flds[0]].append(tuple(flds))
#     seqfile.close()
#     nuc_aligned = ''
#     partial_fasta = open(os.path.join(temp_partial_fasta_loc, 'additional_aligned_nuc_seqs.fasta'), 'w')
#     for cog in list(data.keys()):
#         aa_raw, nu_raw, aa_aln, nm_old2new, nm_new2old = ts.get_raw_alnd_nuke_aa_fastadicts(cog)
#         ks = list(aa_aln.keys())
#         ks_old = list(map(lambda x: nm_new2old[ks[x]], range(len(ks))))
#
#         nuc_tgt_alignment_len=nucleo_alnmt_lens[cog]
#         for rec in data[cog]:
#             myind=int(rec[1])
#             offset=int(rec[4])-1
#             aarw=aa_raw[ks_old[myind]]
#             nurw=nu_raw[ks_old[myind]]
#             aa_putative=nuc2aa(nurw[offset:])
#             aa_aligned=aa_aln[ks[myind]]
#             res = run_nw(aa_putative, aarw)
#             s1=res[7]
#             s2=res[8]
#             nuc_aligned = get_aligned_nuc_sequence(aarw, nurw, aa_aligned, myind, s1,s2, nuc_tgt_alignment_len)
#             partial_fasta.write('>%s\n' % ks[myind])
#             partial_fasta.write(nuc_aligned)
#             partial_fasta.write('\n')
#     partial_fasta.close()
#     # print('>%s' % ks[myind])
#     # print(nuc_aligned[:30] + '(total length: %s)' % len(nuc_aligned))


def make_aligned_nuc_fasta(cog, show_progress=False):
    '''
    Runs the script for a single COG to get the main fasta codon-subbed file.
    :param cog:
    :return:
    '''
    aa_raw, nu_raw, aa_aln, nm_old2new, nm_new2old = get_raw_alnd_nuke_aa_fastadicts(cog)
    ks = list(aa_aln.keys())
    ks_old = list(map(lambda x: nm_new2old[ks[x]], range(len(ks))))

    # output file:
    fasta_fpath = os.path.join(tippref_fold, 'COGS_aligned_nuc', 'aln_%s.fna' % cog)
    fasta_dict = dict.fromkeys(ks,'')
    # fasta_f = open(fasta_fpath, 'w')

    # keeping track of the sequence status for each sequence.
    # seq_status = np.zeros(len(ks), dtype=np.uint8)
    # seq_status = dict.fromkeys(ks,None)
    seq_status = []
    aa_alignment_len = len(aa_aln[ks[0]])
    nuc_tgt_alignment_len = 3 * aa_alignment_len + 3  # includes room for stop codon
    # cog_num_taxa = len(ks)
    st_tic = datetime.datetime.now()
    tic = st_tic
    ln_ct=0
    rogue_seqs = []
    seqs_iter = map(lambda i: (aa_raw[ks_old[i]], nu_raw[ks_old[i]], aa_aln[ks[i]], i), range(len(ks)))

    for sqs in seqs_iter:
        # first get status to make sure it's ok...
        # sqs = (aa_raw[ks_old[i]], nu_raw[ks_old[i]], aa_aln[ks[i]], i)
        # seq_stat = get_seq_status_and_alignment(sqs[0], sqs[1])
        ln=len(sqs[1])
        if ln<10000:
            seq_stat = get_seq_status_and_alignment_scorebased(sqs[0], sqs[1])
        else:
            seq_stat = None

        if seq_stat is None:
            rogue_seqs.append(get_putative_to_rawaa_details(sqs[0], sqs[1], ks[sqs[3]], cog, ks_old[sqs[3]]))
            del fasta_dict[ks[sqs[3]]]
            continue
        seq_status.append(tuple([ks[sqs[3]],]) + tuple(seq_stat[2:]))  # record keeping vector


        # then get sequence, write it, and note the time if needed.
        nuc_aligned = get_aligned_nuc_sequence(*sqs, seq_stat[0], seq_stat[1], nuc_tgt_alignment_len, seq_stat[3])
        fasta_dict[ks[sqs[3]]]=nuc_aligned

        ln_ct += 1
        if show_progress and ln_ct % 500 == 0:
            toc=datetime.datetime.now()-tic
            toc_st = datetime.datetime.now()-st_tic
            tic=datetime.datetime.now()
            print('# seqs: %d\tLast 500: %0.5f\tTotal: %0.5f' % (ln_ct, td2flt(toc),td2flt(toc_st)))

    tot = datetime.datetime.now()-st_tic
    print("writing cog %s to a file, took %s seconds" % (cog, tot))
    write_to_fasta(fasta_fpath, fasta_dict)
    # fasta_f.close()
    write_list_to_file(list(map(lambda x: '%s\t%s\t%s\t%s\t%s\t%s\t%s' % x, seq_status)),
                       os.path.join(tippref_fold, 'COGS_aligned_nuc', 'seq_status_%s.txt' % cog))

    if len(rogue_seqs)>0:
        write_list_to_file(rogue_seqs, os.path.join(tippref_fold, 'COGS_aligned_nuc', 'rogue_seqs_%s.txt' % cog))

def get_putative_to_rawaa_details(aa, nu, seqname, cog, old_seqname):
    res1 = run_nw(nuc2aa(nu), aa)
    res2 = run_nw(nuc2aa(nu[1:]), aa)
    res3 = run_nw(nuc2aa(nu[2:]), aa)
    return ','.join(['%s',]*24) % ((seqname, old_seqname,cog) + tuple(res1[0:7]) + tuple(res2[0:7]) + tuple(res3[0:7]))


def testing_open_seq_status(cog):
    f=open(os.path.join(tippref_fold, 'COGS_aligned_nuc', 'seq_status_%s.txt' % cog),'r')
    ss = {}
    for ln in f:
        a=ln.strip().split('\t')
        ss[a[0]]=tuple(map(int,a[1].split(',')))
    f.close()


    return ss


if __name__=='__main__':
    # make_aligned_nuc_fasta('COG0085')
    # make_aligned_nuc_fasta(cogs[0])
    p = multiprocessing.Pool(20)
    p.map(make_aligned_nuc_fasta, cogs)