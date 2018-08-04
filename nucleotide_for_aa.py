from phylogeny_utilities.utilities import *
import subprocess as sp
import sys
from multiprocessing import Pool

codon_lookup = {'GCT':'A','GCC':'A','GCA':'A','GCG':'A','CGT':'R','CGC':'R','CGA':'R','CGG':'R',
                    'AGA':'R','AGG':'R','AAT':'N','AAC':'N','GAT':'D','GAC':'D','TGT':'C','TGC':'C',
                    'CAA':'Q','CAG':'Q','GAA':'E','GAG':'E','GGT':'G','GGC':'G','GGA':'G','GGG':'G',
                    'CAT':'H','CAC':'H','ATT':'I','ATC':'I','ATA':'I','TTA':'L','TTG':'L','CTT':'L',
                    'CTC':'L','CTA':'L','CTG':'L','AAA':'K','AAG':'K','ATG':'M','TTT':'F','TTC':'F',
                    'CCT':'P','CCC':'P','CCA':'P','CCG':'P','TCT':'S','TCC':'S','TCA':'S','TCG':'S',
                    'AGT':'S','AGC':'S','ACT':'T','ACC':'T','ACA':'T','ACG':'T','TGG':'W','TAT':'Y',
                    'TAC':'Y','GTT':'V','GTC':'V','GTA':'V','GTG':'V','TAA':'(stop)','TGA':'(stop)',
                    'TAG':'(stop)','TRA':'(stop)','TAR':'(stop)'}
reverse_lookup = {  'F': 'TT*',
                    'L': '*T*',
                    'I': 'AT*',
                    'M': 'ATG',
                    'V': 'GT*',
                    'S': 'TC*',
                    'P': 'CC*',
                    'T': 'AC*',
                    'A': 'GC*',
                    'Y': 'TA*',
                    'H': 'CA*',
                    'Q': 'CA*',
                    'N': 'AA*',
                    'K': 'AA*',
                    'D': 'GA*',
                    'E': 'GA*',
                    'C': 'TG*',
                    'W': 'TGG',
                    'R': '*G*',
                    'G': 'GG*' #Serine has been assumed here to be one of the TC* codongs but it asly has AG*
                }

# def string_aa_to_nukes(aa_aln, nuc_raw):
#     out_str = ''
#     len_aa = len(aa_aln.replace('-',''))
#     err_ct = 0
#     # if len(nuc_raw)/3==len_aa+1:
#     #     stop_missing = False
#     # else:
#     #     print 'lengths of strings do not add up: AA length = %s, nuke length = %s (%s)' % (len_aa, len(nuc_raw), len(nuc_raw)/3)
#     #     if len(nuc_raw)/3==len_aa and nuc_raw[0:3].upper() in ['ATG','GTG']:
#     #         print '\tassuming stop codon is missing, moving on...'
#     #         stop_missing = True
#     #     else:
#     #         print 'not plausibly a stop codon issue, terminating...'
#     #         return False, 0
#     position = 0
#     for i in range(len(aa_aln)):
#         if aa_aln[i]=='-':
#             out_str += '---'
#         else:
#             nuc_temp = nuc_raw[(position*3):(position*3+3)].upper()
#             if nuc_temp in codon_lookup.keys() and codon_lookup[nuc_temp]<>aa_aln[i].upper():
#                 err_ct += 1
#                 # print 'amino acid at position %s (%s) does not match codon %s' % (position, aa_aln[i],nuc_temp)
#                 # return False
#             out_str += nuc_temp
#             position += 1
#     if not stop_missing:
#         last_codon=nuc_raw[-3:].upper()
#         if codon_lookup[last_codon]<>'(stop)':
#             err_ct += 1
#             print 'final codon (%s) is not a stop codon (err ct %s)' % (last_codon, err_ct)
#             # return False
#         out_str += last_codon
#     else:
#         out_str += '---'
#     # if err_ct > 0:
#     #     print 'error count: %s' % err_ct
#     return out_str, err_ct
def string_aa_to_nukes(aa_aln, nuc_raw):
    out_str = ''
    pseudo_gene = get_pseudo_nuke_string_from_aa(aa_aln.replace('-',''))
    p = sp.Popen(['/projects/tallis/nute/code/misc_c_utilities/nw'], stdin=sp.PIPE, stdout=sp.PIPE, stderr=sp.PIPE)
    instr = pseudo_gene + '\n' + nuc_raw + '\n'
    ps_out = p.communicate(input=instr)
    ps_out_l = ps_out[0].strip().split('\t')

    pseudo_gene_aln = ps_out_l[0]
    nuke_aln = ps_out_l[1]
    # aa_position = 0
    pseudo_gene_position=0
    while (pseudo_gene_aln[pseudo_gene_position]=='-'):
        pseudo_gene_position += 1

    def get_pseudo_gene_next(curr_pos):
        ct = 0
        marker = curr_pos
        while (ct < 3):
            if pseudo_gene_aln[marker]!='-':
                ct += 1
            marker += 1
        return marker

    for i in range(len(aa_aln)):
        if aa_aln[i]=='-':
            out_str += '---'
        else:
            fr = pseudo_gene_position
            to = get_pseudo_gene_next(fr)
            out_str += nuke_aln[fr:to]
            pseudo_gene_position = to

    if len(nuke_aln)>pseudo_gene_position:
        out_str += nuke_aln[(pseudo_gene_position-1):]
    return out_str




def get_pseudo_nuke_string_from_aa(aminostring):
    out_str = ''
    for i in range(len(aminostring)):
        if aminostring[i] in reverse_lookup.keys():
            out_str += reverse_lookup[aminostring[i]]
        else:
            out_str += '***'
    return out_str

def make_nuke_alignment_from_aa(aminofile, nucfile, nuc_aligned_out):
    mypool=Pool(20)
    aa = read_from_fasta(aminofile)
    nuc_raw = read_from_fasta(nucfile)
    if len(aa.keys())<>len(nuc_raw.keys()):
        print 'Sequence counts do not match between files: %s in AA file, %s in nucleotide file...' % (len(aa.keys()),len(nuc_raw.keys()))

    mm_names =len(set(aa.keys()).symmetric_difference(set(nuc_raw.keys())))
    if mm_names<>0:
        # print 'Sequence names are not exactly the same between files: %s names in only one file...terminating'
        intersect_keys = list(set(aa.keys()).intersection(set(nuc_raw.keys())))
        print 'Sequence names are not exactly the same between files. Output will match the intersection of the two.'
        out_nuc = dict.fromkeys(intersect_keys)
    else:
        out_nuc = dict.fromkeys(aa.keys())
    seq_ct = len(out_nuc.keys())
    print 'total of %s sequences' % seq_ct
    ct = 0
    arglist = []
    for k in out_nuc.keys():
        # ec = 0
        arglist.append((aa[k],nuc_raw[k],k))
        # newstr = string_aa_to_nukes(aa[k],nuc_raw[k])
    print 'starting pool....'
    res = mypool.map(string_aa_to_nukes_intermediate,arglist)
    print '\t\t...done running pool'

    for i in range(len(arglist)):
        k=arglist[i][2]
        if res[i] == False:
            print 'error occured on key %s...terminating' % k
            return False
        out_nuc[k] = res[i]
        # ct +=1
        # if ct % 10000 == 0:
        #     print '%s done' % ct

    maxlen = max(map(len,out_nuc.values()))
    print "max length: %s" % maxlen
    for i in out_nuc.keys():
        leni = len(out_nuc[i])
        out_nuc[i] += '-'*(maxlen-leni)
    print "min length: %s" % min(map(len,out_nuc.values()))
    # assert min(map(len,out_nuc.values()))==max(map(len,out_nuc.values())), 'final fasta does not all have same values'
    write_to_fasta(nuc_aligned_out, out_nuc)

def string_aa_to_nukes_intermediate(args):
    return string_aa_to_nukes(args[0],args[1])

def trim_fastas(fna, faa):
    cnu = read_from_fasta(fna)
    caa = read_from_fasta(faa)
    removal = []
    ct = 0
    for i in cnu.keys():
        if len(cnu[i])>1000000:
            removal.append(i)
            ct += 1
    if ct>0:
        print 'removing %s sequences: %s' % (ct,removal)
        ks = cnu.keys()
        for j in removal:
            ks.pop(ks.index(j))
        write_to_fasta(fna, cnu, ks)
        write_to_fasta(faa, caa, ks)

def split_bimodal_aligments():
    splits={'COG0012':1120,
            'COG0016':1020,
            'COG0018':1710,
            'COG0048':None,
            'COG0049':None,
            'COG0052':800,
            'COG0080':460,
            'COG0081':None,
            'COG0085':3900,
            'COG0087':600,
            'COG0088':630,
            'COG0090':800,
            'COG0091':500,
            'COG0092':775,
            'COG0093':380,
            'COG0094':550,
            'COG0096':None,
            'COG0097':None,
            'COG0098':550,
            'COG0099':400,
            'COG0100':400,
            'COG0102':None,
            'COG0103':430,
            'COG0124':1400,
            'COG0172':1325,
            'COG0184':300,
            'COG0185':300,
            'COG0186':300,
            'COG0197':430,
            'COG0200':400,
            'COG0201':1400,
            'COG0202':850,
            'COG0215':1500,
            'COG0256':400,
            'COG0495':2750,
            'COG0522':610,
            'COG0525':None,
            'COG0533':1500,
            'COG0541':1500,
            'COG0552':1375}
    for k in splits.keys():
        if splits[k] is None:
            continue
        sp = int(float(splits[k])/3.)
        fa = read_from_fasta(k + '.faa')
        fa_sm = {}
        fa_lg = {}
        for seq in fa.keys():
            if len(fa[seq])>sp:
                fa_lg[seq] = fa[seq]
            else:
                fa_sm[seq] = fa[seq]
        write_to_fasta('bimodal_split/' + k + '_large.faa',fa_lg)
        write_to_fasta('bimodal_split/' + k + '_small.faa', fa_sm)
        del fa
        del fa_sm
        del fa_lg

if __name__=='__main__':
    if len(sys.argv)<4:
        aafi = '-h'
        nucfi = ''
        outfi = ''
    else:
        aafi = sys.argv[1]
        nucfi = sys.argv[2]
        outfi = sys.argv[3]

    if '-h' in [aafi, nucfi, outfi]:
        print '''
        Script to convert an aligned amino acid fasta file to an aligned nucleotide file on a codon-by-codon basis
        given corresponding raw files. Does NOT check that codons and amino acids match perfectly.

        usage: python nucleotide_for_aa.py <amino_acid_fasta> <nucleotide_fasta> <output_fasta>'''
        sys.exit(0)
    else:
        make_nuke_alignment_from_aa(aafi,nucfi,outfi)
        # print '%s\n%s\n%s' % (aafi, nucfi, outfi)