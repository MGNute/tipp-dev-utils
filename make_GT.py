'''
Makes the ground truth taxonomy assignments for every sequence in the TIPP references. This
is mostly for scoring and evaluating performance.
'''
from phylogeny_utilities.utilities import *
import os
cogs=[i for i in os.listdir('.') if i[0:3]=='COG']
ranks = ['tax_name','superkingdom','phylum','class','order','family','genus','species']
lnstr = '\t'.join(['%s',] * 10) + '\n'

cogs = ['COG0012','COG0016','COG0018','COG0048','COG0049','COG0052',
            'COG0080','COG0081','COG0085','COG0087','COG0088','COG0090',
            'COG0091','COG0092','COG0093','COG0094','COG0096','COG0097',
            'COG0098','COG0099','COG0100','COG0102','COG0103','COG0124',
            'COG0172','COG0184','COG0185','COG0186','COG0197','COG0200',
            'COG0201','COG0202','COG0215','COG0256','COG0495','COG0522',
            'COG0525','COG0533','COG0541','COG0552']

work = '/projects/tallis/nute/work/metagenomics/tippstar/making_ground_truth_seq_locations'
ncbi = '/projects/tallis/nute/data/ncbi-2017'
gbff_loc = '/projects/tallis/nute/work/metagenomics/tippstar/2017'

def get_files_to_accn_mapping():
    '''

    :return:
    '''
    src = '/projects/tallis/nute/work/metagenomics/tippstar/making_ground_truth_seq_locations/assemblies_in_both_data_sources.txt'
    sfi = open(src,'r')
    sfi.readline()
    accn_to_files = {}
    for ln in sfi:
        if len(ln.strip())>0:
            a=ln.strip().split('\t')
            accn_to_files[a[0]] = {'assembly': a[1], 'file': a[2], 'id': int(a[3])}

    sfi.close()
    return accn_to_files




def add_cog_to_groundtruth_file(cog,gtfile):
    '''
    This is for making the mapping between COG sequences and ground truth NCBI taxonomy
    IDS
    :param cog:
    :param gtfile:
    :return:
    '''

    smap = get_list_from_file(os.path.join(cog,'species.mapping'))
    smap.pop(0)
    smap = map(lambda x: x.split(','),smap)
    smap = dict(map(lambda x: (x[0] + '_' + cog,int(x[1])),smap))

    atf = open(os.path.join(cog,'all_taxon.taxonomy'),'r')
    headers=atf.readline()
    headers=headers.strip()[1:-1].split('\",\"')
    a=headers.pop(0)
    inds = map(lambda x: headers.index(x),ranks)

    atlkp = {}
    for ln in atf:
        tax=ln.strip()[1:-1].split('\",\"')
        atlkp[int(tax[0])]=tuple(map(lambda x: tax[1:][x],inds))

    atf.close()
    ct=0
    for sq in smap.keys():
        v = smap[sq]
        ct+=1
        gtfile.write(lnstr % ((sq,v)+atlkp[v]))
    print ('wrote %s lines for cog %s' % (ct,cog))


