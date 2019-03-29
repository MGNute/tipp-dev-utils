import os, sys
from phylogeny_utilities.utilities import *
from collections import Counter
import numpy as np
import scipy.stats
#from score import *

# fa=read_from_fasta('pasta.fasta')
# owners=[k for k in fa.keys() if fa[k] == dupe_seqs[0]]


taxrank_colsdict = {
    'superkingdom': 2,
    'phylum': 7,
    'class': 11,
    'order': 16,
    'family': 20,
    'genus': 26,
    'species': 32
}

taxrank=[2, 7, 11, 16, 20, 26, 32]

cols=['root', 'root_', 'superkingdom', 'superkingdom_', 'superkingdom__', 'superkingdom___', 'superkingdom____',
      'phylum', 'phylum_', 'phylum__', 'subphylum', 'class', 'class_', 'class__', 'class___', 'subclass', 'order',
      'order_', 'order__', 'suborder', 'family', 'family_', 'family__', 'family___', 'subfamily', 'tribe', 'genus',
      'genus_', 'subgenus', 'species_group', 'species_subgroup', 'species_subgroup_', 'species', 'species_',
      'species__', 'species___', 'subspecies', 'subspecies_', 'subspecies__', 'subspecies___']

rank_cols=dict(zip(cols,range(1,len(cols)+1)))

def get_taxonomy_table(folder=None):
    # tax_f = open('/projects/tallis/nute/data/ncbi-2017/taxonomy_folder/taxonomy.table', 'r')
    # tax_f = open('/projects/tallis/nute/data/ncbi-2017/taxonomy_folder/taxonomy.table.with_GT', 'r')
    if folder is None:
        tax_f = open('taxonomy.table', 'r')
    else:
        tax_f = open(os.path.join(folder,'taxonomy.table'), 'r')
    tax_f.readline()
    taxa_dict = {}
    for ln in tax_f:
        a = ln.strip()[1:-1].split('\",\"')
        taxa_dict[int(a[0])] = a[1:]
    tax_f.close()
    return taxa_dict

def fill_parent_vector(mytaxid, tax_dict=None, specmap=None):
    if tax_dict is None:
        tax_dict = get_taxonomy_table()
    taxid = mytaxid
    new_taxid = None
    new_lev = None
    parent_vector = [-1, ] * 40
    lev = rank_cols[tax_dict[taxid][1]]
    parent_vector[lev - 1] = taxid
    while new_taxid != 1:
        new_taxid = int(tax_dict[taxid][0])             #also sharp
        new_lev = rank_cols[tax_dict[new_taxid][1]]   # the level, i.e. a number
        # fill in all the levels until new_lev
        assert new_lev < lev, 'somehow the rank went down'
        for i in range(lev - new_lev):
            parent_vector[lev - i - 1] = taxid
            parent_vector[new_lev - 1] = new_taxid
        taxid = new_taxid
        lev = new_lev
    return parent_vector

def get_taxonomic_diversity_for_sequence(sq, fa, line_lead = '', tax_tab=None, specmap=None):
    '''
    Finds the distribution and calculates summary statistics of the set of taxa in a fasta
    that shares a sequence.
    :param sq: (string). a seqeunce whose detailed lineage are looking for.
    :param fa: a dictionary object representing a FASTA file.
    :return: Returns the percentage of the sequences taht belong to the most dominant taxon id,
    for each of the 7 major taxonomic levels. Also returns some basic metadata.
    '''

    # find the taxa that share this sequence:
    owners = [k for k in fa.keys() if fa[k]==sq]
    taxa_np = np.ones((len(owners),40),dtype=np.int32)*-1
    num_owners = len(owners)

    # debugging
    # print(str(num_owners))
    # print("length of owners vector: %s" % len(owners))
    # print(str(owners))

    for i in range(len(owners)):
        taxa_np[i,:] = fill_parent_vector(int(specmap[owners[i]]), tax_tab,specmap)

        print(line_lead + '...getting taxon info for genome %d of %d' % (i, num_owners),end='\r',flush=True)
        # if verbose and i > cut:
        #     print('finished with %s owners out of %s' % (cut,len(owners)))
        #     cut += cut_gap
    print(line_lead + '...DONE...                                  ', end='\r', flush=True)

    ranks = ['superkingdom','phylum','class','order','family','genus','species']
    ranks_2char = {'superkingdom':'sk','phylum':'ph','class':'cl','order':'or','family':'fa','genus':'ge','species':'sp'}

    dominant_taxa = dict.fromkeys(ranks,None)
    for k in ranks:
        col=taxrank_colsdict[k]
        md = scipy.stats.mode(taxa_np[:,col])
        dominant_taxa[k] = (md[0][0], md[1][0], md[1][0] / num_owners)  # 3-tuple with (value, count, % of column)
        dominant_taxa[ranks_2char[k]+'_str'] = ranks_2char[k] + ':(%d, %d, %.2f)' % dominant_taxa[k]

    dominant_taxa['rep_seq_name'] = owners[0]
    dominant_taxa['num_owners'] = num_owners

    return dominant_taxa


def review_fasta_file_with_duplicates(fasta_path, debug = False, initial_only = False):
    '''
    Basically what we want to know is if a a big group of identical sequences are in here, do they belong to the
    same (basic) organism.
    :param fasta_path:
    :return:
    '''
    print('...opening file and gathering initial data...',end='\r',flush=True)
    fp_full = os.path.abspath(fasta_path)
    #debugging
    # print('\t\tfp_full: %s' % fp_full)

    fp_fold, fp_file = os.path.split(fp_full)
    # print('\t\tfp_fold: %s' % fp_fold)
    # print('\t\tfp_file: %s' % fp_file)

    assert ('species.mapping' in os.listdir(fp_fold)), "The location of the alignment must also have the species mapping"
    specmap = get_dict_from_file(os.path.join(fp_fold,'species.mapping'), delimiter=",")
    # print('specmap has length: %s' % len(specmap))


    tax_tab = get_taxonomy_table(fp_fold)

    fa=read_from_fasta(fasta_path)
    c=Counter(list(fa.values()))
    dupe_seqs = [k for k in c.keys() if c[k]>1]
    metadata= {
        'full_fasta_path': fp_full,
        'num_seqs':len(list(fa.keys())),
        'num_cols':len(list(fa.values())[0]),
        'num_unique_seqs':len(c),
        'num_dupe_seqs':len(dupe_seqs)
    }

    # Sort the dupe list so largest dupes are first:
    dupe_seqs.sort(reverse=True, key=lambda x: c[x])
    dupe_cts = list(map(lambda x: c[x], dupe_seqs))
    print('                                                ', end='\r', flush=True)

    cutoff = 50     # if a sequence doesn't show up at least this many times, we're not
                    # going to worry about whether multiple taxa show up.

    num_seqs_to_check = sum([1 for k in dupe_cts if k > cutoff])
    results = {}
    print('Fasta File: %s' % fp_full)
    print('    -# seqs:        %5d,    -# cols:        %s' % (metadata['num_seqs'], metadata['num_cols']))
    print('    -# unique seqs: %5d,    -# duplicates:  %s' % ( metadata['num_unique_seqs'],  metadata['num_dupe_seqs']))
    print('    -cutoff: %s copies        -# over cutoff: %s' % (cutoff, num_seqs_to_check))
    if initial_only is True:
        print('\t(option for initial information only was given, stopping here)')
        sys.exit(0)
    print('\n')
    print(' ###)  #copies:,   dominant taxon info (Tax ID, Freq (ct), % of Tot)...')
        #'%4d)  %5d cp,'


    for i in range(num_seqs_to_check):
        # print(' ###) #copies:,
        ln_lead = '%4d) %5d cp,' % (i,dupe_cts[i])
        dt = get_taxonomic_diversity_for_sequence(dupe_seqs[i], fa, ln_lead, tax_tab,  specmap)
        tax_str = ', '.join([dt['ph_str'], dt['cl_str'], dt['or_str'], dt['fa_str'], dt['ge_str'], dt['sp_str']])
        print(ln_lead + ' ' + tax_str)
        if debug and i >=1:
            break


if __name__=='__main__':
    fp = sys.argv[1]
    if '-i' in sys.argv or '--initial_only' in sys.argv:
        review_fasta_file_with_duplicates(fp, False, True)
        sys.exit(0)
    review_fasta_file_with_duplicates(fp,False)
