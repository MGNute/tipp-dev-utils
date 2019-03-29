'''
Miscellaneous functions created when building the 2019 reference packages.
'''

__author__ = 'Michael'
import os, sys
import numpy as np
import dendropy
from phylogeny_utilities.utilities import *

#
# Some relevant locations & variables
#
tippref_fold='/projects/tallis/nute/data/tippref_2018'
outfold=os.path.join(tippref_fold,'COGS_sequence_metadata')

# lookup files (see notes later):
file_ids = os.path.join(tippref_fold, 'ncbi_refseq_intermediate_files', 'file-list-lookup.txt')
seq_nm_pos=os.path.join(tippref_fold, 'ncbi_refseq_intermediate_files', 'ncbi_all_sequences_name_pos.txt')

# this one has a set of name-maps for each COG, mapping from the old name ('<file-number>_<ncbi-protein-id>')
#   to the new name 's##_COGXXXX' where the ## is just a numeric identifier unique to the COG.
cogs_name_map_fold=os.path.join(tippref_fold,'COGS_name_map')
cogs = ['COG0012','COG0016','COG0018','COG0048','COG0049','COG0052',
        'COG0080','COG0081','COG0087','COG0088','COG0090','COG0085',
        'COG0091','COG0092','COG0093','COG0094','COG0096','COG0097',
        'COG0098','COG0099','COG0100','COG0102','COG0103','COG0124',
        'COG0172','COG0184','COG0185','COG0186','COG0197','COG0200',
        'COG0201','COG0202','COG0215','COG0256','COG0495','COG0522',
        'COG0525','COG0533','COG0541','COG0552']

# These three will be various lookups imported from files
#   whose values get assigned in the function that follows.
fll = None
fll_rev = None
file_to_assembly_accn = None
old_seq_to_segID = {}

fa_path='pasta.fasta'
tr_path='pasta.tree'
tr_out_path='pasta.tree.dedup'
dedup_data_out_path='dedup_name_map.json'
import os, json

def make_dedup_tree_and_mappingJSON(fa_path, tr_path, tr_out_path, dedup_data_out_path):
    '''

    :param fa_path:
    :param tr_path:
    :param tr_out_path:
    :param dedup_data_out_path:
    :return:
    '''
    fa=read_from_fasta(fa_path)
    tr=dendropy.Tree.get(path=tr_path, schema='newick', preserve_underscores=True)
    fa_n2c, fa_eq = get_fasta_duplicate_datastruct(fa)
    tax_labs = list(map(lambda x: fa_eq[x]['members'][0], fa_eq.keys()))
    tax_lab_to_class = dict(map(lambda x: (x, fa_n2c[x]),tax_labs))
    tr_dedup = tr.extract_tree_with_taxa_labels(tax_labs)
    tr_dedup.write(path=tr_out_path, schema='newick')
    dd_f = open(dedup_data_out_path, 'w')
    json.dump({
            'fasta_names_to_equiv_class':fa_n2c,
            'fasta_equivalence_class_definitions':fa_eq,
            'deduped_name_to_equivalence_class': tax_lab_to_class
            }, dd_f)
    dd_f.close()






def make_lookups():
    print("\t...making file id lookups")
    # it was necessary for brevity to assign a numeric identifier to each of the refseq files that
    #   were used to create the reference. The mapping is contained in the file below. For purposes
    #   later, we will import this both as a forward and a reverse lookup:

    fll = get_dict_from_file(file_ids)
    fll_rev = dict(map(lambda x: (x[1],x[0]), fll.items()))

    # it will also be convenient to have a mapping from the file number to the NCBI assembly
    #   accession ID:
    file_to_assembly_accn = dict(map(lambda x: (x,fll[x][:15]),fll.keys()))

    # This one is a gigantic lookup. It maps '<file-number>_<ncbi-protein-id>' to the NCBI segment
    #   from which it came. Prints a status update every million lines.
    print("\t...making old_seq_to_SegID lookup")
    seq_nm_pos_f = open(seq_nm_pos,'r')
    ct = 0
    for ln in seq_nm_pos_f:
        a=ln.strip().split('\t')
        old_seq_to_segID[fll_rev[a[0]] + '_' + a[1]] = a[2]
        ct +=1
        if ct % 1000000 == 0:
            print("%s M done     " % int(ct/1000000), end='\r')

    print("")
    print("Done")
    seq_nm_pos_f.close()

'''
Since there were so many sequence duplicates in the COGS, I was curious to know how much
    of the cross-validation holdout sets would be potentially still in the training data due
    to sequence redundancy. The two functions below were used to do some of that analysis.
    
NOTE: A couple other useful functions for this are contained in the file taxonomic_diversity.py
'''
# first get some lookups
test_ass_accns = None
f_assemb_dupes = None
def checking_dupe_effects_make_lookups():
    '''

    :return:
    '''
    # make DB of dupe counts
    assembs='/projects/tallis/nute/work/metagenomics/tipp2new/assembly_crossval_sets.txt'
    test_assembs=get_dict_from_file(assembs)
    global test_ass_accns
    global f_assemb_dupes
    test_ass_accns=list(test_assembs.keys())

    pa_assemb_dupes = os.path.join(outfold,'all_assembly_dupe_cts.txt')
    f_assemb_dupes = open(pa_assemb_dupes,'w')
    lct = f_assemb_dupes.write('assembly\tcog\tseq\teq_class\tcopynum\tsegID\n')

def checking_dupe_effects_make_cog_metadata(mycog):
    '''
    This function created a metadata file that I still have, although I'm not sure why
    it was useful at the time.
    :param mycog:
    :return:
    '''
    global test_ass_accns
    global f_assemb_dupes
    out_fn = os.path.join(outfold, mycog + '.txt')
    print("creating sequence metadata file at: %s" % out_fn)
    #
    print("\t...making name map lookup")
    nm_path = os.path.join(cogs_name_map_fold, mycog + '.txt')
    name_map = get_dict_from_file(nm_path, keysfirst=False)
    #
    fasta_dict = read_from_fasta(os.path.join(tippref_fold, 'refpkg', mycog + '.refpkg', 'pasta.fasta'))
    dupes_seq_nm_to_class, dupes_eq_classes = get_fasta_duplicate_datastruct(fasta_dict)
    assemb_to_seqs = {}
    #
    f_out = open(out_fn, 'w')
    lct = f_out.write('seq\tassemblyAccn\tSegmentID\tEquivClass\n')
    #
    for k in fasta_dict.keys():
        file_num = name_map[k].split('_')[0]
        assemb = file_to_assembly_accn[file_num]
        seg = old_seq_to_segID[name_map[k]]
        eq = dupes_seq_nm_to_class[k]
        lct = f_out.write('%s\t%s\t%s\t%s\n' % (k, assemb, seg, eq))
        if assemb in assemb_to_seqs:
            assemb_to_seqs[assemb].append(k)
        else:
            assemb_to_seqs[assemb] = [k, ]
    #
    print("\nDone")
    f_out.close()
    not_in_list = []

    for aa in test_ass_accns:
        if aa not in assemb_to_seqs:
            # print("assembly %s is not included in cog %s" % (aa, mycog))
            not_in_list.append(aa)
            continue
        ns = len(assemb_to_seqs[aa])
        for i in range(ns):
            s = assemb_to_seqs[aa][i]
            ec = dupes_seq_nm_to_class[s]
            cn = dupes_eq_classes[ec]['copynum']
            sg = old_seq_to_segID[name_map[s]]
            lct = f_assemb_dupes.write('%s\t%s\t%s\t%s\t%s\t%s\n' % (aa, mycog, s, ec, cn, sg))
            # dupes_eq_classes[dupes_seq_nm_to_class[assemb_to_seqs[a][0]]]['copynum']

    print("%s test accessions not in cog %s:" % (len(not_in_list), mycog))
    print(not_in_list)
    for aa in not_in_list:
        lct = f_assemb_dupes.write('%s\t%s\t%s\t%s\t%s\t%s\n' % (aa, mycog, None, None, None, None))
    print("\n********\n")

'''
All of the code below was created to generate the holdout datasets for the 2019 references for testing.
'''
# some relevant locations:
tipp2_fold = '/projects/tallis/nute/work/metagenomics/tipp2new'
full_refpkg ='/projects/tallis/nute/data/tipp-2018'
holdout = os.path.join(tippref_fold,'holdout')

def CV_get_seqs_to_exclude():
    '''
    For each holdout set and each COG, find the sequences that must be excluded from the CV training
        set.

    Creates a dictionary of dictionaries of lists of sequences:
    { 0: {  'COG0012': [s1, s2, ...],
            'COG0016': [....], ...
            },
      1: {  'COG0012': [s1, s2, ...],
            'COG0016': [....], ...
            },
      2: ...,
      3: ...
    }
    :return:
    '''
    excl= {
        0: dict.fromkeys(cogs, None),
        1: dict.fromkeys(cogs, None),
        2: dict.fromkeys(cogs, None),
        3: dict.fromkeys(cogs, None),
    }
    for i in range(4):
        for j in cogs:
            excl[i][j] = []
    list_file = os.path.join(tipp2_fold, 'ref_crossval_seqs_to_remove.txt')
    list_file_f = open(list_file, 'r')
    for ln in list_file_f:
        a=ln.strip().split('\t')
        excl[int(a[0])][a[1]].append(a[2])
    return excl

def CV_copy_ref_for_holdout(cog, fold, excl_d=None):
    '''
    NOTE: THIS FUNCTION IS INCOMPLETE!!
    Makes a copy of the reference package that excludes the holdout data
    :param cog:
    :param fold:
    :param excl_d: dictionary created in the function above.
    :return:
    '''
    if excl_d is None:
        excl_d = get_seqs_to_exclude()
    fold=str(fold)
    orig=os.path.join(full_refpkg, cog + '.refpkg')
    newf=os.path.join(holdout,'f'+fold, cog + '.refpkg')

    #first copy the easy stuff, which requires no changing:
    easy_files = ['all_taxon.taxonomy', 'pasta.hmm', 'pasta.taxonomy.RAxML_info','pasta.tree.RAxML_info']
    for f in easy_files:
        od = os.path.join(orig,f)
        nw = os.path.join(newf,f)
        os.system('cp %s %s' % (od, nw))

    # TODO: Add code to make the pasta copy and the tree copies:

def CV_make_pasta_subalignment(cog,part, excl_d=None):
    '''
    Makes the pasta.fasta files for the CV references.
    :param cog:
    :param part:
    :param excl_d:
    :return:
    '''
    if excl_d is None:
        excl_d = get_seqs_to_exclude()
    part=str(part)
    seqs2excl = excl_d[part][cog]
    origf=os.path.join(full_refpkg, cog + '.refpkg')
    newf=os.path.join(holdout,'f'+fold, cog + '.refpkg')

    old_fasta = read_from_fasta(os.path.join(origf,'pasta.fasta'))
    old_seq2eq, old_eq = get_fasta_duplicate_datastruct(old_fasta)

    # (1) - write out the new pasta.fasta
    new_fasta = dict(filter(lambda x: False if x[0] in seqs2excl else True, old_fasta.items()))
    print ('Extracting PASTA subalignment for %s, part %s' % (cog,part))
    new_fasta = remove_all_blank_columns_utils(new_fasta)
    write_to_fasta(os.path.join(newf, 'pasta.fasta'), new_fasta)

    # (2) - update the duplicate data for the reduced structure
    fasta_reduced=read_from_fasta('pasta.fasta.reduced.fasta')
    new_eq = old_eq.copy()
    for s in seqs2excl:
        new_eq[old_seq2eq[s]]['members'].remove(s)
        new_eq[old_seq2eq[s]]['copynum']-=1

def test_codon_aa_matching(cog):
    '''
    This function is unfortunately where the wheels came off in february 2019. Some strange
    behavior prompted me to double check that the translation from the AA alignment to the
    nucleotide alignment was done correction, and sure enough it was not.
    :param cog:
    :return:
    '''
    faa = read_from_fasta(tippref_fold+'/COGS_aligned_aa/aln_%s.faa' % cog)
    fan = read_from_fasta(tippref_fold+'/COGS_aligned_aa/aln_%s.fna' % cog)
    aatax, aafanp = fasta_dict_to_nparray(faa)
    tax, fanp = fasta_dict_to_nparray(fan)
    print('shapes.. aa: (%s,%s)\t' % aafanp.shape,end='')
    print('shapes.. nu: (%s,%s)\t' % fanp.shape)
    return aatax, aafanp, tax, fanp

#
#   for command line execution
#
def main():
    for cg in cogs:
        make_cog_metadata(cg)

if __name__=='__main__':
    main()