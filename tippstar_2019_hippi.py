'''
Functions written in the process of evaluating the hippi experiments. The three main jobs that these functions
have to do is 1) read and interpret the ground truth data from Nidhi, 2) read the nhmmer outputs generated when
HIPPI ran and summarize and interpret those, 3) same as (2), but for BLAST.

The last job is to marry all those together
'''


#
# 4/1/19
#
import os, sys, math
from phylogeny_utilities.utilities import *
from .hippi_nhmmer_dispatcher import *

# Folders & Locations:
tipp2='/projects/tallis/nute/work/metagenomics/tipp2new'    # /projects/tallis/nute/work/metagenomics/tipp2new
gt_fold=os.path.join(tipp2,'ground_truth')                  # /projects/tallis/nute/work/metagenomics/tipp2new/ground_truth
gt_subfolds=['gt_454','gt_ill100','gt_ill250','gt_pacbio']
blast_results_fold_fulldata = os.path.join(tipp2, 'blast_results')   # /projects/tallis/nute/work/metagenomics/tipp2new/blast_results
blast_results_fold_ho = os.path.join(tipp2, 'blast_results_ho')   # /projects/tallis/nute/work/metagenomics/tipp2new/blast_results_ho
# blast_results_fold = blast_results_fold_ho

f2a=os.path.join(tipp2,'file_to_assembly_map.txt')          # /projects/tallis/nute/work/metagenomics/tipp2new/file_to_assembly_map.txt
hippi_results_fold_fulldata = os.path.join(tipp2, 'hippi_results', 'all_results')    # /projects/tallis/nute/work/metagenomics/tipp2new/hippi_results/all_results
hippi_results_fold_ho = os.path.join(tipp2, 'hippi_results_ho', 'all_results')    # /projects/tallis/nute/work/metagenomics/tipp2new/hippi_results_ho/all_results
# hippi_results_fold = hippi_results_fold_ho

'''
A "qs_id" as it is used here refers to the ID of a set of query sequences. Specifically in 
this usage it is something like "s454_f0_seg1", which would refer to the query sequences from
the 454 simulation, fold 0, segment 1 (each fold is broken into segments for size reasons).
'''
qs_ids = [i[:-4] for i in os.listdir(hippi_results_fold_fulldata) if i[0]=='s']
qs_ids_454 = [i for i in qs_ids if i[:4]=='s454']
qs_ids_ill100 = [i for i in qs_ids if i[:7]=='sill100']
qs_ids_ill250 = [i for i in qs_ids if i[:7]=='sill250']
qs_ids_pacbio = [i for i in qs_ids if i[:7]=='spacbio']

def get_gt_from_nidhi_results(fipath):
    '''
    For a single assembly and modality, it pulls a short list of simulated reads that match
    the marker genes of interest.

    Here, fipath is the full path to one of the 'assignment' files from Nidhi, where each line
    represents a single read. The function scans through and if the assignment is one of the
    COGS of interest, it adds that read and its assignment to the dictionary.
    :param fipath:
    :return: dictionary object as:

        { '<readID>' : ( UniprotOrFetchMG, cogID),
            ...
        }

    '''
    f=open(fipath, 'r')
    fo, fi = os.path.split(fipath)
    if 'pacbio' in fi:
        read_ind = 4
        cog_ind = 14
    else:
        read_ind = 3
        cog_ind = 12

    res={}
    for ln in f:
        a=ln.strip().split('\t')
        if a[cog_ind] in ['COG0088','COG0090','COG0094']:
            res[a[read_ind]] = (a[cog_ind - 1],a[cog_ind])
    f.close()
    return res


def get_full_gt_lookup_from_subfold(subfold):
    '''

    :param subfold: should be one of the strings in the set "gt_subfolds" given earlier. One option
    for each of the 4 sequencing platforms.
    :return all_res: Dictionary keyed by '<accn>_<readID>'. Values are also tuples that come from the
    values in the dictionary output by the function above. So in total, we have:

        {   '<accn-1>_<read-1>' : ( UniprotOrFetchMG, cogID, <accn-1>, <read-1> ),
            '<accn-1>_<read-2>' : ( UniprotOrFetchMG, cogID, <accn-1>, <read-2> ),
            ...
            '<accn-1>_<read-{k_1}>' : ( UniprotOrFetchMG, cogID, <accn-1>, <read-{k_1}> ),
            '<accn-2>_<read-1>' : ( UniprotOrFetchMG, cogID, <accn-2>, <read-1> ),
            '<accn-1>_<read-2>' : ( UniprotOrFetchMG, cogID, <accn-2>, <read-2> ),
            ...
            '<accn-2>_<read-{k_2}>' : ( UniprotOrFetchMG, cogID, <accn-2>, <read-{k_2}> ),
            ....
        }

    OLD STRUCTURE OF THE OUTPUT (changed 5/3):
        { assembly_accession_1 : {  read_1 : ('both', <cog_a>) ,
                                    read_2 : ('both', <cog_b>) ,
                                    ... } ,
          assembly_accession_2 : {  read_1 : ('both', <cog_a>) ,
                                    read_2 : ('both', <cog_b>) ,
                                    ... } ,
          ...
        }
    '''
    fold=os.path.join(gt_fold, subfold)
    sequencer = subfold.split('_')[1]   # one of: '454', 'ill100', 'ill250', 'pacbio'
    delstr = '_' + sequencer + '_assignment.txt'
    f_list = os.listdir(fold)
    all_res = {}
    for gt_file in f_list:
        accn = gt_file.replace(delstr,'')
        gt_file_results = get_gt_from_nidhi_results(os.path.join(fold,gt_file))
        for read_id in gt_file_results.keys():
            my_key = accn + '_' + read_id
            all_res[my_key] = gt_file_results[read_id] + (accn, read_id)
    return all_res



def get_file_to_assembly_map():
    '''
    At some point I gave a sequential name to all the individual files that the simulated reads
    were generated in (e.g. 'fq00001', 'fq00002', etc... This is just a shortcut to retrieve the
    mapping dictionary for those.
    :return:
    '''
    return get_dict_from_file(f2a)

def make_nhmmer_results_file(outfold, qs_id):
    '''
    Wrapper for read_and_combine_results() within the nhmmer dispatcher.

    Combs through a folder of raw HMMER outputs from a HIPPI run and generates the
    particular input set required to run hippi_nhmmer_dispatcher.read_and_combine_results().
    It then runs that function.
    :param outfold:
    :param qs_id:
    :return:
    '''
    res_l = [i for i in os.listdir(os.path.join(outfold,qs_id)) if i[-4:]=='.txt']
    res = list(map(lambda x: (1, os.path.join(outfold, qs_id, x)), res_l))
    read_and_combine_results(res, outfold, qs_id)

def read_blast_results(qsid, f2a=None, holdout=True):
    '''
    Reads a blast output and returns a dictionary with the results.

    :param qsid: same as qs_id elsewhere.
    :param f2a: Just avoids IO when reading a bunch of these.
    :return: dictionary object as:

        { '<accn>_<readID>' : ( accn, readID, cogID, matchingSequenceName, bitScore, eValue),
            ...
        }

    '''
    if holdout:
        blast_results_fold = blast_results_fold_ho
    else:
        blast_results_fold = blast_results_fold_fulldata

    fpa = os.path.join(blast_results_fold, qsid + '.txt')
    if f2a is None:
        f2a = get_file_to_assembly_map()
    f=open(fpa,'r')
    res={}
    for ln in f:
        a=ln.strip().split('\t')
        accn=f2a[a[0][:7]]
        read=a[0][8:]
        key=accn + '_' + read
        seq = a[1]
        cog = a[1].split('_')[1]
        res[key] = (accn, read, cog, seq, float(a[11]), float(a[10]))
    f.close()
    return res

def read_hippi_results(qsid, f2a=None, holdout=True):
    '''
    Reads one of the *combined* hippi results files that is located in the "all_results"
    subfolder of the HIPPI run. <qsid> is used because it expects the following directory
    structure:

        hippi_results/
            all_results/
                <qsid_1>.txt
                <qsid_2>.txt
                ...
            <qsid_1>/
                *_nhmmer_tbl.txt
            <qsid_2>/
                *_nhmmer_tbl.txt
            ...

    :param qsid:
    :param f2a:
    :return: dictionary object as:

        { '<accn>_<readID>' : ( accn, readID, cogID, bitScore, eValue, NumberOfHits),
            ...
        }

    '''
    if holdout:
        hippi_results_fold = hippi_results_fold_ho
    else:
        hippi_results_fold = hippi_results_fold_fulldata

    fpa = os.path.join(hippi_results_fold, qsid + '.txt')
    if f2a is None:
        f2a = get_file_to_assembly_map()
    f=open(fpa,'r')
    f.readline()
    res={}
    for ln in f:
        a=ln.strip().split('\t')
        accn=f2a[a[0][:7]]
        read=a[0][8:]
        key=accn + '_' + read
        cog=a[1].replace('hmmbuild_input_fam_','')[:7]
        if math.log(float(a[3]),10)<-1.0:
            res[key] = (accn, read, cog, float(a[2]), float(a[3]), int(a[5]))
    f.close()
    return res


if __name__=='__main__':
    all_hippi = {}
    for i in qs_ids_454:
        all_hippi.update(read_hippi_results(i))

    all_blast = {}
    for i in qs_ids_454:
        all_blast.update(read_blast_results(i))

    all_gt=get_full_gt_lookup_from_subfold('gt_454')


    # Here, since ground truth comes to us as a two-layered dictionary, it has to be
    #   converted into the "<accn>_<readID>" key format that the hippi/blast dicts use.

    # all_gt_keys=[]
    # k=list(all_gt.keys())  # list of accessions
    # all_gt_keys+=list(map(lambda x: k[0] + '_' + x, all_gt[k[0]].keys()))
    # len(all_gt_keys)
    # for i in k[1:]:
    #     all_gt_keys+=list(map(lambda x: i + '_' + x, all_gt[i].keys()))
    #