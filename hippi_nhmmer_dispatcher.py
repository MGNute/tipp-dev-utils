import os, multiprocessing, subprocess, shutil, datetime, sys, re
from functools import reduce

def get_all_hmms(hmm_repo):
    return list(map(lambda x: (x.replace('hmmbuild_fam_','')).replace('.hmm',''), os.listdir(hmm_repo)))

def run_nhmmer(query_seqs, hmm_file, out_fold, name):
    # invoc = '/projects/tallis/nute/code/hmmer-3.2.1/src/nhmmer --qhmm --tformat fasta --dna --noali -o %s --tblout %s %s %s'

    out_log = os.path.join(out_fold, name + '_nhmmer_out.log')
    out_tbl = os.path.join(out_fold, name + '_nhmmer_tbl.txt')
    # invoc_pop = invoc % (out_log, out_tbl, hmm_file, query_seqs )
    invoc_pop = ['/projects/tallis/nute/code/hmmer-3.2.1/src/nhmmer', '--qhmm', '--tformat', 'fasta', '--dna',
                 '--noali', '-o', out_log, '--tblout', out_tbl, hmm_file, query_seqs]

    sp = subprocess.Popen(invoc_pop, universal_newlines=True)
    rc = sp.wait()
    return (rc, out_tbl)

def run_all_nhmmers(qseqs, temp_folder, np, hmm_repo):

    qseqs_shm = os.path.join('/dev/shm', os.path.basename(qseqs))
    shutil.copyfile(qseqs, qseqs_shm)
    print('qseqs copied to: %s' % qseqs_shm)
    qseq_id = os.path.basename(qseqs).replace('.fasta','')
    out_folder = os.path.join(temp_folder, qseq_id)
    if not os.path.isdir(out_folder):
        os.mkdir(out_folder)

    hmms = get_all_hmms(hmm_repo)

    argslist=[]
    for i in hmms:
        a=None
        a=out_folder
        hmmfile = os.path.join(hmm_repo, 'hmmbuild_fam_%s.hmm' % i)
        argslist.append((qseqs_shm, hmmfile, a, i))

    p = multiprocessing.Pool(np)
    res=p.starmap(run_nhmmer, argslist)
    return res

def read_and_combine_results(res, outfold, qs_id):
    '''
    Here 'res' is the list of tuples returned by run_all_hmmers. In the initial setups, these three args are
    all part of what is typically the same folder, but separating them creates flexibility in case that is
    not the case. For example, the second entry of each element of "res" *should* be of the following form:

        <outfold>/<sequencer>_<crossfold>_<segment>/<hmm_model_name>_nhmmer_tbl_txt
                  |-------------------------------|
                        \equiv "qs_id"

    where the expression in the second subfolder name just spells out the contents of qs_id. It will then
    dump to an assumed existing folder "<outfold>/all_results"

    :param res: list object with each entry in the form of a tuple as:
                    ( return_code, hmmer_table_output_path ), ...
    :param outfold: Overall folder containing the output of an experiment.
                TODO: should verify that the subfolder "all_results" of "outfold" exists and is writable.
    :param qs_id: which set of query sequences the run was on. Usually
    :return:
    '''

    hmmer_tbl_col_defs = {'seqname':(0,37),
                            'targ accession':(37,48),
                            'hmmname':(48,81),
                            'quer accession':(81,91),
                            'hmmfrom':(91,99),
                            'hmm to':(99,107),
                            'alifrom':(107,115),
                            'ali to':(115,123),
                            'envfrom':(123,131),
                            'env to':(131,140),
                            'sq len':(140,147),
                            'strand':(147,155),
                            'E-value':(155,164),
                            'score':(164,171),
                            'bias':(171,177)}

    seq_results = {}

    for r in res:
        if r[0] != 1:
            print('%s returned code %s' % (r[1],r[0]))
        f = open(r[1],'r')
        # f = open(r, 'r')
        for ln in f:
            if ln[0] == '#':
                continue
            a = list(filter(lambda x: x !='',ln.split(' ')))
            try:
                if a[0] not in seq_results:
                    seq_results[a[0]] = [(a[2],float(a[13]), float(a[12]), float(a[14])),]  #(score, e-val, bias)
                else:
                    seq_results[a[0]].append((a[2], float(a[13]), float(a[12]), float(a[14]))) # (score, e-val, bias)
            except:
                print(ln.strip())
                print(a)
                sys.exit(0)
        f.close()

    final_f = open(os.path.join(outfold, 'all_results',qs_id + '.txt'),'w')
    final_f.write('seq_id\tbest_hmm\tscore\te-val\tbias\tnum_hits\n')
    for (k,v) in seq_results.items():
        if len(v)==1:
            max = v[0]
        else:
            max = reduce(lambda x, y: x if x[1]>y[1] else y, v)
        final_f.write('%s\t%s\t%s\t%s\t%s\t%s\n' % tuple([k,] + list(max) + [len(v),]))
    final_f.close()
    # return seq_results


def help():
    str_help = '''
hippi_nhmmer_dispatcher.py: Runs nhmmer for a group of hmms in parallel on a single fasta file. 

usage:  python3 hippi_nhmmer_dispatcher.py -n <# procs> -q <query_seqs> -temp <temp_folder> -pkg <pkg_fold>

'''

def get_args():
    nproc = int(sys.argv[sys.argv.index('-n')+1])
    print('running with %d processors' % nproc)
    qs = sys.argv[sys.argv.index('-q')+1]
    temp_fold = sys.argv[sys.argv.index('-temp')+1]
    pkg = sys.argv[sys.argv.index('-pkg') + 1]
    print("qs:\t%s" % qs)
    print("temp:\t%s" % temp_fold)
    print("pkg:\t%s" % pkg)
    if not os.path.isdir(temp_fold):
        print('Not Found: %s' % temp_fold)
        raise FileNotFoundError
    if not os.path.isdir(pkg):
        print('Not Found: %s' % pkg)
        raise FileNotFoundError
    if not os.path.isfile(qs):
        print('Not Found: %s' % qs)
        raise FileNotFoundError
    return (qs, temp_fold, nproc, pkg)

if __name__=='__main__':
    args=get_args()
    run_all_nhmmers(*args)


    # qseq_id = os.path.basename(args[0]).replace('.fasta', '')