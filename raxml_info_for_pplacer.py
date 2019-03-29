'''
Creates a RAxML_info file in a format that pplacer can digest. Unfortunately giving
    pplacer a file like this is deprecated but maintained for compatibility, but the
    RAxML_info file structure has changed and must be translated to the old form to
    use in TIPP.

    This basically just takes an old output file as a dummy and does some regex work
    to replace the relevant spots with corresponding data from the new file. That
    means that the result is basically a fake RAxML_info file with real data, so it
    shouldn't be interpreted beyond its intended use.
'''
import re, os, sys, shutil

model_fi=os.path.join(os.path.split(os.path.abspath(__file__))[0],'raxml_info_model.txt')


def make_empty_subvals():
    subvals = {'ALNPATTERNS' : None,
            'PERCENTGAPS' : None,
            'RAXMLCALLEDAS' : None,
            'INFERENCETIME' : None,
            'CATLIKELIHOOD' : None,
            'EXECINFOWRITTEN' : None,
            'EXECUTIONTIME' : None,
            'PIA' : None,
            'PIC' : None,
            'PIG' : None,
            'PIT' : None,
            'ALPHA' : None,
            'ACGTR' : None,
            'AGGTR' : None,
            'ATGTR' : None,
            'CGGTR' : None,
            'CTGTR' : None,
            'GTGTR' : None}
    return subvals

def populate_subvals(raxml_info_optimized_file,st=False):
    regex_gtr = ''
    regex_gtr = regex_gtr + 'alpha: (?P<ALPHA>\d+\.\d+[E\-*\d]*)\n'
    regex_gtr = regex_gtr + 'Tree-Length: \d+\.\d+[E\-*\d]\n'
    regex_gtr = regex_gtr + 'rate A <-> C: (?P<ACGTR>\d+\.\d+[E\-*\d]*)\n'
    regex_gtr = regex_gtr + 'rate A <-> G: (?P<AGGTR>\d+\.\d+[E\-*\d]*)\n'
    regex_gtr = regex_gtr + 'rate A <-> T: (?P<ATGTR>\d+\.\d+[E\-*\d]*)\n'
    regex_gtr = regex_gtr + 'rate C <-> G: (?P<CGGTR>\d+\.\d+[E\-*\d]*)\n'
    regex_gtr = regex_gtr + 'rate C <-> T: (?P<CTGTR>\d+\.\d+[E\-*\d]*)\n'
    regex_gtr = regex_gtr + 'rate G <-> T: (?P<GTGTR>\d+\.\d+[E\-*\d]*)\n'
    regex_pi = ''
    regex_pi = regex_pi + 'freq pi\(A\): (?P<PIA>\d+\.\d+[E\-*\d]*)\n'
    regex_pi = regex_pi + 'freq pi\(C\): (?P<PIC>\d+\.\d+[E\-*\d]*)\n'
    regex_pi = regex_pi + 'freq pi\(G\): (?P<PIG>\d+\.\d+[E\-*\d]*)\n'
    regex_pi = regex_pi + 'freq pi\(T\): (?P<PIT>\d+\.\d+[E\-*\d]*)\n'

    subvals = make_empty_subvals()

    if not st:
        test_opt_f = open(raxml_info_optimized_file,'r')
        test_opt_str = test_opt_f.read()
        test_opt_f.close()
    else:
        test_opt_str = raxml_info_optimized_file


    # GTR parameters
    gtr_re = re.compile(regex_gtr)
    gtrvals = gtr_re.search(test_opt_str)
    pi_re = re.compile(regex_pi)
    pivals = pi_re.search(test_opt_str)

    if gtrvals is None or pivals is None:
        if gtrvals is None:
            print ('no match found for the gtr regex')
        if pivals is None:
            print ('no match found for the pivals')
        sys.exit(0)

    subvals['PIA'] = pivals.group('PIA')
    subvals['PIC'] = pivals.group('PIC')
    subvals['PIG'] = pivals.group('PIG')
    subvals['PIT'] = pivals.group('PIT')
    subvals['ALPHA'] = gtrvals.group('ALPHA')
    subvals['ACGTR'] = gtrvals.group('ACGTR')
    subvals['AGGTR'] = gtrvals.group('AGGTR')
    subvals['ATGTR'] = gtrvals.group('ATGTR')
    subvals['CGGTR'] = gtrvals.group('CGGTR')
    subvals['CTGTR'] = gtrvals.group('CTGTR')
    subvals['GTGTR'] = gtrvals.group('GTGTR')

    aln_re_str = 'Alignment Patterns: (?P<ALNPATTERNS>\d+)'
    aln_re = re.compile(aln_re_str)
    aln = aln_re.search(test_opt_str)
    try:
        subvals['ALNPATTERNS'] = aln.group('ALNPATTERNS')
    except:
        print ('Number of alignment patterns not found')

    pctgaps_str = 'Proportion of gaps and completely undetermined characters in this alignment: (?P<PERCENTGAPS>\d+\.\d+)'
    pctgaps_re = re.compile(pctgaps_str)
    pctgaps = pctgaps_re.search(test_opt_str)
    try:
        subvals['PERCENTGAPS'] = pctgaps.group('PERCENTGAPS')
    except:
        print ('Percentage gaps not found')

    calledas_str = 'RAxML was called as follows:\n.*\n(?P<RAXMLCALLEDAS>.*)\n'
    calledas_re = re.compile(calledas_str)
    calledas = calledas_re.search(test_opt_str)
    try:
        subvals['RAXMLCALLEDAS'] = calledas.group('RAXMLCALLEDAS')
    except:
        print ('RAxML called-as string not found')

    execfile_str = 'Execution Log File written to:\s*(?P<EXECINFOWRITTEN>.*)\n'
    execfile_re = re.compile(execfile_str)
    execfile = execfile_re.search(test_opt_str)
    try:
        subvals['EXECINFOWRITTEN'] = execfile.group('EXECINFOWRITTEN')
    except:
        print ('Execution info written to not found')

    return subvals

def make_new_raxml_info_file(outpath, modelpath, subvals):
    outfile = open(outpath,'w')

    modelfi = open(modelpath,'r')
    modelstr = modelfi.read()
    modelfi.close()

    for k in subvals.keys():
        if subvals[k] is not None:
            modelstr = modelstr.replace('xxx' + k + 'xxx', subvals[k])

    outfile.write(modelstr)
    outfile.close()

def print_help():
    print ('''
This script makes a new RAxML_info file that is suitable for use with pplacer from another RAxML_info file
    in a differnt format. Currently it is set up to convert files created by running raxml 8.2.9 with the
    '-f -e ' model to optimize branch lengths and model parameters.

    usage: python raxml_info_for_pplacer.py <raw_raxml_file> <destination_file>
    
        <destionation_file> is optional. Omitting it will cause a backup called <raw_raxml_file>.bkp to
            be generated in the same folder as the original, and then the <raw_raxml_file> to be 
            replaced.
        (Note that a previous version of this file also required the model file as input, but now it is provided
        as part of the tipp_dev_utils repo.)
''')

def backup_old_file(pth):
    fo, fi = os.path.split(pth)
    bkp = fi + '.bkp'
    shutil.copy(pth,os.path.join(fo,bkp))


if __name__=='__main__':
    if '-h' in sys.argv or '--help' in sys.argv:
        print_help()
        sys.exit(0)
    else:
        old_fi = sys.argv[1]
        # model_fi = sys.argv[2]
        if len(sys.argv)>2:
            new_fi = sys.argv[2]
        else:
            print('NOTE: no second argument was given, so we will make these changes in place and save a '
                  'backup in the same folder as the original.')
            backup_old_file(old_fi)
            new_fi = old_fi
        subvals = populate_subvals(old_fi)
        make_new_raxml_info_file(new_fi,model_fi,subvals)
        sys.exit(0)