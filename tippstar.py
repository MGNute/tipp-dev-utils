__author__ = 'Michael'

import sys
sys.path.append('/projects/tallis/nute/code/') #phylogeny_utilities')
# print sys.path

from phylogeny_utilities.utilities import *

testfile = '/projects/tallis/nute/data/ncbi/genomes/data/Yersinia_enterocolitica_palearctica_105_5R_r__uid63663/NC_015224.gbk'

import re, os, multiprocessing

# errsfile = open('/projects/tallis/nute/data/gi-errlines.txt','w')

def make_giant_cds_dna_fasta():
    prefix = '/projects/tallis/nute/data/ncbi-2017/'
    f_prefix = '/projects/tallis/nute/work/metagenomics/tippstar/2017/gbff_files/'
    filelist = open(prefix + 'file-list.txt','r')
    filelist_lookup = open(prefix + 'file-list-lookup.txt','w')
    filelist_dict = {}
    big_cds = open(prefix + 'ncbi_all_sequences.faa','w')
    big_dna = open(prefix + 'ncbi_all_sequences.fna','w')
    big_trerrs = open(prefix + 'ncbi-translation-errors.txt','w')
    big_lerrs = open(prefix + 'ncbi-length-errrors.txt','w')

    counter = 0

    files = get_list_from_file(prefix + 'file-list.txt')
    print "got file list..."
    filesfull = []
    for i in files:
        filesfull.append(f_prefix + i)
    print "created full list of files (filesfull)"

    p = multiprocessing.Pool(24)
    for i in range(83):
        mn = i * 1000
        if (i+1) * 1000 > len(filesfull):
            mx = len(filesfull)
        else:
            mx = (i+1) * 1000

        results = p.map(bp_genbank_get_CDS_dict,filesfull[mn:mx])
        print "got CDS dict for i=%s" % i

        for fi in results:
            counter += 1
            # if counter % 10000 == 0:
            #     print '\t%s files done' % counter
            record_file_results(big_cds, big_dna, big_lerrs, big_trerrs, counter, fi, filelist_dict, filelist_lookup)
        print '\tDone writing for i=%s'

    filelist_lookup.close()
    big_cds.close()
    big_dna.close()
    big_trerrs.close()
    big_lerrs.close()


def record_file_results(big_cds, big_dna, big_lerrs, big_trerrs, counter, fi, filelist_dict, filelist_lookup):
    cds = fi[0]
    dna = fi[1]
    sten = fi[2]
    trerrs = fi[3]
    lerrs = fi[4]
    fn = fi[5]
    # fn = prefix + fi.strip()
    filelist_lookup.write(str(counter) + '\t' + fn + '\n')
    # print str(counter) + '\t' + fi.strip()
    filelist_dict[counter] = fn
    # cds, dna, sten, trerrs, lerrs = bp_genbank_get_CDS_dict(fn)
    for locus in cds.keys():
        big_cds.write('>' + str(counter) + '_' + locus + '\n')
        big_cds.write(cds[locus] + '\n')
    for locus in dna.keys():
        big_dna.write('>' + str(counter) + '_' + locus + '\n')
        big_dna.write(str(dna[locus]) + '\n')
    for trerr in trerrs:
        big_trerrs.write(str(counter) + '_' + trerr + '\n')
    for lerr in lerrs:
        big_lerrs.write(str(counter) + '_' + lerr + '\n')

        # if counter>3:
        #     break


def get_file_taxid_map():
    sourcefile='/projects/tallis/nute/data/ncbi/folder-taxid-map.txt'
    strRegEx = '^(?P<filename>.*?):.*taxon:(?P<tid>\d+)\"'
    pat = re.compile(strRegEx)
    # if testst != None:
    #     result = pat.search(testst)
    #     print result.group('filename')
    #     print result.group('tid')

    fi = open(sourcefile,'r')
    outdict = {}
    for i in fi:
        if len(i)>1:
            result = pat.search(i)
            try:
                f = '/projects/tallis/nute/data/ncbi/' + result.group('filename')
                tid = result.group('tid')
                outdict[f]=tid
            except:
                print i

    fi.close()
    return outdict

def get_file_index_map():
    filename = '/projects/tallis/nute/data/ncbi/file-list-lookup.txt'
    myf = open(filename,'r')
    outdict = {}
    for i in myf:
        if len(i)>1:
            a=i.strip().split('\t')
            id = int(a[0])
            fn = a[1]
            outdict[id]=fn
    myf.close()
    return outdict

def make_species_mapping():
    sequence_list='/projects/tallis/nute/data/ncbi/ncbi_sequence_names.txt'
    species_mapping_location = '/projects/tallis/nute/data/ncbi/species.mapping'
    seqnames = get_list_from_file(sequence_list)
    indexToFile = get_file_index_map()
    fileToTaxid = get_file_taxid_map()
    species_mapping = open(species_mapping_location,'w')
    species_mapping.write('seqname,tax_id\n')
    k=0

    for i in seqnames:
        seq,num = file_from_seqname(i.strip())
        taxid = fileToTaxid[indexToFile[num]]
        species_mapping.write(seq + ',' + taxid + '\n')
        # k+=1
        # if k>1000:
        #     break

    species_mapping.close()

def read_species_mapping_to_dict():
    fi = open('/projects/tallis/nute/data/ncbi/species.mapping','r')
    sm = {}
    for i in fi:
        if len(i)>1:
            a=i.strip().split(',')
            sm[a[0]]=i.strip()
    fi.close()
    return sm

def make_species_mapping_per_cog():
    prefix = '/projects/tallis/nute/data/ncbi/motu_all/'
    outpref = '/projects/tallis/nute/data/tipp/tipp/refpkg_new/'
    cog_list = get_list_from_file(prefix + 'cog-list.txt')
    sm = read_species_mapping_to_dict()

    for i in cog_list:
        print i
        outfile = open(outpref + i + '/species.mapping','w')
        taxa = get_list_from_file(prefix + i + '-taxa.txt')
        for j in taxa:
            outfile.write(sm[j]+'\n')
        outfile.close()
        del taxa

def file_from_seqname(strSeqName):
    seq = strSeqName[1:]
    file_num = int(seq[0:seq.index('|')])
    return seq, file_num


def get_gi(str_test):
    strRegex='data/(?P<location>.*?).gbk:VERSION\s+(?P<ver>\S+)\s+GI:(?P<gi>\S+)'
    pat = re.compile(strRegex)
    result = pat.search(str_test)
    try:
        l= result.group('location')
        v= result.group('ver')
        g= result.group('gi')
        return (l,v,g)
    except:
        errsfile.write(str_test)
        print "ERROR: %s" %str_test[0:100]
        return (None, None, None)


def replace_sequence_header_perl(loc,ver,gi):
    fileloc = '/projects/tallis/nute/data/krakendb/library/Bacteria/' + loc + '.fna'
    try:
        assert os.path.exists(fileloc)==True
    except:
        print "%s is not a valid file" % loc
        return None

    # fasta=open(fileloc,'r')
    # header = fasta.readline()
    # newheader = header.replace(ver,ver + '|')[1:]
    args = (ver, gi, ver, fileloc)
    line = 'perl -pi -e \'s/%s/gi|%s|ref|%s|/g\' %s' %args
    return line

def split_files_by_mode():
    split_locs={'COG0012':1120,
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

def rename_fastaheaders_in_krakendb():
    '''
    Function to generate a script required to rename the fasta files created by
    Biopython from the genbank files. The fasta sequence headers must
    conform to the specs required by kraken.

    :return:
    '''
    mapsfile = '/projects/tallis/nute/data/gi-maps.txt'
    scriptfile = '/projects/tallis/nute/data/rename-krakendb-sequences.sh'
    mf = open(mapsfile,'r')
    sf = open(scriptfile,'w')

    for ln in mf:
        (i,j,k) = get_gi(ln)
        if i <> None:
            cmd = replace_sequence_header_perl(i,j,k)
            sf.write(cmd + '\n')

    mf.close()
    sf.close()

# seq = enumerate(SeqIO.parse(testfile,"genbank"))
# ind, rec = seq.next()

if __name__=='__main__':
    pass
    # (lo, ve, gi) = get_gi(teststr2)
    # print replace_sequence_header_perl(lo, ve, gi)

    # rename_fastaheaders_in_kradendb()
    make_giant_cds_dna_fasta()
    # make_species_mapping()
    # make_species_mapping_per_cog()