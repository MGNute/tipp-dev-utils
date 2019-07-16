'''
This was written in an effort to make a TIPP reference out of the SILVA data, unfortunately
    that ended in failure because mapping from SILVA to NCBI Taxonomy is not well-defined.
'''
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
import re, json, sys, os
import datetime
from multiprocessing import Pool
import subprocess
from ncbi_querying import *

from phylogeny_utilities.utilities import *

file_lookup = get_dict_from_tab_delimited_file('/projects/tallis/nute/data/ncbi-2017/file-list-lookup.txt')
strRE = 'derived from\s+(?P<insdn>[A-Za-z0-9]+)\.'
pat = re.compile(strRE)
file_ids = get_list_from_file('/projects/tallis/nute/data/ncbi-2017/all_cog_uniq_file_ids.txt')

def get_file_CDS_to_INSDN_dict(id):
    file_path = file_lookup[id]

    si = SeqIO.parse(file_path,'genbank')

    cds_to_insdn = {}

    for sr in si:
        res = pat.search(sr.annotations['comment'])
        if res is not None:
            acc = res.group('insdn')
        else:
            acc = None
        for feat in sr.features:
            if feat.type=='CDS' and 'protein_id' in feat.qualifiers.keys():
                k = id + '_pt:' + feat.qualifiers['protein_id'][0]
                cds_to_insdn[k] = acc

    return cds_to_insdn

def make_cog_silva_taxonid_map(cog):
    wd = '/projects/tallis/nute/data/ncbi-2017'
    cog_id = cog.replace('.txt','')
    silva_tax = os.path.join(wd,'silva_tax','tax_slv_ssu_128.acc_taxid_dedupe_tabd')
    std = get_dict_from_tab_delimited_file(silva_tax)
    cog_insdn = os.path.join(wd,'COG_insdn_maps',cog)
    ci = get_dict_from_tab_delimited_file(cog_insdn)

    outmap = open(os.path.join(wd,'COG_silva_maps',cog),'w')
    for k in ci.keys():
        taxon = None
        try:
            taxon = std[ci[k]]
        except:
            taxon = None
        outmap.write('%s,%s\n' % (k,taxon))
    outmap.close()





# def get_ncbi_assembly_id(aid):
#     """
#     Writen by Nute on 4/10/18.
#     Takes NCBI Accession ID and returns an NCBI Assembly Accession if one exists
#
#     Parameters
#     -----------
#     aid : string
#           NCBI Accession ID
#
#     Returns
#     --------
#     assemb_id : string
#           NCBI Assembly Accession ID
#     """
#     p = subprocess.Popen(["curl", "-s",
#                           "https://www.ncbi.nlm.nih.gov/nuccore/" + aid],
#                          stdout=subprocess.PIPE)
#     (content, error) = p.communicate()
#     uidtext = "/assembly\?LinkName=nuccore_assembly&amp;from_uid=(?P<uid>[0-9]+)"
#     assembtext = "<p class=\"title\"><a href=\"/assembly/(?P<assid>[A-Za-z0-9\._]+)\""
#     uidre = re.compile(uidtext)
#     uidmatch = uidre.search(content)
#     assembre = re.compile(assembtext)
#     assembmatch = assembre.search(content)
#     if assembmatch is not None and assembmatch('asssid'):
#         return assembmatch.group('assid')
#     elif uidmatch is not None and uidmatch.group('uid'):
#         uid = uidmatch.group('uid')
#         p = subprocess.Popen(["curl", "-s",
#                               "https://www.ncbi.nlm.nih.gov/assembly?LinkName=nuccore_assembly&from_uid=" + uid],
#                              stdout=subprocess.PIPE)
#         (content2, error) = p.communicate()
#         assembmatch2 = assembre.search(content2)
#         if assembmatch2 is not None and assembmatch2.group('assid'):
#             return assembmatch2.group('assid')
#     return None

def get_file_CDS_to_ncbi_taxid_dict(id):
    file_path = file_lookup[id]
    pa, fi = os.path.split(file_path)

    si = SeqIO.parse(file_path,'genbank')

    cds_to_ncbi = {}

    for sr in si:
        acc = None
        if sr.features[0].type != 'source':
            print ("%s - record: %s does not have a source" % (fi, sr.id))
            continue
        elif 'db_xref' not in sr.features[0].qualifiers.keys():
            print( "%s - record: %s does not have field called db_xref" % (fi, sr.id))
            continue
        else:
            for dx in sr.features[0].qualifiers['db_xref']:
                if dx[0:5]=='taxon':
                    acc = dx.replace('taxon:','')
                    continue
            if acc is None:
                print ("%s - record: %s does not have a db_xref entry for taxon" % (fi, sr.id))
                continue


        for feat in sr.features:
            if feat.type=='CDS' and 'protein_id' in feat.qualifiers.keys():
                k = id + '_pt:' + feat.qualifiers['protein_id'][0]
                cds_to_ncbi[k] = (acc,sr.name)

    return cds_to_ncbi

def get_all_files_in_list(file_list):
    my_ids = get_list_from_file(file_list + '.txt')
    all_dict = {}
    # ct = 0
    for i in my_ids:
        # ct +=1
        # all_dict[i]=get_file_CDS_to_INSDN_dict(i)
        all_dict[i] = get_file_CDS_to_ncbi_taxid_dict(i)
        # if ct % 25 == 0:
        #     print "%s done" % ct

    # outf = open(file_list + '.json','w')
    outf = open(file_list + '_ncbi.json', 'w')
    json.dump(all_dict,outf)
    outf.close()
    print ('%s done' % file_list)

# def get_full_cds_dict():
#     str_name = '/projects/tallis/nute/data/ncbi-2017/partial_cog_uniq/all_cog_uniq_file_ids_'
#
#     full = {}
#
#     for i in range(320):
#         print str(i)
#         # myf = open(str_name + str(i) + '.json', 'r')
#         myf = open(str_name + str(i) + '_ncbi.json','r')
#         a = json.load(myf)
#         myf.close()
#         full.update(a)
#         del a
#
#     # js = open('/projects/tallis/nute/data/ncbi-2017/full_cds_dict.json','w')
#     js = open('/projects/tallis/nute/data/ncbi-2017/full_cds_dict_ncbi.json', 'w')
#     json.dump(full,js)
#     js.close()
#     return full

def get_full_cds_dict():
    # pa = '/projects/tallis/nute/data/ncbi-2017/full_cds_dict.json'
    pa = '/projects/tallis/nute/data/ncbi-2017/full_cds_dict_ncbi.json'
    myf = open(pa,'r')
    full = json.load(myf)
    myf.close()
    # print "loaded full_cds_dict.json"
    print ("loaded full_cds_dict_ncbi.json")
    return full

def make_cog_insdn_map(cog_file, insdn_lookup):
    fas = read_from_fasta(cog_file)
    out_file = cog_file.replace('.faa','.txt')
    out_file = out_file.replace('COGS','COG_insdn_maps')
    of = open(out_file, 'w')

    for i in fas.keys():
        a = i.strip().split('_')
        fi_ind = a[0]
        try:
            insdn = insdn_lookup[fi_ind][i]
        except:
            insdn = '(Not Found)'
        of.write('%s\t%s\n' % (i,insdn))
    of.close()

def make_cog_ncbi_map(cog_file, ncbi_lookup):
    fas = read_from_fasta(cog_file)
    out_file = cog_file.replace('.faa','.txt')
    out_file = out_file.replace('COGS','COG_ncbi_maps')
    of = open(out_file, 'w')

    for i in fas.keys():
        a = i.strip().split('_')
        fi_ind = a[0]
        try:
            ncbi_acc = ncbi_lookup[fi_ind][i][1]
            taxon = ncbi_lookup[fi_ind][i][0]
        except:
            ncbi_acc = '(Not Found)'
            taxon = '(Not Found)'
        # of.write('%s,%s,%s\n' % (i,ncbi_acc,taxon))
        of.write('%s\t%s\n' % (i,taxon))
    of.close()

def make_file_seq_accn_taxid_lookup():
    '''
    Added by Nute on 4/10/2018 to figure out which accessions map to which taxa
    :return:
    '''
    ncbi_lkp = get_full_cds_dict()
    outf = open('/projects/tallis/nute/data/ncbi-2017/file_seq_accn_taxid.txt','w')
    outf.write('file_id\tseqname\tncbi_tax_id\tncbi_accn_id\n')

    fct = 0
    lct = 0
    for fid in ncbi_lkp.keys():
        for seqid in ncbi_lkp[fid].keys():
            tax_accn = ncbi_lkp[fid][seqid]
            outf.write('%s\t%s\t%s\t%s\n' % (fid, seqid, tax_accn[0], tax_accn[1]))
            lct +=1
        fct +=1
        if fct % 1000 == 0:
            print (fct)
    outf.close()
    print ('wrote %s lines' % lct)

def make_accn_taxid_lookup():
    '''
    Added by Nute on 4/10/2018 to figure out which accessions map to which taxa
    :return:
    '''
    ncbi_lkp = get_full_cds_dict()
    outf = open('/projects/tallis/nute/data/ncbi-2017/accn_taxid_lookup_2.txt','w')
    outf.write('ncbi_accession_id\tnum_taxon_ids\tncbi_taxon_id\n')
    accns = {}
    fct = 0
    for fid in ncbi_lkp.keys():
        for seqid in ncbi_lkp[fid].keys():
            tax_accn = ncbi_lkp[fid][seqid]
            if tax_accn[1] not in accns.keys():
                accns[tax_accn[1]] = set([tax_accn[0]])
            else:
                accns[tax_accn[1]].add(tax_accn[0])
        fct +=1
        if fct % 1000 ==0:
            print ('file ct:\t%s' % fct)
    print ('final step ....')
    accnct = 0
    for i in accns.keys():
        outf.write('%s,%s,%s\n' % (i, len(foo), ','.join(list(accns[i]))))
        accnct+=1
        if accnct % 20000 == 0:
            print ('accn ct:\t%s' % accnct)
    outf.close()

if __name__=='__main__':
    # get_all_files()
    # print datetime.datetime.now()
    # ncbi_lkp = get_full_cds_dict()
    # print datetime.datetime.now()
    #
    # for line in sys.stdin:
    #     print line.strip()
        # make_cog_insdn_map(line.strip(), insdn_lkp)
        # make_cog_ncbi_map(line.strip(), ncbi_lkp)
        # make_cog_silva_taxonid_map(line.strip())


    # p = Pool(16)
    # l = []
    # for i in range(320):
    #     l.append('/projects/tallis/nute/data/ncbi-2017/partial_cog_uniq/all_cog_uniq_file_ids_' + str(i))
    # p.map(get_all_files_in_list,l)

    # make_file_seq_accn_taxid_lookup()
    make_accn_taxid_lookup()