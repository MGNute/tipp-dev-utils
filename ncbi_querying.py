'''
Getting IDs and stuff from NCBI.
'''
from Bio import Entrez
Entrez.email = 'nute2@illinois.edu'
import re

def get_ncbi_assembly_id(aid):
    handle2 = Entrez.efetch(db="nucleotide", id=aid, rettype="gb", retmode="text")
    return_str = handle2.read()
    re_str = 'Assembly: (?P<assemb>[A-Za-z0-9\._]+)'
    assemb_re = re.compile(re_str)
    assid = assemb_re.search(return_str)
    if assid is not None and assid.group('assemb'):
        return assid.group('assemb')
    else:
        return None

"""
The part below contains functions for getting NCBI taxonomy information from
NCBI accession IDs.

Written by Erin K. Molloy (emolloy2@illinois.edu) in September 2016.
"""

def get_ncbi_taxon_id(aid):
    """
    (Written by Erin)
    Takes NCBI Accession ID and returns NCBI Taxon ID.

    Parameters
    -----------
    aid : string
          NCBI Accession ID

    Returns
    --------
    xid : string
          NCBI Accession ID
    """
    import subprocess
    pad = 100

    p = subprocess.Popen(["curl", "-s",
                          "https://www.ncbi.nlm.nih.gov/nuccore/" + aid],
                         stdout=subprocess.PIPE)
    (content, error) = p.communicate()

    start = "ORGANISM="
    end = "&amp;"
    s = content.find(start) + len(start)
    e = s + content[s:s+pad].find(end)
    return content[s:e]


def get_ncbi_taxonomy(aids, ofile):
    """
    Takes NCBI Accession ID and returns NCBI taxonomic lineage.

    Parameters
    ----------
    xid : list of strings
          NCBI Accession IDs


    Returns
    -------
    df : dataframe
         NCBI taxonomic lineage
    """
    from ete3 import NCBITaxa
    import pandas as pd
    of = open(ofile, 'w')

    ncbi = NCBITaxa()

    # keys = ["accid", "taxid",
    # "superkingdom", "phylum", "class",
    # "order", "family", "genus", "species"]
    keys = ['accid', 'taxid']

    rows = []
    # print "starting the list of accession ids"
    ct = 0
    for aid in aids:
        xid = get_ncbi_taxon_id(aid)

        # lineage = ncbi.get_lineage(xid)
        # ranks = ncbi.get_rank(lineage)
        # names = ncbi.get_taxid_translator(lineage)

        tax = {}
        tax["accid"] = aid
        tax["taxid"] = xid
        of.write('%s,%s\n' % (aid, xid))
        # for k in keys[2:]:
        # try:
        # i = ranks.values().index(k)
        # except ValueError:
        # i = -1

        # if i == -1:
        # tax[k] = "NA"
        # else:
        # tax[k] = names[ranks.keys()[i]].replace(' ', '_')
        # rows.append(tax)
        # if ct % 100 ==0 :
        # ct +=1
        # print "%s done." % ct

    of.close()
    # df = pd.DataFrame(rows, columns=keys)
    # return df


if __name__ == "__main__":
    import sys

    ifile = str(sys.argv[1])
    ofile = str(sys.argv[2])

    with open(ifile, 'r') as f:
        aids = f.readlines()
        aids = [aid.rstrip('\n') for aid in aids]

    get_ncbi_taxonomy(aids, ofile)
    # df = get_ncbi_taxonomy(aids)
    # df.to_csv(ofile, header=True, index=False)

