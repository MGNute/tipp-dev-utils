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