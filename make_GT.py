from phylogeny_utilities.utilities import *
import os
cogs=[i for i in os.listdir('.') if i[0:3]=='COG']
ranks = ['tax_name','superkingdom','phylum','class','order','family','genus','species']
lnstr = '\t'.join(['%s',] * 10) + '\n'

def add_cog_to_groundtruth_file(cog,gtfile):

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
    print 'wrote %s lines for cog %s' % (ct,cog)


