__author__ = 'Michael'

prefix = 'C:/Users/Michael/Grad School Stuff/Research/Phylogenetics/Metagenomics/tippstar/'

def make_pbs_files(cog):
    pbstmp = open(prefix + 'pbs/COGtemp.pbs','rU')
    pbs = pbstmp.read()
    pbs = pbs.replace('mycog',cog)
    newpbs = open(prefix + 'pbs/' + cog + '_pasta.pbs','w')
    newpbs.write(pbs)
    newpbs.close()
    # myshell = open(prefix + 'submit_pastas.sh','a')
    # myshell.write('qsub ' + cog + '_pasta.pbs\n')
    # myshell.close()

def main():
    cogs = ['COG0012','COG0016','COG0018','COG0048','COG0049','COG0052',
            'COG0080','COG0081','COG0085','COG0087','COG0088','COG0090',
            'COG0091','COG0092','COG0093','COG0094','COG0096','COG0097',
            'COG0098','COG0099','COG0100','COG0102','COG0103','COG0124',
            'COG0172','COG0184','COG0185','COG0186','COG0197','COG0200',
            'COG0201','COG0202','COG0215','COG0256','COG0495','COG0522',
            'COG0525','COG0533','COG0541','COG0552']
    for i in cogs:
        make_pbs_files(i)

if __name__=='__main__':
    main()