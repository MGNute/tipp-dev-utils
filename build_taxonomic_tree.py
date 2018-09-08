#!/lusr/bin/python
'''
Reads in a taxonomy file, outputs a taxonomic tree in newick format
Created on Jan 31, 2012

python2.7 build_taxonomic_tree.py ../data/taxonomy/all_taxon.taxonomy taxonomic.tree
python2.7 build_taxonomic_tree.py ../data/taxonomy/all_taxon.taxonomy taxonomic.tree ../data/taxonomy/species.mapping

@author: namphuon
'''
from dendropy import Tree, Node
import sys
import os


if __name__ == '__main__':
  taxonomyFile = sys.argv[1]
  speciesList = sys.argv[2]
  taxonomyTree = sys.argv[3]

  species = {}
  lines = open(speciesList,'r')
  for line in lines:
    species[line.strip()] = line.strip()
    
  lines = open(taxonomyFile,'r')
  header = lines.next()
  nodes_dict = {}

  #Read first line, root node
  line = lines.next()
  results = line.strip().split(',')
  tree = Tree()  
  root = Node()
  root.__dict__['label'] = results[0].replace("\"","")
  nodes_dict[results[0].replace("\"","")] = root
  
  prune = ['1']
  
  #Add root node to tree
  tree.__dict__['seed_node'].add_child(root)
  for line in lines:
    results = line.strip().split(',')
    node = Node();
    node.__dict__['label'] = results[0].replace("\"","")
    nodes_dict[results[0].replace("\"","")] = node
    nodes_dict[results[1].replace("\"","")].add_child(node)
    if results[0].replace("\"","") not in species:
      prune.append(results[0].replace("\"",""))
      
  for taxa in prune:
    nodes_dict[taxa].label=''
  tree.delete_outdegree_one_nodes()
      
  output = open(taxonomyTree, 'w')
  output.write(str(tree) + ";");
  
  output.close()
  lines.close()