
'''collapseLowSupportBranches.py verion 1.0 copyright (C) Joel Sharbrough 2020

USAGE:

	python collapseLowSupportBranches.py treeFile supportNum

Notes:
	
	-Treefile should be in newick format
	-supportNum should be an integer [50]
	
	
'''

import dendropy
import sys
import random
def collapseLowSupportBranches(treeFile,supportNum=50):
    tree = dendropy.Tree.get(path=treeFile,schema='newick')
    for e in tree.postorder_edge_iter():
        if e.head_node.label != None:
            if int(e.head_node.label) < supportNum:
                e.collapse()
    fileSplit = treeFile.split('.')
    outfile = open(fileSplit[0] + '.collapsed.trees','w')
    outfile.write(tree.as_string(schema='newick'))
    outfile.close()

helpStatement = '\ncollapseLowSupportBranches.py verion 1.0 copyright (C) Joel Sharbrough 2020\n\nUSAGE:\n\n\tpython collapseLowSupportBranches.py treeFile supportNum\n\nNotes:\n\n\t-Treefile should be in newick format\n\t-supportNum should be an integer [50]\n\n'
if len(sys.argv) == 1:
	sys.stderr.write(helpStatement)
elif sys.argv[1] == 'help' or sys.argv[1] == 'h':
	sys.stderr.write(helpStatement)
elif len(sys.argv) == 3:
	collapseLowSupportBranches(sys.argv[1],int(sys.argv[2]))
elif len(sys.argv) == 2:
	collapseLowSupportBranches(sys.argv[1])
else:
	sys.stderr.write(helpStatement)