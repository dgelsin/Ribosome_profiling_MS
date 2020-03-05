#!/usr/bin/env python


import sys


f=sys.argv[1]

g=sys.argv[2]

#h=sys.argv[3]

keywords = set()
with open(f) as list_file:
    for line in list_file:
        if line.strip():
            keywords.add(line.strip())

#with open(g) as master_file:
#	for line in master_file:
#		if set(line.split()) & keywords:
#			print line.rstrip()

#original
with open(g) as master_file:
	for line in master_file:
		if set(line.split()) & keywords:
			print line.rstrip()

#			search_results.write(line)
#    with open(h, 'w') as search_results:
        