#!/usr/bin/env python3

'''
Title: blastRecip.py

Date: 01.29.2019

Author: Mariam Miari

List of functions and methods: open-split-print-startswith-rstrip-frozenset

Run by : The program is to be run this way: ./blastRecip.py concproteome.blastp output2

Input: concProteome.blastp

Output: see "output2". 4816 lines were obtained.

Descrption and Procedure : Queries and targets are extracted from the blast
output and added to a dictionary (myDict) as keys and values, respectively.Each
of queries and targets were separated and later split into 2 columns.If the query
of one species was the target of another, and vice-versa, then they are orthologous
thus they were printed (check lines 49-53). To make sure that the entries are
unique, frozenset function was used where by each set will contain other subsets.
'''
myDict={}
set1 = set()
with open("concProteome.blastp") as recip:
    for line in recip:
        #Save Query in the variable "query"
        if line.startswith("Query="):
            query= line.rstrip().split(" ")[1]
            #for the purpose of obtaining the queries and targets as two species, each in a column.
            query_prefix = query.split("_")[0]
            ignore = False # define ignore variable for later use.
            #Save Target in the variable "target"
        elif line.startswith(">"):
            target= line.rstrip()[1:]
            #for the purpose of obtaining the queries and targets as two species, each in a column.
            target_prefix = target.split("_")[0]
            #Ignore the 100% similarity sequences
            if query == target or query_prefix == target_prefix:
                continue
            elif not ignore:
                #save only the best hit
                ignore = True
                besthit = target
                #Creating a dictionary with query as keys and besthits as values.
                myDict[query] = besthit
with open("output2", "w") as out:
    #loop over the keys (queries) to find the orthologous sequences.
    for query in myDict:
        besthit = myDict[query]#
        if besthit in myDict:
            #to get the reciprocals.
            if query == myDict[besthit]:
                #pick only the unique pairs
                if frozenset({query,besthit}) in set1:
                    continue
                #if query.startswith("Pax") and besthit.startswith("Phch"): #Another way of doing it.
                    #continue
                else:
                    print(query, besthit, sep ="\t", file=out)
                    set1.add(frozenset({query,besthit}))
