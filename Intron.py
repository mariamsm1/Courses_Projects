#Finding the intron size using the exon start and end columns in a gtf file.
-----------------------------
#!/usr/bin/env python3

myDict = {}
myIntron = {}
with open("genemark.gtf", "r") as scaf:
	for line in scaf:
		line=line.rstrip().split("\t")
		if line[2] == "exon":
			scaffold = line[0]
			gene_id = line[8].split(" ")[1].replace('\"', '').replace(';', '')
			start=line[3]
			end = line[4]
			key = scaffold+','+gene_id
			coord = [start,end]
			if key in myDict.keys():
				myDict[key].append(coord)
			else:
				myDict[key] = [coord]
print(myDict)
for gene in myDict.keys():
	for pos1,pos2 in zip(myDict[gene],myDict[gene][1:]): # we write it twice because we need the elements of both lists in the keys. [1:] means the second list in the big list
        Intronstart,Intronend = pos1[1],pos2[0] #so here i take the end of first list and start of second list in the zipped positions.
        intron_coord = [Intronstart,Intronend]
        if gene in myIntron.keys():
            myIntron[gene].append(intron_coord)
        else:
            myIntron[gene]= [intron_coord]
intronsum = 0
intron_num = 0
for gene in myIntron:
    for coord in myIntron[gene]:
        print(gene, coord[0],coord[1], int(coord[1])-int(coord[0]))
        intronsum += int(coord[1])-int(coord[0])
        intron_num += 1
print(intronsum/intron_num)
