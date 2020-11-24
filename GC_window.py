#This project aimed at calculating and plotting GC content in a DNA sequence in a given user-defined window and step sizes.

---------------------------
#!/usr/bin/env python3

import sys
import matplotlib.pyplot as plt

myfragment = []
mylist = []

filename = input('Please enter your fasta file: ')
fasta = open(filename,'r')
myfragment = []
mylist = []
for line in fasta:
    if not line.startswith('>'):
        seq = line.rstrip()
        myfragment.append(seq)
sequence = ''.join(myfragment)#after all lines are appended they got joined --> outside the for loop.

windowSize = int(input('Enter your window size: '))
if windowSize > len(sequence):
    print('**ERROR** winSize must not be larger than sequence length.', file = sys.stderr)
    print(int(input('Try again! Enter your window size which is less than {}'.format(filename))))

step = int(input('Enter your step size: '))
if step > windowSize:
    print('**ERROR** step must not be larger than winSize.', file = sys.stderr)
    print(int(input('Enter your step size again: ')))
    sys.exit()

for i in range(0,len(sequence)- windowSize+1,step):# i go out of the first for loop because i need it to continue joining all the sequences in the file then i work on the window.
#location of the start of window is i. (at first iteration starts from 0)
    window = sequence[i:i+windowSize]
    GC_count = round((window.count('G') + window.count('C')) /len(window) * 100,4)
    mylist.append(GC_count)

# the last number in the range is excluded.

plt.figure(figsize=(15,5)) # im opening a new figure over the previous one. thats why when i keep plt.subplots() then im making a new figure over that one. so i must remove fig, ax = ..etc.
plt.plot([i for i in range(0,len(sequence)-windowSize+1,step)],mylist) #i represents the start of the window. # this syntax creates a list of 'i's which are the start of each window.
# ax.set_ylim([0,100])
# ax.set_xlim([0,50000])

plt.title('GC content over BRCA1 gene sequence')
plt.xlabel('bps')
plt.ylabel('GC(%)')
plt.savefig("gc_plot.png")

fasta.close()
