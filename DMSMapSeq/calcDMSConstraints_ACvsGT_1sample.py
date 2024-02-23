# ----------------------------------------------------------------
# Input: pkl files, sample name in command
# Output: TSV file of DMS reactivities, in SHAPE reactivity format
# Reactivity calculation:
#  1) Mutation rate at each position = # of mutations/# of reads
#     AC and GT mutation rates are calculated and stored separately.
#
#
# ----------------------------------------------------------------
#!/usr/bin/env python3

import numpy as np
import pickle
import sys
from Bio import SeqIO
import os
import time
import csv

start = time.time()

sample = sys.argv[1] #Denatured, DMS1, DMS2, DMS3, or noDMS
path = os.getcwd()
#print(path)

pool_names = []
seqs = []

with open(path+"/short_long_MS2.fa","r") as fasta_h:
    for i in SeqIO.parse(fasta_h,"fasta"):
        pool_names.append(i.id)
        seqs.append(str(i.seq))

#print(len(seqs))
#print(seqs[1])

timepoint0 = time.time()
with open(path+"/"+sample+"__barcode_Bar2_1_mut_cov_PE.pkl","rb") as f:
    _, _, mut, cov, _ = pickle.load(f)

print(len(mut[0])) # should be 300

timepoint1 = time.time()

low_cov_seq_count = 0
for i in range(len(seqs)): # Loop through all sequences in the pool

    threshold = 400
    seq = seqs[i][19:290]
    pass_threshold = [True]*len(seq)
    #print(seq)

    cover = cov[i][19:290]

    for j in np.arange(len(seq)): # Loop through all nt in a sequence
        if cover[j] < threshold:
            pass_threshold[j] = False

    if np.average(cover) < threshold:
        low_cov_seq_count += 1
        continue

    # Calculate raw reactivities for all replicates, setting low coverage and G/T positions as nan
        
    mutat = mut[i][19:290]
        
    reactsAC = []
    reactsGT = []
    for j in np.arange(len(seq)):
        if seq[j] in ["A","a","C","c"]: # AC nucleotides
            reactsGT.append(np.nan)
            if not pass_threshold[j]:
                reactsAC.append(np.nan)
            else:
                raw_react = mutat[j]/cover[j]
                reactsAC.append(raw_react)
        else: # GT nucleotides
            reactsAC.append(np.nan)
            if not pass_threshold[j]:
                reactsGT.append(np.nan)
            else:
                raw_react = mutat[j]/cover[j]
                reactsGT.append(raw_react)

    name = pool_names[i]

    # Write positions and reactivities to a csv file (SHAPE reactivity format)
    namecol = [name]*len(seq)
    pos = np.arange(1, len(seq)+1)

    c = open(path+f"/raw_reacts/{sample}_rawreactsAC.tsv","a")
    writer = csv.writer(c, delimiter="\t")
    writer.writerows(zip(namecol,pos,reactsAC))
    c.close()
    
    d = open(path+f"/raw_reacts/{sample}_rawreactsGT.tsv","a")
    writer = csv.writer(d, delimiter="\t")
    writer.writerows(zip(namecol,pos,reactsGT))
    d.close()

timepoint2 = time.time()

print("Time taken to parse fasta:"+str(timepoint0-start))
print("Time taken to unpickle:"+str(timepoint1-timepoint0))
print("Time taken to loop all sequences:" +str(timepoint2-timepoint1))
print("low coverage sequence count:" + str(low_cov_seq_count))
        

