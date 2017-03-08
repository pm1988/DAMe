#! /usr/bin/env python
###############################################################
# Python script to plot the per sample length distribution and#
# make a summary file with the counts from each pool.         #
###############################################################

import numpy as np
import argparse

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Make a sample summary file and plot length metrics.")
    parser.add_argument("-f", "--fasta", type=str, dest="fasta", metavar="FastaFile", required=True, help="Filtered reads fasta file.")
    parser.add_argument("-p", "--psinfo", type=str, dest="psinfo", metavar="PSInfoFile", required=True, help="PS Information file.")
    parser.add_argument("-n", "--numpools", type=int, dest="pools", metavar="NumPools", required=True, help="Number of pools.")
    args = parser.parse_args()


## Read and process the ps info file to get tag-pool combo
psinfo = open(args.psinfo)
sample_pool = {}
for line in psinfo:
    toks = line.strip().split()
    if toks[0] in sample_pool:
        sample_pool[toks[0]][toks[1]+"-"+toks[2]] = int(toks[3])-1
    else:
        sample_pool[toks[0]] = {}
        sample_pool[toks[0]][toks[1]+"-"+toks[2]] = int(toks[3])-1
psinfo.close()
    
## Read and process the fasta file.
fasta = open(args.fasta)
pool_lens = [ {} for x in xrange(args.pools)]
sample_pool_counts = {}
total_count = 0

for line in fasta:
    if line[0] == '>':
        (sample, tagcombs, counts) = line[1:].strip().split()
        if sample not in sample_pool_counts:
            tagnames = tagcombs.split("_")[0].split(".")
            temp = [int(x) for x in counts.split("_")]
            sample_pool_counts[sample] = [0]*args.pools
            for tagname, cnt in zip(tagnames, temp):
                if tagname in sample_pool[sample]:
                    sample_pool_counts[sample][sample_pool[sample][tagname]] = cnt
            total_count = np.sum(temp)
        else:
            tagnames = tagcombs.split("_")[0].split(".")
            temp = [int(x) for x in counts.split("_")]
            for tagname, cnt in zip(tagnames, temp):
                if tagname in sample_pool[sample]:
                    sample_pool_counts[sample][sample_pool[sample][tagname]] += cnt
            total_count = np.sum(temp)
    else:
        seqlen = len(line.strip())
        for tagname, cnt in zip(tagnames, temp):
            if tagname in sample_pool[sample]:
                curpool = sample_pool[sample][tagname]
                if seqlen not in pool_lens[curpool]:
                    pool_lens[curpool][seqlen] = cnt
                else:
                    pool_lens[curpool][seqlen] += cnt

fasta.close()

######## Print out per sample summary
outfile = open("SampleVsPools.summary.txt", "w")
outfile.write("Sample")
for i in xrange(args.pools):
    outfile.write("\tPool"+str(i+1))
    
outfile.write("\n")

for sample in sample_pool_counts:
    outfile.write(sample)
    for pool in xrange(args.pools):
        if pool in sample_pool[sample].values():
            outfile.write("\t"+str(sample_pool_counts[sample][pool]))
        else:
            outfile.write("\t--")
    outfile.write("\n")

outfile.close()
######## Do the plotting
import matplotlib
matplotlib.use('Agg')
from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.pyplot as plt

pdf = PdfPages('SequenceLengthDistribution.pdf')
colors = ['firebrick', 'slateblue', 'seagreen', 'peru', 'teal', 'darkorchid', 'goldenrod', 'yellowgreen']

pool_lens_combined = {}
for pool in xrange(args.pools):
    for seqlen in pool_lens[pool]:
        if seqlen in pool_lens_combined:
            pool_lens_combined[seqlen] += pool_lens[pool][seqlen]
        else:
            pool_lens_combined[seqlen] = pool_lens[pool][seqlen]
#    print(pool_lens_combined)

plt.figure()
plt.bar(pool_lens_combined.keys(), pool_lens_combined.values(), color='slategray', width=1, log=True)
plt.title('Sequence length distribution: Pools combined')
plt.ylabel('Number of sequences')
plt.xlabel('Length of sequences')
plt.yscale('log')
pdf.savefig()  # saves the current figure into a pdf page
plt.close()

for pool in xrange(args.pools):
    plt.figure()
    if not pool_lens[pool]:
        print "Warning: Pool %d is empty. Make sure that this is correct."%(pool)
        continue
    plt.bar(pool_lens[pool].keys(), pool_lens[pool].values(), color=colors[pool%8], width=1, log=True)
    plt.title('Sequence length distribution: Pool '+str(pool+1))
    plt.ylabel('Number of sequences')
    plt.xlabel('Length of sequences')
    plt.yscale('log')
    pdf.savefig()  # saves the current figure into a pdf page
    plt.close()
    
pdf.close()
