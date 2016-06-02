#! /usr/bin/env python

#######################################################
# File: assessClusteringParameters.py                 #
# Author: Shyam Gopalakrishnan                        #
# Date: 12th February 2016                            #
# Description: This script runs sumatra and sumaclust #
# for multiple values of identity cutoff and abundance#
# proportions to generate plots that can used to      #
# choose values for these parameters.                 #
#######################################################

import argparse
import subprocess
import numpy as np
import sys

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Generate plots to explore the parameter space for OTU clustering parameters.")
    parser.add_argument("-i", "--inFasta", dest="fasta", type=str, metavar="InputFasta", required=True, help='Input fasta file for cluster.')
    parser.add_argument("-mint", "--minIdentityCutoff", dest="mint", type=float, metavar="minT", required=False, help='Minimum identity cutoff for sumatra/sumaclust', default=0.7)
    parser.add_argument("-minR", "--minAbundanceRatio", dest="minR", type=float, metavar="minR", required=False, help='Minimum abundance ratio for sumatra/sumaclust', default=0.5)
    parser.add_argument("-step", "--stepSize", dest="step", type=float, metavar="stepSize", required=False, help='Step size for t and R', default=0.01)
    parser.add_argument("-t", "--threads", dest="threads", type=int, metavar="threads", required=False, help='Number of threads to use', default=4)
    parser.add_argument("-o", "--out", dest="out", type=str, metavar="outfile", required=False, help="Output file pdf", default="")
    args = parser.parse_args()
    if args.out == "":
        args.out = "ClusteringParms.t_"+str(args.mint)+".R_"+str(args.minR)+".step_"+str(args.step)+".pdf"

## First run sumatra to generate the file with the pairwise information.
lines = subprocess.check_output(("sumatra -p "+ str(args.threads) +" -a -t "+str(args.mint)+" "+args.fasta).split())
lines = lines.split('\n')
idents = [float(x.strip().split()[-1]) for x in lines if x.strip() != ""]

ts = [ round(x, 2) for x in np.arange(args.mint, 1.0001, args.step) ]
Rs = [ round(x, 2) for x in np.arange(args.minR, 1.0001, args.step) ]
##Rs = Rs[::-1]
outs = []
print "Running Sumaclust ..."
print "0% done.",
totnum = 0.0
increment = 100.0/(len(ts)*len(Rs))
for t in ts:
    for R in Rs:
        if t > 0.999:
            lines = subprocess.check_output(("sumaclust -a -t 0.9999 -R "+str(R)+" -e -f "+args.fasta).split(), stderr=subprocess.STDOUT)
        else:
            lines = subprocess.check_output(("sumaclust -a -t "+str(t)+" -R "+str(R)+" -e -f "+args.fasta).split(), stderr=subprocess.STDOUT)
        lines = lines.split("\n")
#        print t, R, "->", lines[-2].strip().split()[-3]
        outs.append(int(lines[-2].strip().split()[-3]))
        totnum += increment
        print "\r%.2lf%% done."%totnum,
        sys.stdout.flush()
print "\n",
outs = np.array(outs)
outmean = int(np.mean(outs))
outmedian = int(np.median(outs))
outs[outs > outmean] = outmean
outs.resize((len(ts), len(Rs)))
#print outs

import matplotlib
matplotlib.use('Agg')

from matplotlib.backends.backend_pdf import PdfPages
from matplotlib import cm
import matplotlib.pyplot as plt

pdffile = PdfPages(args.out)

plt.hist(idents, 20, facecolor='green', alpha=0.75)
plt.grid(True)
plt.xlabel('Similarity measure')
plt.ylabel('Number of pairs')
plt.title('Distribution of pairwise similarity measures')
pdffile.savefig()

plt.figure()
cax = plt.pcolormesh(np.log10(outs), cmap=cm.rainbow)
cb = plt.colorbar(cax, ticks=np.log10([np.min(outs), outmean, outmedian]))
cb.ax.set_yticklabels([np.min(outs), "> "+str(outmean), outmedian])
plt.ylabel('Similarity (t)')
plt.xlabel('Abundance proportion (R)')
plt.xticks(np.arange(len(Rs))+0.5, Rs, rotation='vertical')
plt.yticks(np.arange(len(ts))+0.5, ts)
plt.xlim((0, len(Rs)))
plt.ylim((0, len(ts)))
plt.title('Number of OTUs (log-scale) vs clustering parameters')
pdffile.savefig()

counts = []
for R in Rs:
    lines = subprocess.check_output(("sumaclust -a -t 0.97 -R "+str(R)+" -e -f "+args.fasta).split(), stderr=subprocess.STDOUT)
    lines = lines.split("\n")
    counts.append(int(lines[-2].strip().split()[-3]))
plt.figure()
plt.plot(Rs, counts, "r-", linewidth=2)
plt.plot(Rs, counts, "ro")
plt.axis([np.min(Rs)*0.9, np.max(Rs)*1.1, 0, np.max(counts)*1.1])
plt.xlabel('Abundance proportion')
plt.ylabel('Number of OTUs found')
plt.title('Number of OTUs vs R for t=0.97')
pdffile.savefig()

counts = []
for t in ts:
    if t > 0.999:
        lines = subprocess.check_output(("sumaclust -a -t 0.9999 -R 0.8 -e -f "+args.fasta).split(), stderr=subprocess.STDOUT)
    else:
        lines = subprocess.check_output(("sumaclust -a -t "+str(t)+" -R 0.8 -e -f "+args.fasta).split(), stderr=subprocess.STDOUT)
    lines = lines.split("\n")
    counts.append(int(lines[-2].strip().split()[-3]))
plt.figure()
plt.plot(ts, counts, "b-", linewidth=2)
plt.plot(ts, counts, "bo")
plt.axis([np.min(ts)*0.9, np.max(ts)*1.1, 0, np.max(counts)*1.1])
plt.xlabel('Similarity cutoff')
plt.ylabel('Number of OTUs found')
plt.title('Number of OTUs vs t for R=0.8')
pdffile.savefig()

pdffile.close()
