#! /usr/bin/env python
#####################################################
# Helper script to parse the SummaryCounts.txt file #
# by the used, 1-used and unused tag combinations.  #
#####################################################

import sys
import argparse

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Pretty sort summary counts using info file')
    parser.add_argument("-p", "--psinfo", dest="psinfo", type=str, metavar="PSInfoFile", required=True, help='PS info file')
    parser.add_argument("-l", "--pool", dest="pool", type=str, metavar="Poolname", required=True, help='Name of the pool')
    parser.add_argument("-s", "--summary", dest="summary",type=str, metavar="SummaryCountsFile", required=True, help='Summary counts file')
    parser.add_argument("-o", "--outfile", dest="output", type=str, metavar="OutputFile", required=True, help='Output file')
    args = parser.parse_args()

### Parse the psinfo file to get all tags in there for this pool
psfile = open(args.psinfo)
usedtags = []
tagpair = []

for line in psfile:
    toks = line.strip().split()
    if len(toks) != 4:
        print "Incorrect number of columns in ps info file in line:'", line
        sys.exit(1)
    if toks[3] != args.pool:
        ## not the correct pool, so just skip the line
        continue
    tp = (toks[1], toks[2])
    if tp not in tagpair:
        tagpair.append(tp)
    if toks[1] not in usedtags:
        usedtags.append(toks[1])
    if toks[2] not in usedtags:
        usedtags.append(toks[2])

psfile.close()

## Parse the summary counts file and output the stuff up, split into combinations
## where both tags are used, only one tag is used or neither is used.
sumfile = open(args.summary)
sumfile.readline() ## discard header
usedpair = []
bothused = []
oneused = []
unused = []

countsTots = [0,0,0,0]
countsUniq = [0,0,0,0]

for line in sumfile:
    toks = line.strip().split()
    tp = (toks[0], toks[1])
    if tp in tagpair:        
        usedpair.append(line)
        countsTots[0] += int(toks[3])
        countsUniq[0] += int(toks[2])
    elif toks[0] in usedtags and toks[1] in usedtags:
        bothused.append(line)
        countsTots[1] += int(toks[3])
        countsUniq[1] += int(toks[2])
    elif toks[0] in usedtags or toks[1] in usedtags:
        oneused.append(line)
        countsTots[2] += int(toks[3])
        countsUniq[2] += int(toks[2])
    else:
        unused.append(line)
        countsTots[3] += int(toks[3])
        countsUniq[3] += int(toks[2])

sumfile.close()

pctTots = [str(round(x*100.0/sum(countsTots),2)) for x in countsTots]
pctUniq = [str(round(x*100.0/sum(countsUniq),2)) for x in countsUniq]

outfile = open(args.output, "w")
outfile.write("                                                    \tTotal seqs\tTotal unique seqs\t% Total seqs\t%Total unique seqs\n")
outfile.write("Tag combinations where the tag pair was used        \t"+str(countsTots[0])+"\t"+str(countsUniq[0])+"\t"+pctTots[0]+"\t"+pctUniq[0]+"\n")
outfile.write("Tag combinations where both tags used\n")
outfile.write("but not in this combination                         \t"+str(countsTots[1])+"\t"+str(countsUniq[1])+"\t"+pctTots[1]+"\t"+pctUniq[1]+"\n")
outfile.write("Tag combinations where only one of the tags was used\t"+str(countsTots[2])+"\t"+str(countsUniq[2])+"\t"+pctTots[2]+"\t"+pctUniq[2]+"\n")
outfile.write("Tag combinations where neither tag was used         \t"+str(countsTots[3])+"\t"+str(countsUniq[3])+"\t"+pctTots[3]+"\t"+pctUniq[3]+"\n")
outfile.write("Total                                               \t"+str(sum(countsTots))+"\t"+str(sum(countsUniq))+"\t100.00\t100.00\n\n")

outfile.write("Tag combinations where the tag pair was used.\n")
outfile.write("---------------------------------------------\n")
outfile.write("Tag1Name\tTag2Name\tNumUniqSeqs\tSumTotalFreq\n")
outfile.write("".join(usedpair))
outfile.write("\n")
outfile.write("Tag combinations where both tags used - but not in this combination.\n")
outfile.write("--------------------------------------------------------------------\n")
outfile.write("Tag1Name\tTag2Name\tNumUniqSeqs\tSumTotalFreq\n")
outfile.write("".join(bothused))
outfile.write("\n")
outfile.write("Tag combinations where only one of the tags was used.\n")
outfile.write("-----------------------------------------------------\n")
outfile.write("Tag1Name\tTag2Name\tNumUniqSeqs\tSumTotalFreq\n")
outfile.write("".join(oneused))
outfile.write("\n")
outfile.write("Tag combinations where neither tag was used.\n")
outfile.write("--------------------------------------------\n")
outfile.write("Tag1Name\tTag2Name\tNumUniqSeqs\tSumTotalFreq\n")
outfile.write("".join(unused))



outfile.close()

        
