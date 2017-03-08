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

for line in sumfile:
    toks = line.strip().split()
    tp = (toks[0], toks[1])
    if tp in tagpair:
        usedpair.append(line)
    elif toks[0] in usedtags and toks[1] in usedtags:
        bothused.append(line)
    elif toks[0] in usedtags or toks[1] in usedtags:
        oneused.append(line)
    else:
        unused.append(line)

sumfile.close()

outfile = open(args.output, "w")
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

        
