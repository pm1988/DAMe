##Import libraries

import re
import os
import sys
from string import maketrans

##Define the functions

def RC(seq):
	seq=seq[::-1] # reverse seq
	transtab=maketrans('ACGTMRWSYKVHDB', 'TGCAKYWSRMBDHV')
	return seq.translate(transtab)


def readTags(tags, TAGS):
	file=open(tags)
	line= file.readline()
        tagLength = -1
	while line:
		line=line.strip().split()
		#Initialize it in the dictionary
		if not TAGS.has_key(line[1]): 
			TAGS[line[1]]=[] 
		#Fill it
		TAGS[line[1]].append(line[0])
		TAGS[line[1]].append(RC(line[0]))
                if tagLength == -1:
                        tagLength = len(line[0])
                elif tagLength == 0:
                        pass
                elif tagLength != len(line[0]):
                        tagLength = 0
		#Read next line
		line= file.readline()   
	file.close()
	return (TAGS, tagLength)


def readPrimers(primers, PRIMERS, AMBIG):
	file=open(primers)
	line=file.readline()   
	while line:
		line=line.rstrip().split()
		#Initialize it in the dictionary
		if not PRIMERS.has_key(line[0]): 
			PRIMERS[line[0]]=[[],[]] 
		# Get the regExp equivalents
		Frc=RC(line[1])
		Rrc=RC(line[2])
		F=line[1]
		R=line[2]
		for key in AMBIG:
			Frc=re.sub(key,AMBIG[key], Frc)
			Rrc=re.sub(key,AMBIG[key], Rrc)
			F=re.sub(key,AMBIG[key], F)
			R=re.sub(key,AMBIG[key], R)
		#Fill it
		PRIMERS[line[0]][0].append(F)   ## Possible on A side. Forward normal (regExp)
		PRIMERS[line[0]][0].append(R)   ## Possible on A side. Reverse normal (regExp)
		PRIMERS[line[0]][1].append(Frc)   ## Possible on A side. RC(Forward) (regExp)
		PRIMERS[line[0]][1].append(Rrc)   ## Possible on A side. RC(Reverse) (regExp)
		#Read next line
		line= file.readline()   
	file.close()
	return PRIMERS


def findHammingDistance(seq1, seq2):
        """This function finds the difference between two given
        sequences. If the sequences are different lengths, then
        it returns a large difference.
        """
        if len(seq1) != len(seq2):
                return(50)
        dist = 0
        for index in range(len(seq1)):
                dist += (seq1[index] != seq2[index])
        return dist

def findBestTag(seq, TAGS, reverse=False):
        """This function finds the best tag for the sequence
        irrespective of the length. If there are multiple tags
        with the same number of mismatches, then it returns the
        longest one as the best match.
        """
        bestDist = 100
        bestTag = ""
        bestLength = 0
        numBest = 0
        if reverse:
                index = 1
        else:
                index = 0
        for tag in TAGS:
                curtaglen = len(TAGS[tag][0])
                dist = findHammingDistance(seq, TAGS[tag][index])
                if dist < bestDist or (dist == bestDist and curtaglen > bestLength):
                        bestDist = dist
                        bestTag = tag
                        bestLength = curtaglen
                        numBest = 1
                elif dist == bestDist and curtaglen == bestLength:
                        numBest += 1
        if bestDist != 0 or numBest != 1:
                return []
        else:
                return [bestTag]


def GetPiecesInfo(line, PRIMERS, TAGS, keepPrimersSeq, tagLength, keepLongest=False):
        ## Key change made on 2 June 2016 Shyam
        ## the length of the tags must be the same for
        ## all the tags.
        ## IF NOT, then the program will __not work__
        ## Change 2 on 10th July 2016 - tags of different
        ## lengths will work if you allow the best tag
        ## to be chosen based on length of match.
        ## Another key change - 10th July 2016 - now the tag
        ## primer combination does not have to be at the start
        ## of the read. You are allowed to have some space
        ## between start of read and start of tag+primer.
	for key in PRIMERS:
		#For a forward seq
		primIniPos= [(m.start(0), m.end(0)) for m in re.finditer(PRIMERS[key][0][0], line)]
                ### primIniPos contains the start and end positions of forward primer match in the sequence ++ FR
		if len(primIniPos)>0: # If it's a forward read and there are no seq errs in the primer seq
			if keepPrimersSeq: 		
				primIniPosPrim=primIniPos[0][0]   
				primIniPosTags=primIniPos[0][0]
			else:
				primIniPosPrim=primIniPos[0][1] 
				primIniPosTags=primIniPos[0][0]
			primFinPos=[(m.start(0), m.end(0)) for m in re.finditer(PRIMERS[key][1][1], line)]
			if len(primFinPos)>0: #If there's no seq error in the primer seq
				if keepPrimersSeq: 		
					primFinPosPrim=primFinPos[0][1]   
					primFinPosTags=primFinPos[0][1]
				else:
					primFinPosPrim=primFinPos[0][0] 
					primFinPosTags=primFinPos[0][1]
				PrimerName=key
				between=line[primIniPosPrim:primFinPosPrim]  
				if len(between)==0:    #If the seq only contained primers and no amplified barcode
					return [1]
				else:
					#Figure out the names of the tags
                                        if tagLength == 0:
                                                tagName1 = findBestTag(line[:primIniPosTags], TAGS)
                                                tagName2 = findBestTag(line[primFinPosTags:], TAGS, True)
                                        else:
                                                tagName1= [tagName for tagName in TAGS if TAGS[tagName][0]==line[(primIniPosTags-tagLength):primIniPosTags]]
                                                tagName2= [tagName for tagName in TAGS if TAGS[tagName][1]==line[primFinPosTags:(primFinPosTags+tagLength)]]
					if len(tagName1)>0 and len(tagName2)>0:
						tagName1=tagName1[0]
						tagName2=tagName2[0]
					else: #If there was a seq error in the tags
						return [1]
			else:
				return [1]
		else:  #For a reverse seq
			primIniPos= [(m.start(0), m.end(0)) for m in re.finditer(PRIMERS[key][0][1], line)]
                        ### Now the primIniPos contains the start and end positions of the RC primer match in the seq ++ RC
			if len(primIniPos)>0:
				if keepPrimersSeq: 		
					primIniPosPrim=primIniPos[0][0]   
					primIniPosTags=primIniPos[0][0]
				else:
					primIniPosPrim=primIniPos[0][1] 
					primIniPosTags=primIniPos[0][0]
				primFinPos=[(m.start(0), m.end(0)) for m in re.finditer(PRIMERS[key][1][0], line)]
				if len(primFinPos)>0: #If there's no seq error in the primer seq
					if keepPrimersSeq: 		
						primFinPosPrim=primFinPos[0][1]   
						primFinPosTags=primFinPos[0][1] 
					else:
						primFinPosPrim=primFinPos[0][0] 
						primFinPosTags=primFinPos[0][1] 
					PrimerName=key
					between=line[primIniPosPrim:primFinPosPrim]  
					if len(between)==0:
						return [1]
					else:
						between=RC(between)  ## Since this is a reverse seq, reverse complement the sequence
						#Figure out the names of the tags
						# here the tags used are reversed because the seq was reversed. Tag1 is tag2 and tag2 is tag1
					        #Figure out the names of the tags
                                                if tagLength == 0:
                                                        tagName1 = findBestTag(line[:primIniPosTags], TAGS)
                                                        tagName2 = findBestTag(line[primFinPosTags:], TAGS, True)
                                                else:
                                                        tagName2= [tagName for tagName in TAGS if TAGS[tagName][0]==line[(primIniPosTags-tagLength):primIniPosTags]]  
                                                        tagName1= [tagName for tagName in TAGS if TAGS[tagName][1]==line[primFinPosTags:(primFinPosTags+tagLength)]]  
						if len(tagName1)>0 and len(tagName2)>0:
							tagName1=tagName1[0]
							tagName2=tagName2[0]
						else: #If there was a seq error in the tags
							return [1]
				else:
					return [1]
			else:
				return [1]
	Info=[tagName1, tagName2 , PrimerName, between]
	return Info



def FillHAP(HAP, tagName1, tagName2 , PrimerName, between):
	#Add tags to HAP if it doesn't have it already.  
	tagHapKey="_".join([tagName1, tagName2])
	if not HAP.has_key(tagHapKey):   ##key is the Tag1_Tag2 
		HAP[tagHapKey]=[[],[],{}] #tag1, tag2, Seq[Freq, PrimerName]  
		HAP[tagHapKey][0]=tagName1
		HAP[tagHapKey][1]=tagName2
	#Fill it
	#Set the seq and its counter if it doesn't have it, alse add it to its seq key counter
	if not HAP[tagHapKey][2].has_key(between): 
		HAP[tagHapKey][2][between]=[1]
		HAP[tagHapKey][2][between].append(PrimerName)
	else:
		HAP[tagHapKey][2][between][0]+=1
	return HAP



def PrintSortedCollapsedCountedSeqs(HAP):
	for TagComb in HAP:
		out=open("%s.txt"%(TagComb), "w")
		tagName1=HAP[TagComb][0]
		tagName2=HAP[TagComb][1]
		for Seq in HAP[TagComb][2]:
			a="\t".join([HAP[TagComb][2][Seq][1], tagName1, tagName2, str(HAP[TagComb][2][Seq][0]), Seq])
			out.write("%s\n"%(a))
		out.close()



def PrintSummaryFile(HAP):
	out=open("SummaryCounts.txt", "w")
	a="\t".join( ("#tagName1", "tagName2", "NumUniqSeqs", "SumTotalFreq" ))
	out.write("%s\n"%(a))
	for TagComb in HAP:
		tagName1=HAP[TagComb][0]
		tagName2=HAP[TagComb][1]
		NumUniqSeqs=len(HAP[TagComb][2])
		SumTotalFreq=0
		for Seq in HAP[TagComb][2]:
			SumTotalFreq= SumTotalFreq + HAP[TagComb][2][Seq][0]
		a="\t".join( (tagName1, tagName2, str(NumUniqSeqs), str(SumTotalFreq)))
		out.write("%s\n"%(a))
	out.close()




