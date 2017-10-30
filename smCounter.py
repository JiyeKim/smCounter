# smCounter: barcode-aware variant caller
# Chang Xu. 23May2016; online version of 10APR2017
import os
import datetime
import subprocess
import math
import operator
import argparse
import random
import multiprocessing
import traceback
from collections import defaultdict

# 3rd party modules
import pysam
import scipy.stats

# global contants (note that multiprocessing imports this .py file, so do not do much work outside functions)
pcr_error = 1e-6
pcr_no_error = 1.0 - 3e-5

atgc = ('A', 'T', 'G', 'C')

#-------------------------------------------------------------------------------------
# function to calculate posterior probability for each barcode. 
#-------------------------------------------------------------------------------------
def calProb(oneBC, mtDrop):
    # oneBC: read information of one barcode
    # oneBC: {'3LMDR:00198:00085': [['G', 0.00039810717055349735, None]]}

    print 'calProb()-----'
    print 'oneBC:{}, mtDrop:{}'.format(oneBC, mtDrop)

    outDict = defaultdict(float)
    if len(oneBC) <= mtDrop:  # number of reads <= mtDrop
        outDict['A'] = 0.0
        outDict['T'] = 0.0
        outDict['G'] = 0.0
        outDict['C'] = 0.0
    else:
        prodP = defaultdict(float)
        cnt = defaultdict(int)
        tmpOut = defaultdict(float)
        rightP = 1.0
        sumP = 0.0
        pcrP = defaultdict(float)

        # set ATGC count = 0
        for char in atgc:
          cnt[char] = 0

        # get unique bases. Make sure uniqBaseList contains 4 members, unless the barcode already contains more than or equal to 4 bases/indels
        # NOTE: existBase contains only the alleles, including indels, with at least 1 read in the MT. uniqBase may contain more. 
        existBase = set([info[0][0] for info in oneBC.values()])  # 'G'
        uniqBase = set([info[0][0] for info in oneBC.values()]) # 'G'
        if len(uniqBase) < 4:
          for b in atgc:    # ('A', 'T', 'G', 'C')
              if b not in uniqBase:
                  uniqBase.add(b)
                  if len(uniqBase) == 4:
                      # When the uniqBase has 'A', 'T', 'G', 'C', it will be stop.
                      break

        uniqBaseList = list(uniqBase)
        print 'uniqBase:{}'.format(uniqBase)
        print 'uniqBaseList:{}'.format(uniqBaseList)

        # set initial value in prodP to be 1.0
        for b in uniqBaseList:
            prodP[b] = 1.0

        # prodP = {'A': 1.0, 'C': 1.0, 'T': 1.0, 'G': 1.0}

        for info in oneBC.values():
            # info = [['G', 0.00039810717055349735, None]]
            base = info[0][0]   # 'G'
            # prob is the error probability
            prob = info[0][1]  # 0.00039810717055349735
            pairOrder = info[0][2]  # pairOrder; None

            if pairOrder not in [None, 'Paired']:
                prob = 0.1

            # prodP is the probability of no sequencing error for each base
            prodP[base] *= 1.0 - prob
            cnt[base] += 1

            for char in list(uniqBase - set([base])):
                prodP[char] *= prob

            # rightP is the probabilty that there is no sequencing error, hence the alternative alleles come from PCR error
            rightP *= 1.0 - prob

        '''
        print 'prodP:', prodP
        print 'rightP:', rightP
        print 'cnt -->', cnt
        '''

        for char in uniqBaseList:   # ['A', 'T', 'G', 'C']
            ratio = (cnt[char] + 0.5) / (len(oneBC) + 0.5 * len(uniqBaseList))
            pcrP[char] = 10.0 ** (-6.0 * ratio)

        for key in prodP.keys():
            if key in existBase:
                # tmpOut[key] is P(BC|key), or the likelihood of all reads in the barcode, given the true allele is *key*.  
                tmpOut[key] = pcr_no_error * prodP[key] + rightP * min([pcrP[char] for char in pcrP.keys() if char != key])
            else:
                tmpOut[key] = rightP
                for char in existBase:
                    if char != key:
                        tmpOut[key] *= pcrP[char] 

            sumP += tmpOut[key]

        for key in prodP.iterkeys():
            outDict[key] = 0.0 if sumP <= 0 else tmpOut[key] / sumP
    return outDict

#-------------------------------------------------------------------------------------
# convert variant type, reference base, variant base to output format
#-------------------------------------------------------------------------------------
def convertToVcf(origRef,origAlt):
    vtype = '.'
    ref = origRef
    alt = origAlt
    if len(origAlt) == 1:
        vtype = 'SNP'
    elif origAlt == 'DEL':
        vtype = 'SDEL'
    else:
        vals = origAlt.split('|')
        if vals[0] in ('DEL', 'INS'):
            vtype = 'INDEL'
            ref = vals[1]
            alt = vals[2]
    return (ref, alt, vtype)

#-------------------------------------------------------------------------------------
# check if a locus is within or flanked by homopolymer region and/or low complexity region
#-------------------------------------------------------------------------------------
def isHPorLowComp(chrom, pos, length, refb, altb, refGenome):
    # get reference bases for interval [pos-length, pos+length]
    refs = pysam.FastaFile(refGenome)
    chromLength = refs.get_reference_length(chrom)
    pos0 = int(pos) - 1
    Lseq        = refs.fetch(reference=chrom, start=max(0,pos0-length)  , end=pos0).upper()
    Rseq_ref = refs.fetch(reference=chrom, start=       pos0+len(refb) , end=min(pos0+len(refb)+length,chromLength)).upper()
    Rseq_alt = refs.fetch(reference=chrom, start=       pos0+len(altb) , end=min(pos0+len(altb)+length,chromLength)).upper()
    refSeq = Lseq + refb + Rseq_ref
    altSeq = Lseq + altb + Rseq_alt
    # check homopolymer
    homoA = refSeq.find('A'*length) >= 0 or altSeq.find('A'*length) >= 0
    homoT = refSeq.find('T'*length) >= 0 or altSeq.find('T'*length) >= 0
    homoG = refSeq.find('G'*length) >= 0 or altSeq.find('G'*length) >= 0
    homoC = refSeq.find('C'*length) >= 0 or altSeq.find('C'*length) >= 0 
    homop = homoA or homoT or homoG or homoC

    # check low complexity -- window length is 2 * homopolymer region. If any 2 nucleotide >= 99% 
    len2 = 2 * length
    LseqLC    = refs.fetch(reference=chrom, start=max(0,pos0-len2)   , end=pos0).upper()
    Rseq_refLC = refs.fetch(reference=chrom, start=       pos0+len(refb), end=min(pos0+len(refb)+len2,chromLength)).upper() # ref seq
    Rseq_altLC = refs.fetch(reference=chrom, start=       pos0+len(altb), end=min(pos0+len(altb)+len2,chromLength)).upper() # alt seq
    refSeqLC = LseqLC + refb + Rseq_refLC
    altSeqLC = LseqLC + altb + Rseq_altLC
    lowcomp = False
    
    # Ref seq
    totalLen = len(refSeqLC)
    for i in range(totalLen-len2):
        subseq = refSeqLC[i:(i+len2)]
        countA = subseq.count('A')
        countT = subseq.count('T')
        countG = subseq.count('G')
        countC = subseq.count('C')
        sortedCounts = sorted([countA, countT, countG, countC], reverse=True)
        top2Freq = 1.0 * (sortedCounts[0] + sortedCounts[1]) / len2
        if top2Freq >= 0.99:
            lowcomp = True
            break
        
    # If ref seq is not LC, check alt seq
    if not lowcomp:
        totalLen = len(altSeqLC)
        for i in range(totalLen-len2):
            subseq = altSeqLC[i:(i+len2)]
            countA = subseq.count('A')
            countT = subseq.count('T')
            countG = subseq.count('G')
            countC = subseq.count('C')
            sortedCounts = sorted([countA, countT, countG, countC], reverse=True)
            top2Freq = 1.0 * (sortedCounts[0] + sortedCounts[1]) / len2
            if top2Freq >= 0.99:
                lowcomp = True
                break

    return (homop, lowcomp)

#-------------------------------------------------------------------------------------
# filter variants
#-------------------------------------------------------------------------------------
def filterVariants(ref,alt,vtype,origAlt,origRef,usedMT,strongMTCnt,chrom,pos,hpLen,refGenome,MTCnt,alleleCnt,cvg,discordPairCnt,concordPairCnt,reverseCnt,forwardCnt,lowQReads,r1BcEndPos,r2BcEndPos,r2PrimerEndPos,primerDist):
    # init output string
    fltr = ';'
    
    # low coverage filter
    if usedMT < 5:
        fltr += 'LM;' 

    # low number of strong MTs filter
    if strongMTCnt[origAlt] < 2 :
        fltr += 'LSM;'

    # check region for homopolymer or low complexity
    (isHomopolymer,isLowComplexity) = isHPorLowComp(chrom, pos, hpLen, ref, alt, refGenome)
    
    # homopolymer filter
    if isHomopolymer and 1.0 * MTCnt[origAlt] / usedMT < 0.99:
        fltr += 'HP;'

    # low complexity filter
    if isLowComplexity and 1.0 * MTCnt[origAlt] / usedMT < 0.99:
        fltr += 'LowC;'

    # strand bias and discordant pairs filter
    af_alt = 100.0 * alleleCnt[origAlt] / cvg
    pairs = discordPairCnt[origAlt] + concordPairCnt[origAlt] # total number of paired reads covering the pos
    if pairs >= 1000 and 1.0 * discordPairCnt[origAlt] / pairs >= 0.5:
        fltr += 'DP;'
    elif af_alt <= 60.0:
        refR = reverseCnt[origRef]
        refF = forwardCnt[origRef]
        altR = reverseCnt[origAlt]
        altF = forwardCnt[origAlt]
        fisher = scipy.stats.fisher_exact([[refR, refF], [altR, altF]])
        oddsRatio = fisher[0]
        pvalue = fisher[1]
        if pvalue < 0.00001 and (oddsRatio >= 50 or oddsRatio <= 1.0/50):
            fltr += 'SB;'

    # base quality filter. Reject if more than 40% reads are lowQ
    if vtype == 'SNP' and origAlt in alleleCnt.keys() and origAlt in lowQReads.keys():
        bqAlt = 1.0 * lowQReads[origAlt] / alleleCnt[origAlt] 
    else:
        bqAlt = 0.0
    if bqAlt > 0.4:
        fltr += 'LowQ;'

    # random end and fixed end position filters
    if vtype == 'SNP':
        # random end position filter
        endBase = 20  # distance to barcode end of the read
        # R1
        refLeEnd = sum(d <= endBase for d in r1BcEndPos[origRef])  # number of REF R2 reads with distance <= endBase
        refGtEnd = len(r1BcEndPos[origRef]) - refLeEnd           # number of REF R2 reads with distance > endBase
        altLeEnd = sum(d <= endBase for d in r1BcEndPos[origAlt])  # number of ALT R2 reads with distance <= endBase
        altGtEnd = len(r1BcEndPos[origAlt]) - altLeEnd           # number of ALT R2 reads with distance > endBase
        fisher = scipy.stats.fisher_exact([[refLeEnd, refGtEnd], [altLeEnd, altGtEnd]])
        oddsRatio = fisher[0]
        pvalue = fisher[1]
        if pvalue < 0.001 and oddsRatio < 0.05 and af_alt <= 60.0:
            fltr += 'R1CP;'
        # R2
        refLeEnd = sum(d <= endBase for d in r2BcEndPos[origRef])  # number of REF R2 reads with distance <= endBase
        refGtEnd = len(r2BcEndPos[origRef]) - refLeEnd           # number of REF R2 reads with distance > endBase
        altLeEnd = sum(d <= endBase for d in r2BcEndPos[origAlt])  # number of ALT R2 reads with distance <= endBase
        altGtEnd = len(r2BcEndPos[origAlt]) - altLeEnd           # number of ALT R2 reads with distance > endBase
        fisher = scipy.stats.fisher_exact([[refLeEnd, refGtEnd], [altLeEnd, altGtEnd]])
        oddsRatio = fisher[0]
        pvalue = fisher[1]
        if pvalue < 0.001 and oddsRatio < 0.05 and af_alt <= 60.0:
            fltr += 'R2CP;'

        # fixed end position filter
        endBase = primerDist # distance to primer end of the read
        refLeEnd = sum(d <= endBase for d in r2PrimerEndPos[origRef])   # number of REF R2 reads with distance <= endBase
        refGtEnd = len(r2PrimerEndPos[origRef]) - refLeEnd            # number of REF R2 reads with distance > endBase
        altLeEnd = sum(d <= endBase for d in r2PrimerEndPos[origAlt])   # number of ALT R2 reads with distance <= endBase
        altGtEnd = len(r2PrimerEndPos[origAlt]) - altLeEnd            # number of ALT R2 reads with distance > endBase
        fisher = scipy.stats.fisher_exact([[refLeEnd, refGtEnd], [altLeEnd, altGtEnd]])
        oddsRatio = fisher[0]
        pvalue = fisher[1]
        # reject if variant is clustered within 2 bases from primer sequence due to possible enzyme initiation error
        if altLeEnd + altGtEnd > 0:
            if 1.0 * altLeEnd / (altLeEnd + altGtEnd) >= 0.98 or (pvalue < 0.001 and oddsRatio < 1.0/20):
                fltr += 'PrimerCP;'

    # done
    return fltr
    
#-------------------------------------------------------------------------------------
# function to call variants
#-------------------------------------------------------------------------------------
def vc(bamFile, chrom, pos, minBQ, minMQ, mtDepth, rpb, hpLen,
       mismatchThr, mtDrop, maxMT, primerDist, refGenome):

    samfile = pysam.AlignmentFile(bamFile, 'rb')
    idx = 0
    cvg = 0
    bcDict = defaultdict(lambda: defaultdict(list)) 
    allBcDict = defaultdict(list) 
    alleleCnt = defaultdict(int)
    MTCnt = defaultdict(int)
    r0BcEndPos = defaultdict(list)
    r1BcEndPos = defaultdict(list)
    r2BcEndPos = defaultdict(list)
    r2PrimerEndPos = defaultdict(list)
    r0PrimerEndPos = defaultdict(list)
    MT3Cnt = 0
    MT5Cnt = 0
    MT7Cnt = 0
    MT10Cnt = 0
    strongMTCnt = defaultdict(int)
    predIndex = defaultdict(lambda: defaultdict(float))
    finalDict = defaultdict(float)
    r0Cnt = defaultdict(int)
    r1Cnt = defaultdict(int)
    r2Cnt = defaultdict(int)
    forwardCnt = defaultdict(int)
    reverseCnt = defaultdict(int)
    concordPairCnt = defaultdict(int)
    discordPairCnt = defaultdict(int)
    mismatchCnt = defaultdict(float)
    bqSum = defaultdict(int)
    lowQReads = defaultdict(int)

    # set threshold for strongMT based on mtDepth
    if rpb < 1.5:
        smt = 2.0
    elif rpb < 3.0:
        smt = 3.0
    else:
        smt = 4.0
    print 'smt:{}'.format(smt)

    # get reference base
    refseq = pysam.FastaFile(refGenome)
    origRef = refseq.fetch(reference=chrom, start=int(pos)-1, end=int(pos))
    origRef = origRef.upper()

    print 'target chromosome:{}'.format(chrom)
    print 'target position:{}'.format(pos)
    print 'reference base:{}'.format(origRef)

    # pile up reads
    for read in samfile.pileup(region = chrom + ':' + pos + ':' + pos,
        truncate=True, max_depth=1000000, stepper='nofilter'):

        print '\ncoverage at base:{}, number of reads:{}'.format(read.pos, read.n)

        for pileupRead in read.pileups:
            # read ID
            qname = pileupRead.alignment.query_name
            qnameSplit = qname.split(":")
            readid = ':'.join(qnameSplit[:-2])

            print '\n\n*QNAME:{}'.format(qname)
            print 'Query sequence:{}'.format(pileupRead.alignment.query_sequence)

            # barcode sequence
            BC = qnameSplit[-2]
            print 'Barcode(BC):{}'.format(BC)
            # duplex tag - temporary hack from end of readid - should be CC, TT, or NN for duplex runs
            duplexTag = qnameSplit[-3]

            # mapping quality
            mq = pileupRead.alignment.mapping_quality
            print 'mapping quality(mq):{}'.format(mq)

            # get NM tag 
            NM = 0 
            allTags = pileupRead.alignment.tags
            for (tag, value) in allTags:
                if tag == 'NM':  # 
                    NM = value
                    break
            print 'Edit distance tag(NM):{}'.format(NM)
            # NM: number of changes necessary to make this equal the reference, excluding clipping

            # count number of INDELs in the read sequence
            nIndel = 0
            cigar = pileupRead.alignment.cigar
            print 'CIGAR:{}'.format(cigar)
            print 'CIGAR:{}'.format(pileupRead.alignment.cigarstring)

            cigarOrder = 1
            leftSP = 0  # soft clipped bases on the left
            rightSP = 0  # soft clipped bases on the right
            for (op, value) in cigar:
                # 1 for insertion
                if op == 1 or op == 2:
                    nIndel += value 
                if cigarOrder == 1 and op == 4:
                    leftSP = value
                if cigarOrder > 1 and op == 4:
                    rightSP += value
                cigarOrder += 1
            print 'nIndel:{}, cigarOrder:{}, leftSP:{}, rightSP:{}'.format(
                nIndel, cigarOrder, leftSP, rightSP)

            # Number of mismatches except INDEL, including softcilpped sequences 
            mismatch = max(0, NM - nIndel)
            # read length, including softclip
            readLen = pileupRead.alignment.query_length
            # calculate mismatch per 100 bases
            mismatchPer100b = 100.0 * mismatch / readLen if readLen > 0 else 0.0
            print 'mismatch:{}, readLen:{}, mismatch%:{}'.format(
                mismatch, readLen, mismatchPer100b)

            # paired read
            if pileupRead.alignment.is_read1:
                pairOrder = 'R1'
            if pileupRead.alignment.is_read2:
                pairOrder = 'R2'
            if (pileupRead.alignment.is_read1 is False) and (pileupRead.alignment.is_read2 is False):
                # in case of single-end reads 
                pairOrder = None

            # +/- strand
            strand = 'Reverse' if pileupRead.alignment.is_reverse else 'Forward'
            print 'is_read1:{}, is_read2:{}, final pairOrder:{}, strand:{}'.format(
                pileupRead.alignment.is_read1,
                pileupRead.alignment.is_read2,
                pairOrder,
                strand,
            )

            # coverage -- read, not fragment
            cvg += 1

            # 1. check if the site is the beginning of insertion
            print 'pileupRead.query_position:{}'.format(pileupRead.query_position)
            # PileupRead.indel: indel length for the position following the current pileup site
            print 'pileupRead.indel:{}'.format(pileupRead.indel)
            print 'pileupRead.is_del:{}'.format(pileupRead.is_del)
            print 'query alignment length:{}'.format(pileupRead.alignment.query_alignment_length)
            print 'query alignment sequence:{}'.format(pileupRead.alignment.query_alignment_sequence)

            if pileupRead.indel > 0:
                site = pileupRead.alignment.query_sequence[pileupRead.query_position]
                inserted = pileupRead.alignment.query_sequence[
                    (pileupRead.query_position + 1) : \
                    (pileupRead.query_position + 1 +  pileupRead.indel)]
                print 'site:{}, inserted:{}'.format(site, inserted)

                base = 'INS|' + site + '|' + site + inserted
                bq = pileupRead.alignment.query_qualities[pileupRead.query_position]
                bqSum[base] += bq
                # inclusion condition
                incCond = bq >= minBQ and mq >= minMQ and mismatchPer100b <= mismatchThr

                alleleCnt[base] += 1
                mismatchCnt[base] += mismatchPer100b
                if pairOrder == 'R1':
                    r1Cnt[base] += 1
                elif pairOrder == 'R2':
                    r2Cnt[base] += 1
                else:
                    r0Cnt[base] += 1 
                if strand == 'Reverse':
                    reverseCnt[base] += 1
                else:
                    forwardCnt[base] += 1

            # 2. check if the site is the beginning of deletion
            elif pileupRead.indel < 0:
                # get alt base
                site = pileupRead.alignment.query_sequence[pileupRead.query_position]
                # get ref base
                deleted = refseq.fetch(
                    reference=chrom, start=int(pos),
                    end=int(pos) + abs(pileupRead.indel)
                )
                deleted = deleted.upper()
                base = 'DEL|' + site + deleted + '|' + site
                bq = pileupRead.alignment.query_qualities[pileupRead.query_position]
                bqSum[base] += bq

                # inclusion condition
                incCond = bq >= minBQ and mq >= minMQ and mismatchPer100b <= mismatchThr
                alleleCnt[base] += 1
                mismatchCnt[base] += mismatchPer100b
                if pairOrder == 'R1':
                    r1Cnt[base] += 1
                elif pairOrder == 'R2':
                    r2Cnt[base] += 1
                else:
                    r0Cnt[base] += 1
                if strand == 'Reverse':
                    reverseCnt[base] += 1
                else:
                    forwardCnt[base] += 1
                
            # 3. site is not beginning of any INDEL
            else:
                # If the site ifself is a deletion, set quality = minBQ 
                if pileupRead.is_del:
                    base = 'DEL'
                    bq = minBQ
                    bqSum[base] += bq
                    # inclusion condition
                    incCond = bq >= minBQ and mq >= minMQ and mismatchPer100b <= mismatchThr

                # if the site is a regular locus, 
                else: 
                    base = pileupRead.alignment.query_sequence[pileupRead.query_position]
                    # note: query_sequence includes soft clipped bases
                    bq = pileupRead.alignment.query_qualities[pileupRead.query_position]
                    bqSum[base] += bq
                    # count the number of low quality reads (less than Q20 by default) for each base
                    if bq < minBQ:
                        lowQReads[base] += 1

                    # inclusion condition
                    incCond = bq >= minBQ and mq >= minMQ and mismatchPer100b <= mismatchThr

                    if pairOrder == 'R1':
                        # distance to the barcode end in R1; 
                        if strand == 'Reverse':
                            distToBcEnd = pileupRead.alignment.query_alignment_length - (pileupRead.query_position - leftSP)
                        else:
                            distToBcEnd = pileupRead.query_position - leftSP
                        if incCond:
                            r1BcEndPos[base].append(distToBcEnd)
                        r1Cnt[base] += 1
                    elif pairOrder == 'R2':
                        # distance to the barcode and/or primer end in R2. Different cases for forward and reverse strand
                        if strand == 'Reverse':
                            distToBcEnd = pileupRead.query_position - leftSP
                            distToPrimerEnd = pileupRead.alignment.query_alignment_length - (pileupRead.query_position - leftSP)
                        else:
                            distToBcEnd = pileupRead.alignment.query_alignment_length - (pileupRead.query_position - leftSP)
                            distToPrimerEnd = pileupRead.query_position - leftSP
                        if incCond:
                            r2BcEndPos[base].append(distToBcEnd)
                            r2PrimerEndPos[base].append(distToPrimerEnd)
                        r2Cnt[base] += 1
                    else:  # same as R2
                        if strand == 'Reverse':
                            distToBcEnd = pileupRead.query_position - leftSP
                            distToPrimerEnd = pileupRead.alignment.query_alignment_length - (pileupRead.query_position - leftSP)
                        else:
                            distToBcEnd = pileupRead.alignment.query_alignment_length - (pileupRead.query_position - leftSP)
                            distToPrimerEnd = pileupRead.query_position - leftSP
                        if incCond:
                            r0BcEndPos[base].append(distToBcEnd)
                            r0PrimerEndPos[base].append(distToPrimerEnd)
                        r0Cnt[base] += 1

                    if strand == 'Reverse':
                        reverseCnt[base] += 1
                    else:
                        forwardCnt[base] += 1

                    print 'distToBcEnd:{}'.format(distToBcEnd)
                    # print 'r2BcEndPos:{}'.format(r2BcEndPos)
                    # print 'r2PrimerEndPos:{}'.format(r2PrimerEndPos)
                    print 'r0BcEndPos:{}'.format(r0BcEndPos)
                    print 'r0PrimerEndPos:{}'.format(r0PrimerEndPos)

                alleleCnt[base] += 1
                mismatchCnt[base] += mismatchPer100b

            print 'base:{}'.format(base)
            print 'is_reverse:{}'.format(pileupRead.alignment.is_reverse)
            print 'strand:{}'.format(strand)

            print '=== incCond:{} ==='.format(incCond)
            '''
            print '\t--bq >= minBQ ? {} >= {}, {}'.format(bq, minBQ, bq >= minBQ)
            print '\t--mq >= minMQ ? {} >= {}, {}'.format(mq, minMQ, mq >= minMQ)
            print '\tmismatch% <= mismatchThr ? {} <= {}, {}'.format(
                mismatchPer100b, mismatchThr, mismatchPer100b <= mismatchThr)
            '''
            # print 'bqSum:{}'.format(bqSum)
            # print 'alleleCnt:{}'.format(alleleCnt)
            # print 'mismatchCnt:{}'.format(mismatchCnt)
            # print 'r1Cnt:{}, r2Cnt:{}, r0Cnt:{}'.format(r1Cnt, r2Cnt, r0Cnt)
            # print 'reverseCnt:{}'.format(reverseCnt)
            # print 'forwardCnt:{}'.format(forwardCnt)

            # count total number of fragments and MTs
            if readid not in allBcDict[BC]:
                allBcDict[BC].append(readid)

            # decide which read goes into analysis
            if incCond:
                # new read id
                if readid not in bcDict[BC]:
                    prob = pow(10.0, -bq / 10.0)  # base-calling error probabolities, P
                    readinfo = [base, prob, pairOrder]
                    bcDict[BC][readid].append(readinfo)

                # same read id (for paired read, max prob will be remain b/w 2 reads)
                elif base == bcDict[BC][readid][0][0] or base in ['N', '*']:
                    bcDict[BC][readid][0][1] = \
                        max((pow(10.0, -bq / 10.0) , bcDict[BC][readid][0][1]))
                    bcDict[BC][readid][0][2] = 'Paired'
                    if base == bcDict[BC][readid][0][0]:
                        concordPairCnt[base] += 1

                # same read id but different base
                else:
                    del bcDict[BC][readid]
                    discordPairCnt[base] += 1
            else:
                continue

            print '=== bcDict ===\n{}'.format(bcDict)
            print '-'*60
        

        print 'at position {}:{} result---------'.format(chrom, pos)
        print 'Total Barcodes:{}'.format(len(allBcDict))
        print 'Passed Barcodes:{}'.format(len(bcDict))
        print '=' * 60

    print 'final BcDict------'
    print 'allBcDict({}):{}'.format(len(allBcDict), allBcDict.keys())
    print 'bcDict({}):{}'.format(len(bcDict), bcDict.keys())
    
    '''
    allBcDict = {
        'chr10-0-43606627-CAAAACGCAATA': ['3LMDR:00595:00544'],
        'chr10-1-43606567-TATTGCGTTTTG': ['3LMDR:00318:00485',
    }                                     '3LMDR:01078:01158'],
    bcDict = {
        'chr10-0-43606627-CAAAACGCAATA': {
            '3LMDR:00595:00544': [['A', 0.00031622776601683794, None]]},
        'chr10-1-43606567-TATTGCGTTTTG': {
            '3LMDR:00839:00334': [['G', 0.0012589254117941675, None]],
            '3LMDR:01078:01158': [['G', 0.0025118864315095794, None]],
            '3LMDR:00318:00485': [['G', 0.00039810717055349735, None]],
            '3LMDR:01099:01104': [['A', 0.0007943282347242813, None]]}
    }
    '''
    ##### endfor ####

    # total number of MT, fragments, reads, including those dropped from analysis
    allMT = len(allBcDict)  # number of barcodes
    allFrag = sum([len(allBcDict[bc]) for bc in allBcDict])  # number of reads
    print 'barcodes#(allMt):{}, reads#(allFrag):{}'.format(allMT, allFrag)

    # downsampling MTs (not dropped) to args.maxMT

    # maxMT: randomly downsample to X MTs (max number of MTs at any position).
    # if set to 0 (default), maxMT = 2.0 * mean MT depth
    ds = maxMT if maxMT > 0 else int(round(2.0 * mtDepth))

    # MTs used
    usedMT = min(ds, len(bcDict))  # choose minimum barcode number
    print 'downsample(ds, maxMT):{}'.format(ds)
    print 'MTs used:{}'.format(usedMT)

    # done if zero coverage (note hack for 41 blank output fields!)
    # bcDict has no reads
    if usedMT == 0:  # no reads
        out_long = '\t'.join([chrom, pos, origRef] + ['']*41 + ['Zero_Coverage'])
        return out_long
        
    if len(bcDict) > ds:
        random.seed(pos)
        bcKeys = random.sample(bcDict.keys(), ds)  # choose random X MTs among barcodes (X=ds)
        print 'downsampled random barcodes:{}'.format(bcKeys)
    else:  # no need to downsample
        bcKeys = bcDict.keys()
    usedFrag = sum([len(bcDict[bc]) for bc in bcKeys])  # number of downsampled reads
    print 'downsampled reads#(usedFrag):{}'.format(usedFrag)

    totalR1 = sum(r1Cnt.values())
    totalR2 = sum(r2Cnt.values())
    totlaR0 = sum(r0Cnt.values())
    print 'totlaR0:{}'.format(totlaR0)

    for bc in bcKeys:
        print "P to Q, and find max base for {}".format(bc)

        bcProb = calProb(bcDict[bc], mtDrop)

        '''
        bcProb = {
             'A': 0.0009066429730636374,
             'INS|G|GA': 0.9972800710808093,
             'G': 0.0009066429730636374,
             'T': 0.0009066429730636374
        }
        '''
        # Probability to QualityScore(Index) by barcode
        for char in bcProb.iterkeys():
            x = 1.0 - bcProb[char]
            log10P = -math.log10(x) if x > 0.0 else 16.0
            predIndex[bc][char] = log10P
            finalDict[char] += log10P

        max_base = [x for x in predIndex[bc].keys() if predIndex[bc][x] == max(predIndex[bc].values())]
        print 'max_base:{} in barcode({})'.format(max_base, bc)
        # max_base: ['G']  the base which has max log10 value in one barcode
        
        if len(max_base) == 1:
            cons = max_base[0]
            print "MTCnt + 1 !! now for base {}".format(cons)
            MTCnt[cons] += 1
            if predIndex[bc][cons] > smt:  # check threshold < current index value
                strongMTCnt[cons] += 1

        # Tie in max predIndex is most likely due to single read MT. 

        if len(bcDict[bc]) == 1:  # only one read in barcode
            cons = bcDict[bc].values()[0][0][0]
            print "MTCnt + 1 !! now for base {}".format(cons)
            MTCnt[cons] += 1

        if len(bcDict[bc]) >= 3:
            MT3Cnt += 1
        if len(bcDict[bc]) >= 5:
            MT5Cnt += 1
        if len(bcDict[bc]) >= 7:
            MT7Cnt += 1
        if len(bcDict[bc]) >= 10:
            MT10Cnt += 1

    print 'predIndex:{}'.format(predIndex)
    print 'finalDict:{}'.format(finalDict)
    print 'strongMTCnt:{}'.format(strongMTCnt)
    print 'MTCnt:{}'.format(MTCnt)

    sortedList = sorted(finalDict.items(), key=operator.itemgetter(1), reverse=True)
    print 'sortedList:{}'.format(sortedList)
    maxBase = sortedList[0][0]
    maxPI = sortedList[0][1]
    secondMaxBase = sortedList[1][0]
    secondMaxPI = sortedList[1][1]
    
    # call variants
    origAlt = secondMaxBase if maxBase == origRef else maxBase
    altPI = secondMaxPI if maxBase == origRef else maxPI
    print 'call variants------'
    print 'origRef:{}, origAlt:{}, altPI:{}'.format(origRef, origAlt, altPI)

    # convert from internal smCounter format to format needed for output
    (ref, alt, vtype) = convertToVcf(origRef,origAlt)
    print 'VCF line:{}'.format((ref, alt, vtype))

    # apply filters if PI >= 5 (at least 2 MTs), and locus not in a deletion
    fltr = ';'
    if altPI >= 5 and vtype in ('SNP', 'INDEL'):
        fltr = filterVariants(ref,alt,vtype,origAlt,origRef,usedMT,strongMTCnt,chrom,pos,hpLen,refGenome,MTCnt,alleleCnt,cvg,discordPairCnt,concordPairCnt,reverseCnt,forwardCnt,lowQReads,r1BcEndPos,r2BcEndPos,r2PrimerEndPos,primerDist)
        
    # identify possible bi-allelic variants - top 2 alleles are non-reference and both VMFs >= 45%. Not necessarily passing the filters. 
    mfAlt = 1.0 * MTCnt[maxBase] / usedMT  # MT fraction of the base with the highest PI
    mfAlt2 = 1.0 * MTCnt[secondMaxBase] / usedMT  # MT fraction of the base with the second highest PI

    # conditions to be considered bi-allelic
    if maxBase != origRef and secondMaxBase != origRef and mfAlt >= 0.45 and mfAlt2 >= 0.45:
    
        # convert from internal smCounter format to format needed for output
        origAlt2 = secondMaxBase
        (ref2, alt2, vtype2) = convertToVcf(origRef,origAlt2)

        # apply filters to 2nd variant if PI2 >= 5 (at least 2 MTs), and locus not in a deletion
        fltr2 = ';'
        if secondMaxPI >= 5 and vtype2 in ('SNP', 'INDEL'):
            fltr2 = filterVariants(ref2,alt2,vtype2,origAlt2,origRef,usedMT,strongMTCnt,chrom,pos,hpLen,refGenome,MTCnt,alleleCnt,cvg,discordPairCnt,concordPairCnt,reverseCnt,forwardCnt,lowQReads,r1BcEndPos,r2BcEndPos,r2PrimerEndPos,primerDist)

        # prepare output for bi-allelic variants (if var2 is filtered, regardless of var1, do nothing. output var1 only)
        if fltr == ';' and fltr2 == ';': # both var1 and var2 pass the filters -- this is a bi-allelic variant. var1's statistics (MT, DP, etc) are reported
            alt = alt + ',' + alt2
            vtype = vtype.lower() + ',' + vtype2.lower()
        elif fltr != ';' and fltr2 == ';':   # if var1 is filtered and the var2 passes, then it's a single variant of var2
            alt = alt2
            fltr = fltr2
            origAlt = origAlt2

    # build detailed output vector
    print 'build detailed output-------'
    print 'alleleCnt:{}'.format(alleleCnt)
    print 'coverage(cvg):{}'.format(cvg)

    frac_alt = round((1.0 * alleleCnt[origAlt] / cvg),4)     # based on all reads, including the excluded reads
    frac_A  = round((1.0 * alleleCnt['A'] / cvg),4)
    frac_T  = round((1.0 * alleleCnt['T'] / cvg),4) 
    frac_G  = round((1.0 * alleleCnt['G'] / cvg),4) 
    frac_C  = round((1.0 * alleleCnt['C'] / cvg),4)
    fracs = (alleleCnt['A'], alleleCnt['T'], alleleCnt['G'], alleleCnt['C'], frac_A, frac_T, frac_G, frac_C)
    
    print 'MTCnt:{}'.format(MTCnt)
    print 'usedMT:{}'.format(usedMT)

    MT_f_alt = round((1.0 * MTCnt[origAlt] / usedMT),4)  # based on only used MTs
    MT_f_A  = round((1.0 * MTCnt['A'] / usedMT),4)
    MT_f_T  = round((1.0 * MTCnt['T'] / usedMT),4)
    MT_f_G  = round((1.0 * MTCnt['G'] / usedMT),4)
    MT_f_C  = round((1.0 * MTCnt['C'] / usedMT),4)
    MTs = (MT3Cnt, MT5Cnt, MT7Cnt, MT10Cnt, MTCnt['A'], MTCnt['T'], MTCnt['G'], MTCnt['C'], MT_f_A, MT_f_T, MT_f_G, MT_f_C)
    
    strongMT = (strongMTCnt['A'], strongMTCnt['T'], strongMTCnt['G'], strongMTCnt['C'])
    predIdx = (round(finalDict['A'], 2), round(finalDict['T'], 2), round(finalDict['G'], 2), round(finalDict['C'], 2))

    outvec = [chrom, pos, ref, alt, vtype, cvg, allFrag, allMT, usedFrag, usedMT, round(finalDict[origAlt], 2), alleleCnt[origAlt], frac_alt, MTCnt[origAlt], MT_f_alt, strongMTCnt[origAlt]]
    outvec.extend(fracs)
    outvec.extend(MTs)
    outvec.extend(strongMT)
    outvec.extend(predIdx)
    outvec.append(fltr)
    out_long = '\t'.join((str(x) for x in outvec))
    print 'out_log====>{}'.format(out_long)
    return out_long

#------------------------------------------------------------------------------------------------
# wrapper function for "vc()" - because Python multiprocessing module does not pass stack trace
#------------------------------------------------------------------------------------------------
def vc_wrapper(*args):
     output = vc(*args)
     return output

#------------------------------------------------------------------------------------------------
# global for argument parsing (hack that works when calling from either command line or pipeline)
#------------------------------------------------------------------------------------------------
parser = None
def argParseInit(): # this is done inside a function because multiprocessing module imports the script
    global parser
    parser = argparse.ArgumentParser(description='Variant calling using molecular barcodes', fromfile_prefix_chars='@')
    parser.add_argument('--outPrefix', default=None, required=True, help='prefix for output files')
    parser.add_argument('--bamFile', default=None, required=True, help='BAM file')
    parser.add_argument('--bedTarget', default=None, required=True, help='BED file for target region')
    parser.add_argument('--mtDepth', default=None, required=True, type=int, help='Mean MT depth')
    parser.add_argument('--rpb', default=None, required=True, type=float, help='Mean read pairs per MT')
    parser.add_argument('--nCPU', type=int, default=1 , help='number of CPUs to use in parallel')
    parser.add_argument('--minBQ', type=int, default=20, help='minimum base quality allowed for analysis')
    parser.add_argument('--minMQ', type=int, default=30, help='minimum mapping quality allowed for analysis')
    parser.add_argument('--hpLen', type=int, default=10, help='Minimum length for homopolymers')
    parser.add_argument('--mismatchThr', type=float, default=6.0, help='average number of mismatches per 100 bases allowed')
    parser.add_argument('--mtDrop', type=int, default=0, help='Drop MTs with lower than or equal to X reads.')
    parser.add_argument('--maxMT', type=int, default=0, help='Randomly downsample to X MTs (max number of MTs at any position). If set to 0 (default), maxMT = 2.0 * mean MT depth')
    parser.add_argument('--primerDist', type=int, default=2, help='filter variants that are within X bases to primer')
    parser.add_argument('--threshold', type=int, default=0, help='Minimum prediction index for a variant to be called. Must be non-negative. Typically ranges from 10 to 60. If set to 0 (default), smCounter will choose the appropriate cutoff based on the mean MT depth.')
    parser.add_argument('--refGenome', default = '/qgen/home/rvijaya/downloads/alt_hap_masked_ref/ucsc.hg19.fasta')
    parser.add_argument('--bedTandemRepeats', default = None, required=False, help = 'bed for UCSC tandem repeats')
    parser.add_argument('--bedRepeatMaskerSubset', default = None, required=False, help = 'bed for RepeatMasker simple repeats, low complexity, microsatellite regions')
    parser.add_argument('--bedtoolsPath', default = '/qgen/bin/bedtools-2.25.0/bin/', help = 'path to bedtools')
    parser.add_argument('--runPath', default=None, help='path to working directory')
    parser.add_argument('--logFile', default=None, help='log file')
    parser.add_argument('--paramFile', default=None, help='optional parameter file that contains the above paramters. if specified, this must be the only parameter, except for --logFile.')
    
#--------------------------------------------------------------------------------------
# main function
#--------------------------------------------------------------------------------------
def main(args):
    # log run start
    timeStart = datetime.datetime.now()
    print "smCounter started at " + str(timeStart)

    # if argument parser global not assigned yet, initialize it
    if parser == None:
        argParseInit()
    
    # get arguments passed in via a lambda object (e.g. from upstream pipeline)
    if type(args) is not argparse.Namespace:
        argsList = []
        for argName, argVal in args.iteritems():
            argsList.append("--{0}={1}".format(argName, argVal))
        args = parser.parse_args(argsList)
        
    # get arguments from disk file specified on command line (warning: this silently deletes all actual command line parameters)
    elif args.paramFile != None:
        args = parser.parse_args(("@" + args.paramFile,))
        
    # echo all parameters to the log file
    for argName, argVal in vars(args).iteritems():
        print(argName, argVal)
     
    # change working directory to runDir
    if args.runPath != None:
        os.chdir(args.runPath)

    # make list of loci to call variants
    locList = []
    for line in open(args.bedTarget, 'r'):
        if not line.startswith("track "):
            (chrom, regionStart, regionEnd) = line.strip().split('\t')[0:3]
            for pos in range(int(regionStart), int(regionEnd)):
                locList.append((chrom, str(pos+1)))

    print 'locList ----> {}'.format(locList)
    pool = multiprocessing.Pool(processes=args.nCPU)
    results = [
        pool.apply_async(
            vc_wrapper,
            args=(
                args.bamFile, x[0], x[1], args.minBQ, args.minMQ,
                args.mtDepth, args.rpb, args.hpLen, args.mismatchThr,
                args.mtDrop, args.maxMT, args.primerDist, args.refGenome
            )
         ) for x in locList
    ]
    output = [p.get() for p in results]

    pool.close()
    pool.join()

    # check for exceptions thrown by vc()
    for idx in range(len(output)):
        line = output[idx]
        if line.startswith("Exception thrown!"):
            raise Exception("Exception thrown in vc() at location: " + str(locList[idx]))
    
    # report start of variant filtering
    print("begin variant filtering and output")
    
    # merge and sort RepeatMasker tracks (could be done prior to run)  Note: assuming TRF repeat already merged and esorted!!
    bedExe = args.bedtoolsPath + 'bedtools'
    if args.bedRepeatMaskerSubset is not None:
        bedRepeatMasker = args.outPrefix + '.tmp.repeatMasker.bed'
        cmd = ('{bedExe} merge -c 4 -o distinct -i {repeatbed} | '
               '{bedExe} sort -i - > {out}').format(
                bedExe=bedExe,
                repeatbed=args.bedRepeatMaskerSubset,
                out=bedRepeatMasker
              )
        subprocess.check_call(cmd, shell=True)

    # merge and sort target region
    bedTarget = args.outPrefix + '.tmp.target.bed'
    cmd = ('{bedExe} merge -i {targetbed} | '
           '{bedExe} sort -i - > {out}').format(
            bedExe=bedExe,
            targetbed=args.bedTarget,
            out=bedTarget,
          )
    subprocess.check_call(cmd, shell=True)
    
    # intersect 2 repeats tracks with target region (for tandem repeats)
    if args.bedTandemRepeats is not None:
        cmd = ('{bedExe} intersect -a {repeatbed} -b {targetbed} | '
               '{bedExe} sort -i - > {outprefix}.tmp.target.repeats1.bed').format(
               bedExe=bedExe,
               repeatbed=args.bedTandemRepeats,
               targetbed=bedTarget,
               outprefix=args.outPrefix,
              )
        subprocess.check_call(cmd, shell=True)

    # for simple repeats
    if args.bedRepeatMaskerSubset:
        cmd = ('{bedExe} intersect -a {repeatbed} -b {targetbed} | '
              '{bedExe} sort -i -> {outprefix}.tmp.target.repeats2.bed').format(
               bedExe=bedExe,
               repeatbed=bedRepeatMasker,
               targetbed=bedTarget,
               outprefix=args.outPrefix,
              )
        subprocess.check_call(cmd, shell=True)

    # read in tandem repeat list
    trfRegions = defaultdict(list)
    if args.bedTandemRepeats:
        for line in open(args.outPrefix + '.tmp.target.repeats1.bed', 'r'):
            vals = line.strip().split()
            (chrom, regionStart, regionEnd) = vals[0:3]
            trfRegions[chrom].append((int(regionStart), int(regionEnd), "RepT;"))

    # read in simple repeat, low complexity, satelite list
    rmRegions = defaultdict(list)
    if args.bedRepeatMaskerSubset:
        for line in open(args.outPrefix + '.tmp.target.repeats2.bed', 'r'):
            (chrom, regionStart, regionEnd, typeCodes) = line.strip().split()
            repTypes = []
            for typeCode in typeCodes.split(","):
                if typeCode == 'Simple_repeat':
                    repTypes.append('RepS')
                elif typeCode == 'Low_complexity':
                    repTypes.append('LowC')
                elif typeCode == 'Satellite':
                    repTypes.append('SL')
                else:
                    repTypes.append('Other_Repeat')
            repType = ";".join(repTypes) + ";"
            rmRegions[chrom].append((int(regionStart), int(regionEnd), repType))
        
    '''
    # remove intermediate files
    os.remove(args.outPrefix + '.tmp.target.bed')
    os.remove(args.outPrefix + '.tmp.repeatMasker.bed')
    os.remove(args.outPrefix + '.tmp.target.repeats1.bed')
    os.remove(args.outPrefix + '.tmp.target.repeats2.bed')
    '''

    # set up header columns (Note: "headerAll" must parallel the output of the vc() function.)
    headerAll = (
        'CHROM', 'POS', 'REF', 'ALT', 'TYPE', 'DP', 'FR', 'MT', 'UFR',
        'UMT', 'PI', 'VDP', 'VAF', 'VMT', 'VMF', 'VSM', 'DP_A', 'DP_T',
        'DP_G', 'DP_C', 'AF_A', 'AF_T', 'AF_G', 'AF_C', 'MT_3RPM',
        'MT_5RPM', 'MT_7RPM', 'MT_10RPM', 'UMT_A', 'UMT_T', 'UMT_G',
        'UMT_C', 'UMF_A', 'UMF_T', 'UMF_G', 'UMF_C', 'VSM_A', 'VSM_T',
        'VSM_G', 'VSM_C', 'PI_A', 'PI_T', 'PI_G', 'PI_C', 'FILTER'
    )
    headerVariants = (
        'CHROM', 'POS', 'REF', 'ALT', 'TYPE', 'DP', 'MT', 'UMT',
        'PI', 'THR', 'VMT', 'VMF', 'VSM', 'FILTER'
    )

    # set up hash of variable fields
    headerAllIndex = {}
    for i in range(len(headerAll)):
        headerAllIndex[headerAll[i]] = i
    print 'headerAllIndex--->',  headerAllIndex

    # ALL repeats filter. If MT fraction < 40% and the variant is inside the tandem repeat region, reject.
    for i in range(len(output)):
        outline = output[i]
        lineList = outline.split('\t')
        chromTr = lineList[headerAllIndex['CHROM']]
        altTr = lineList[headerAllIndex['ALT']]
        try:
            posTr = int(lineList[headerAllIndex['POS']])
        except ValueError:
            continue
        try:
            altMtFracTr = float(lineList[headerAllIndex['VMF']])
        except ValueError:
            continue
        try:
            pred = int(float(lineList[headerAllIndex['PI']]))
        except ValueError:
            pred = 0
        
        if pred >= 5 and altTr != 'DEL':
            # check tandem repeat from TRF if MT fraction < 40%
            if altMtFracTr < 40:
                for (locL, locR, repType) in trfRegions[chromTr]:
                    if locL < posTr <= locR:
                        lineList[-1] += repType
                        break

            # check simple repeat, lc, sl from RepeatMasker
            for (locL, locR, repType) in rmRegions[chromTr]:
                if locL < posTr <= locR:
                    lineList[-1] += repType
                    break
        
        lineList[-1] = 'PASS' if lineList[-1] == ';' else lineList[-1].strip(';')
        output[i] = '\t'.join(lineList)

    # VCF header
    header_vcf = \
        '##fileformat=VCFv4.2\n' + \
        '##reference=GRCh37\n' + \
        '##INFO=<ID=TYPE,Number=1,Type=String,Description="Variant type: SNP or INDEL">\n' + \
        '##INFO=<ID=DP,Number=1,Type=Integer,Description="Total read depth">\n' + \
        '##INFO=<ID=MT,Number=1,Type=Integer,Description="Total MT depth">\n' + \
        '##INFO=<ID=UMT,Number=1,Type=Integer,Description="Filtered MT depth">\n' + \
        '##INFO=<ID=PI,Number=1,Type=Float,Description="Variant prediction index">\n' + \
        '##INFO=<ID=THR,Number=1,Type=Integer,Description="Variant prediction index minimum threshold">\n' + \
        '##INFO=<ID=VMT,Number=1,Type=Integer,Description="Variant MT depth">\n' + \
        '##INFO=<ID=VMF,Number=1,Type=Float,Description="Variant MT fraction">\n' + \
        '##INFO=<ID=VSM,Number=1,Type=Integer,Description="Variant strong MT depth">\n' + \
        '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n' + \
        '##FORMAT=<ID=AD,Number=.,Type=Integer,Description="Filtered allelic MT depths for the ref and alt alleles">\n' + \
        '##FORMAT=<ID=VF,Number=1,Type=Float,Description="Variant MT fraction, same as VMF">\n' + \
        '##FILTER=<ID=RepT,Description="Variant in simple tandem repeat region, as defined by Tandem Repeats Finder">\n' + \
        '##FILTER=<ID=RepS,Description="Variant in simple repeat region, as defined by RepeatMasker">\n' + \
        '##FILTER=<ID=LowC,Description="Variant in low complexity region, as defined by RepeatMasker">\n' + \
        '##FILTER=<ID=SL,Description="Variant in micro-satelite region, as defined by RepeatMasker">\n' + \
        '##FILTER=<ID=HP,Description="Inside or flanked by homopolymer region">\n' + \
        '##FILTER=<ID=LM,Description="Low coverage (fewer than 5 MTs)">\n' + \
        '##FILTER=<ID=LSM,Description="Fewer than 2 strong MTs">\n' + \
        '##FILTER=<ID=SB,Description="Strand bias">\n' + \
        '##FILTER=<ID=LowQ,Description="Low base quality (mean < 22)">\n' + \
        '##FILTER=<ID=MM,Description="Too many genome reference mismatches in reads (default threshold is 6.5 per 100 bases)">\n' + \
        '##FILTER=<ID=DP,Description="Too many discordant read pairs">\n' + \
        '##FILTER=<ID=R1CP,Description="Variants are clustered at the end of R1 reads">\n' + \
        '##FILTER=<ID=R2CP,Description="Variants are clustered at the end of R2 reads">\n' + \
        '##FILTER=<ID=PrimerCP,Description="Variants are clustered immediately after the primer, possible enzyme initiation error">\n' + \
        '\t'.join(('#CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO', 'FORMAT', args.outPrefix)) + '\n'

    # set cutoff value for about 20 FP/Mb
    threshold = int(math.ceil(14.0 + 0.012 * args.mtDepth)) if args.threshold == 0 else args.threshold

    # open output files
    outAll = open(args.outPrefix + '.smCounter.all.txt', 'w')
    outVariants = open(args.outPrefix + '.smCounter.cut.txt', 'w')
    outVcf = open(args.outPrefix + '.smCounter.cut.vcf', 'w')

    # write column headers
    outAll.write('\t'.join(headerAll) + '\n')
    outVariants.write('\t'.join(headerVariants) + '\n')
    outVcf.write(header_vcf)
    
    for line in output:
        # write to the detailed output
        outAll.write(line)
        outAll.write("\n")
        
        # unpack text fields
        fields = line.split('\t')
        
        # skip if no PI
        PI = fields[headerAllIndex['PI']]
        if len(PI) == 0:
            continue
        
        # get ALT and prediction index
        ALT = fields[headerAllIndex['ALT']]
        QUAL = str(int(float(PI))) # truncate PI to conform to VCF phred-like tradition

        # write to vcf file and short output
        if int(QUAL) >= threshold and ALT != 'DEL': # if PI > threshold, write to vcf (regardless of filters)
        
            # parse fields needed from main data vector
            CHROM = fields[headerAllIndex['CHROM']]
            POS = fields[headerAllIndex['POS']]
            REF = fields[headerAllIndex['REF']]
            TYPE = fields[headerAllIndex['TYPE']]
            DP = fields[headerAllIndex['DP']]
            MT = fields[headerAllIndex['MT']]
            UMT = fields[headerAllIndex['UMT']]
            VMT = fields[headerAllIndex['VMT']]
            VMF = fields[headerAllIndex['VMF']]
            VSM = fields[headerAllIndex['VSM']]
            FILTER= fields[headerAllIndex['FILTER']]
            THR = str(threshold)
            INFO = ';'.join(('TYPE='+TYPE, 'DP='+DP, 'MT='+MT, 'UMT='+UMT, 'PI='+PI, 'THR='+THR, 'VMT='+VMT, 'VMF='+VMF, 'VSM='+VSM))
            
            # hack attempt to satisfy downstream software - not correct for germline heterozygous, male X, etc, etc, etc 
            alts = ALT.split(",")
            if len(alts) == 2:
                genotype = '1/2'
            elif len(alts) != 1:
                raise Exception("error hacking genotype field for " + alts)
            elif CHROM == "chrY" or CHROM == "chrM":
                genotype = '1'
            elif float(VMF) > 0.95:
                genotype = '1/1'
            else:
                genotype = '0/1'
            REFMT = str(int(UMT) - int(VMT))
            AD = REFMT + "," + VMT
            if len(alts) == 2:
                AD = AD + ",1"  # horrific hack for the 2nd alt

            # output
            FORMAT = 'GT:AD:VF'
            SAMPLE = ":".join((genotype,AD,VMF))
            ID = '.'
            vcfLine = '\t'.join((CHROM, POS, ID, REF, ALT, QUAL, FILTER, INFO, FORMAT, SAMPLE)) + '\n'
            shortLine = '\t'.join((CHROM, POS, REF, ALT, TYPE, DP, MT, UMT, PI, THR, VMT, VMF, VSM, FILTER)) + '\n' 
            outVcf.write(vcfLine)
            outVariants.write(shortLine)
            
            # debug counter for summary
            if TYPE == 'SNP':
                numCalledSnps = 0
            else:
                numCalledIndels = 0

        # parse fields needed from main data vector
        CHROM = fields[headerAllIndex['CHROM']]
        POS = fields[headerAllIndex['POS']]
        REF = fields[headerAllIndex['REF']]
        TYPE = fields[headerAllIndex['TYPE']]
        DP = fields[headerAllIndex['DP']]
        MT = fields[headerAllIndex['MT']]
        UMT = fields[headerAllIndex['UMT']]
        VMT = fields[headerAllIndex['VMT']]
        VMF = fields[headerAllIndex['VMF']]
        VSM = fields[headerAllIndex['VSM']]
        FILTER= fields[headerAllIndex['FILTER']]
        THR = str(threshold)
        INFO = ';'.join((
            'TYPE='+TYPE, 'DP='+DP, 'MT='+MT, 'UMT='+UMT, 'PI='+PI,
            'THR='+THR, 'VMT='+VMT, 'VMF='+VMF, 'VSM='+VSM
        ))
            
        # hack attempt to satisfy downstream software - not correct for germline heterozygous, male X, etc, etc, etc 
        alts = ALT.split(",")
        if len(alts) == 2:
            genotype = '1/2'
        elif len(alts) != 1:
            raise Exception("error hacking genotype field for " + alts)
        elif CHROM == "chrY" or CHROM == "chrM":
            genotype = '1'
        elif float(VMF) > 0.95:
            genotype = '1/1'
        else:
            genotype = '0/1'
        REFMT = str(int(UMT) - int(VMT))
        AD = REFMT + "," + VMT
        if len(alts) == 2:
            AD = AD + ",1"  # horrific hack for the 2nd alt

        # output
        FORMAT = 'GT:AD:VF'
        SAMPLE = ":".join((genotype,AD,VMF))
        ID = '.'
        vcfLine = '\t'.join((CHROM, POS, ID, REF, ALT, QUAL, FILTER, INFO, FORMAT, SAMPLE)) + '\n'
        shortLine = '\t'.join((CHROM, POS, REF, ALT, TYPE, DP, MT, UMT, PI, THR, VMT, VMF, VSM, FILTER)) + '\n' 
        outVcf.write(vcfLine)
        outVariants.write(shortLine)
          
        # debug counter for summary
        if TYPE == 'SNP':
            numCalledSnps    = 0
        else:
            numCalledIndels = 0

    outVcf.close()
    outAll.close()
    outVariants.close()
    
    # log run completion
    timeEnd = datetime.datetime.now()
    print("smCounter completed running at " + str(timeEnd))
    print("smCounter total time: "+ str(timeEnd-timeStart))
    
    # pass threshold back to caller
    return threshold

#----------------------------------------------------------------------------------------------
# pythonism to run from the command line
#----------------------------------------------------------------------------------------------
if __name__ == "__main__":
    # init the argumet parser
    argParseInit()
    
    # get command line arguments
    args = parser.parse_args()
    
    # initialize logger
    import run_log
    run_log.init(args.logFile)
    
    # call main program
    main(args)
