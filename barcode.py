"""Barcode-aware variant calling with smCounter.
"""
import os
import sys

import pandas as pd

# USAGE: barcode.py xxx.bam

def pairwise_alignment(target, real, min_dist=1):
    ''' compare 2 barcodes weather they are same or not.
        return "real" tag if they are same'''
    # print '==================='
    from Bio import pairwise2
    
    if abs(len(target) - len(real)) > min_dist:
        return None
    
    min_score = len(target) - min_dist if len(target) > len(real) else len(real) - min_dist
        
    # print 'Min score:', min_score    
    pairwise2.MAX_ALIGNMENTS = 10
    result = None
    for alignment in pairwise2.align.globalms(target, real, 1, -1, 0, 0):
        # print pairwise2.format_alignment(*alignment)
        if alignment[2] >= min_score:
            result = real
        else:
            continue
    return result

def parse_cigars_get_barcode_pos(cigar, flag):
    import re
    cigar_pattern = '([0-9]+[MIDNSHP])'
    cigar_search = re.compile(cigar_pattern)
    cigar_strings = cigar_search.findall(cigar)
    atomic_cigar_pattern = '([0-9]+)([MIDNSHP])'
    atomic_cigar_search = re.compile(atomic_cigar_pattern)

    if flag == 0: # +strand
        cigar_strings.reverse()
        
    barcode_pos = 0
    for cigar_string in cigar_strings:
        cigar_op_list = atomic_cigar_search.findall(cigar_string)[0]
        if cigar_op_list[1] == 'M':
            break
        else:
            barcode_pos += int(cigar_op_list[0])
            
    return (cigar_strings, barcode_pos)


def get_12n_barcode(row, n=12):
    # print row.name
    if row.FLAG not in [0, 16]:
        return None
    
    strings, barcode_pos = parse_cigars_get_barcode_pos(row.CIGAR, row.FLAG)
    
    # barcode is located upstream when read is - strand
    if barcode_pos == 0:
        return None
    barcode = row.SEQ[:barcode_pos] if row.FLAG == 16 else row.SEQ[-barcode_pos:]
    
    # cut barcode sequence 12 base
    if len(barcode) > 12:
        if row.FLAG == 0: # + strand
            barcode = barcode[:n]
        elif row.FLAG == 16: # - strand
            barcode = barcode[-n:]
            
    # check barcode length (default: > 7bp)
    minimum_barcode_length = n - 5
    if len(barcode) < minimum_barcode_length:
        return None
    return barcode


def compare_and_merge_barcode(target, barcodes_df, min_dist):
    ''' return matched barcodes in realgroup with given target '''
    if target.is_real is True:
        return target.name
    elif target.modified is not None:
        return target.modified
    
    modified = None
    for realb in barcodes_df[barcodes_df.is_real == True].itertuples():
        modified = pairwise_alignment(target.name, realb.Index, min_dist)
        if modified:
            ''' check the number of reads is large enough to merge '''
            if target.reads_no * 6 <= realb.reads_no:
                total_reads = target.reads_no + realb.reads_no
                barcodes_df.set_value(realb.Index, 'reads_no', total_reads)
                barcodes_df.set_value(target.name, 'reads_no', 0)
                break
            else:
                modified = None
    return modified


def get_new_barcode(row):
    original_barcode = row.BARCODE if row.BARCODE else 'undefined'
    clustered_barcode = row.FINAL_BARCODE if row.FINAL_BARCODE else original_barcode
    
    new_id = '{ori_id}:{chrno}-{strand}-{position}-{cluster_bar}:{ori_bar}'.format(
        ori_id=row.name,
        chrno=row.RNAME,
        strand='1' if row.FLAG == 16 else '0',
        position=row.POS,
        cluster_bar=clustered_barcode,
        ori_bar=original_barcode,
    )
    print(new_id)
    return new_id


def get_final_barcode(row, final_barcodes):
    if row.BARCODE is None:
        return None
    return final_barcodes.loc[row.BARCODE].modified


def _run_cluster_barcode(bamfile):
    prefix = os.path.splitext(bamfile)[0]
    samfile = "{}.sam".format(prefix)

    # convert bam to sam
    os.system("samtools view -G 4 {bamfile} > {samfile}".format(**locals()))
    
    # barcode extraction
    df = pd.read_table(
        samfile, index_col=0, header=None,
        usecols=list(range(0, 11)),
        names = ['QNAME', 'FLAG', 'RNAME', 'POS',
                 'MAPQ', 'CIGAR', 'RNEXT', 'PNEXT',
                 'TLEN', 'SEQ', 'QUAL'],
    )
    df['BARCODE'] = df.apply(get_12n_barcode, axis=1)
    df_dropna = df.dropna(subset=['BARCODE'])

    print 'total reads:', df.shape
    print 'total reads after dropna:', df_dropna.shape

    # grouping
    read_groups = df_dropna.sort_values(
        ['FLAG'], ascending=False).groupby(by=['BARCODE'])
    count = read_groups.agg('count')
    count = count.query('BARCODE != 0')
    count.to_csv("{prefix}_read_groups_sample_by_barcode.csv".format(**locals()))
    print("all groups:", count.index.size)
    print("1 read group:", count.query('FLAG == 1').index.size)
    print("over 1 read group:", count.query('FLAG > 1').index.size)
    print count.FLAG.describe()

    # Barcode clustering
    # step 1
    print(' --- step1 --- ')
    ''' sorting by number of reads in groups '''
    barcodes = count.sort_values(by='FLAG', ascending=False)[['FLAG']]
    barcodes.columns = ['reads_no']
    real_barcode_cutoff = int(barcodes.reads_no.max() * 0.05)
    barcodes['is_real'] = barcodes.reads_no >= real_barcode_cutoff
    barcodes.head()
    print 'real barcodes (reads >', real_barcode_cutoff, "):", \
          barcodes.is_real.sum(), '/', barcodes.index.size

    # Step 2
    print(' --- step2 --- ')
    barcodes['modified'] = None
    barcodes['modified'] = barcodes.apply(
        compare_and_merge_barcode, axis=1, barcodes_df=barcodes, min_dist=1)
    barcodes[:10]
    print 'total real barcodes:', barcodes.modified.count(), \
          '/', barcodes.index.size

    # Step 3
    print(' --- step3 --- ')
    barcodes_step3 = barcodes.copy()
    barcodes_step3.loc[
        (barcodes_step3.is_real == False) & \
        (barcodes_step3.reads_no >= 2), 'is_real'] = True
    print barcodes_step3.head()

    # Step 4
    print(' --- step4 --- ')
    barcodes_step4 = barcodes_step3.copy()
    barcodes_step4['modified'] = barcodes_step4.apply(
        compare_and_merge_barcode, axis=1,
        barcodes_df=barcodes_step4, min_dist=2
    )
    barcodes_step4[:10]
    print 'total real barcodes:', barcodes_step4.modified.count(), \
          '/', barcodes_step4.index.size

    # Summary
    a = barcodes_step4.query('is_real == False and reads_no > 0').index.size
    print('undefined barcodes: ', a)
    b = barcodes_step4.query('is_real == True').index.size
    print('real barcodes:', b)
    c = barcodes_step4.query('is_real == False and reads_no == 0').index.size
    print('merged barcodes:', c)
    print('total:', a + b + c)

    final_barcodes = barcodes_step4.copy()
    df['FINAL_BARCODE'] = df.apply(
        get_final_barcode, axis=1, final_barcodes=final_barcodes)
    df['FINAL_QNAME'] = df.apply(get_new_barcode, axis=1)
    
    print 'writing new samfile...'
    tempsamfile = "{}_temp.sam".format(prefix)
    finalsamfile = "{}_final.sam".format(prefix)
    finalbamfile = "{}_final.bam".format(prefix)
    finalbaifile = "{}_final.bai".format(prefix)
    with open(tempsamfile, 'w') as f:
        for line in open(samfile):
            items = line.split('\t')
            original_id = items[0]
            if not isinstance(df.loc[original_id], pd.Series):
                continue

            new_id = df.loc[original_id].FINAL_QNAME
            line = line.replace(original_id, new_id)
            f.write(line)

    os.system(
        'samtools view -H {bamfile} > {finalsamfile} && \
        cat {tempsamfile} >> {finalsamfile}'.format(**locals())
    )
    os.system('samtools view -Sb {finalsamfile} > {finalbamfile}'.format(**locals()))
    os.system('samtools index {finalbamfile} {finalbaifile}'.format(**locals()))
    print('Done--->', finalbamfile)
    return finalbamfile

if __name__ == "__main__":
    bamfile = sys.argv[1]
    _run_cluster_barcode(bamfile)
    print bamfile
