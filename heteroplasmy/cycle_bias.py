#! /bin/python
import os
import sys
import pysam
import numpy as np

def cycle_bias(pos,minor,bam,win=15):
    tot_minor_reads = 0
    minor_in_edge = 0
    samfile = pysam.AlignmentFile(bam,"rb")
    
    for pileupcolumn in samfile.pileup('MT', int(pos)-1, int(pos), stepper='all',max_depth=10000000, mask=False,truncate=True):
        if pileupcolumn.pos==int(pos)-1:
            for pileupread in pileupcolumn.pileups:
                position_in_segment = pileupread.query_position
                aligned_segment = pileupread.alignment
                if position_in_segment is not None:
                
                    segment_minor = aligned_segment.query_sequence[position_in_segment]
                    segment_qual = aligned_segment.query_qualities[position_in_segment]
                    if segment_minor == minor and segment_qual >=30:

                        tot_minor_reads+=1

                        rlen = aligned_segment.query_length
                        
                        edges = range(0,win)+range(rlen-win,rlen+1)
            
                        if position_in_segment in edges:
                            minor_in_edge+=1
    try:
        ratio = minor_in_edge/float(tot_minor_reads)
        
    except:
        ratio = np.nan
    
    
    return ratio


pos,minor,bam = sys.argv[1:4]
print cycle_bias(pos,minor,bam)

