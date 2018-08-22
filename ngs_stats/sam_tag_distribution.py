#! /bin/env python

import pysam
import sys
import numpy as np
import matplotlib.pyplot as plt
plt.switch_backend('agg')

#sys.argv = [0, "test.bam"]

samfile = pysam.AlignmentFile(sys.argv[1], "rb")

notag = 0
value = list()

for read in samfile:
	try:
		# Get NM value
		value.append(read.get_tag("NM")/float(read.query_alignment_length))
	except:
		# Count reads with no NM tag
		notag+=1

#print notag
#print value

data = np.array(value)
plt.hist(data, weights=np.zeros_like(data) + 1. / data.size ) 
plt.text(0.4,0.5, "N reads = %d" % (notag+len(value)))
plt.text(0.4,0.45, "No tag = %d" % notag)
plt.text(0.4,0.4, "min = %d | max = %d | median = %d" % (np.min(value), np.max(value), np.median(value)))
plt.text(0.4,0.35, "Q.25 = %1.f | Q.75 = %.1f | Q.90 = %.1f" % (np.percentile(value,25), np.percentile(value,75), np.percentile(value,90) ) )
plt.ylim([0,1])
plt.xlim([0,1])

plt.savefig("NM-tag_distribution_%s.pdf" % sys.argv[1])
