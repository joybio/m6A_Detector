#!/root/miniconda3/bin/python

"""extract.cds.py """

__date__ = "2019-10-2"
__author__ = "Junbo Yang"
__email__ = "yang_junbo_hi@126.com"
__license__ = "PKU.jia.group"

#imports
import re
import os
import optparse
#import bisect
#sort  package
import linecache
from optparse import OptionParser
import collections
from collections import defaultdict

parser = OptionParser('Usage: %prog ')
parser.add_option('-i','--in',
                dest='input',
                help='input file')
parser.add_option('-o','--out',
                dest='out',
                help='count and percentage file')
(options,args) = parser.parse_args()

data = open(options.input,'r')
out= open(options.out,'w')

for i in data:
	line = i
	i = i.strip().split("\t")
	chrom = i[0]
	start = int(i[1])
	stop = int(i[2])
	cov = i[3]
	while stop - start >= 1:
		out.write(chrom + "\t" + str(start) + "\t" + str(start+1) + "\t" + cov + "\n")
		start += 1
data.close()
out.close()

