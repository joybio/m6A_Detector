#!/root/miniconda3/bin/python
__date__ = "2022-7-4"
__author__ = "Junbo Yang"
__email__ = "yang_junbo_hi@126.com"
__license__ = "PKU.jia.group"

#imports
import re
import os
import optparse
from optparse import OptionParser
import math
import numpy as np
from numpy import *
#from scipy.stats import ttest_ind
import time
from time import strftime

import sys
from sys import argv
parser = OptionParser('Usage: %prog -f reads_start_from_genomecov_XS+ -r reads_start_from_genomecov_XS- -n 3 -l reads_len -o output_file')

parser.add_option('-i','--in',
	dest='input',
	help="input file.")

parser.add_option('-f','--for',
	dest='forward',
	help="depth of forward")
parser.add_option('-r','--rev',
	dest='reverse',
	help='depth of reverse')
parser.add_option('-n','--num',
		dest='num',
		help='number: which column should choose.')
parser.add_option('-a','--arrest',
		dest='arrest',
		help='arrest rate for filter.')
parser.add_option('-o','--out',
	dest='out',
	help='Output annotation file.')
parser.add_option('-l','--filter',
	dest='filter',
	help='filter results; default [Reads end>=3 and pvalue <=0.05 and arrest rate>=0.3] .')

(options,args) = parser.parse_args()

if len(sys.argv) == 1:
	parser.print_help()
	sys.exit(1)

print("INFO {} Thread 1: Start analysis......".format(time.strftime("%Y-%m-%d %H:%M:%S",time.localtime(time.time()))))

filtered_stop = open(options.input,"r")
number = int(options.num)
arrest = float(options.arrest)
forward = open(options.forward,"r")
reverse = open(options.reverse,"r")
filtered = open(options.filter,"w")
forward_end_dict = {}
reverse_end_dict = {}
for i in forward:
	i = i.strip().split("\t")
	chrom = i[0]
	start = int(i[1])-1
	stop = int(i[1])
	key = chrom+"\t"+str(start)+"\t"+str(stop)
	cov = int(i[number])
#   print(cov)	 
	forward_end_dict[key] = cov
forward.close()
#first step: expand the result of bedtools genomecov and construct dict
for i in reverse:
	i = i.strip().split("\t")
	chrom = i[0]
	start = int(i[1])-1
	stop = int(i[1])
	key = chrom+"\t"+str(start)+"\t"+str(stop)
	cov =int(i[number])
	reverse_end_dict[key] = cov
reverse.close()
out = open(options.out,"w")
for i in filtered_stop:
	if i.startswith("#Chrom"):
		i = i.strip()
		out.write(i+"\tarrest rate\tcurrent coverage\tflank coverage\n")
	elif i.startswith("#"):
		out.write(i)
	else:
		line = i.strip()
		i = i.strip().split("\t")
		chrom = i[0]
		start = int(i[1])
		strand = i[4]
		pvalue = float(i[-1])
		Reads_end = int(i[3])
		key = chrom+"\t"+str(start)+"\t"+str(start+1)
		left = chrom+"\t"+str(start-1)+"\t"+str(start)
		right = chrom+"\t"+str(start+1)+"\t"+str(start+2)
		if strand == "+":
			if left in forward_end_dict.keys() and key in forward_end_dict.keys():
        			arrest_rate = round((1-forward_end_dict[key]/(forward_end_dict[left]+0.01)),3)
	        		#arrest_rate = round(1-mean([x/(y+0.01) for x,y in zip(forward_end_dict[key],forward_end_dict[left])]),3)
		        	out.write(line + "\t" + str(arrest_rate) + "\t" + str(forward_end_dict[key]) + "\t" + str(forward_end_dict[left]) + "\n")
			#if all([int(x) >=3 for x in Reads_end]) and pvalue <= 0.05 and arrest_rate >= arrest:
        			if Reads_end >=3 and pvalue <= 0.05 and arrest_rate >= arrest:
	        			#filtered.write(line + "\t" + str(arrest_rate) + "\t" + "|".join([str(x) for x in forward_end_dict[key]]) + "\t" + "|".join([str(x) for x in forward_end_dict[left]])+ "\n")
		        		filtered.write(line + "\t" + str(arrest_rate) + "\t" + str(forward_end_dict[key]) + "\t" + str(forward_end_dict[left]) + "\n")
		else:
			if right in reverse_end_dict.keys() and key in reverse_end_dict.keys():
        			arrest_rate = round((1-reverse_end_dict[key]/(reverse_end_dict[right]+0.01)),3)
	        		#arrest_rate = round(1-mean([x/(y+0.01) for x,y in zip(reverse_end_dict[key],reverse_end_dict[right])]),3)
		        	#out.write(line + "\t" + str(arrest_rate) + "\t" + "|".join([str(x) for x in reverse_end_dict[key]]) + "\t" + "|".join([str(x) for x in reverse_end_dict[right]]) + "\n")
        			out.write(line + "\t" + str(arrest_rate) + "\t" + str(reverse_end_dict[key]) + "\t" + str(reverse_end_dict[right]) + "\n")
			#if all([int(x) >=3 for x in Reads_end]) and pvalue <= 0.05 and arrest_rate >= arrest:
        			if Reads_end >=3 and pvalue <= 0.05 and arrest_rate >= arrest:
	        			#filtered.write(line + "\t" + str(arrest_rate) + "\t" + "||".join([str(x) for x in reverse_end_dict[key]]) + "\t" + "|".join([str(x) for x in reverse_end_dict[right]])+ "\n")
		        		filtered.write(line + "\t" + str(arrest_rate) + "\t" + str(reverse_end_dict[key]) + "\t" + str(reverse_end_dict[right]) + "\n")
out.close()
print("INFO {} DONE.".format(time.strftime("%Y-%m-%d %H:%M:%S",time.localtime(time.time()))))

