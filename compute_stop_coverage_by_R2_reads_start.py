#!/root/miniconda3/bin/python

"""
#1) Create a summary that shows potential m6A sites of reads start or end\n
#2) number of Reads start or end must be [score]-fold greater than the average number of 15 nt flanking region \n
#3) if there is more than 5 continuouse positon don't have reads end in either flanking orientation. the continuouse positon and the left position in the flanking region will not take into account\n
#4) if the number of reads start/end is zero, but the position is not one of the five continuouse positon.The reads start or end will be assigned to one.\n
#5) the number of reads end of position i (m6A) must be no less than (args) in the IP sample\n
"""

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
from scipy import stats
import numpy as np
from numpy import *
#from scipy.stats import ttest_ind
import time
from time import strftime
import sys
from sys import argv

parser = OptionParser('Usage: %prog -f reads_start_from_genomecov_XS+ -r reads_start_from_genomecov_XS- -n 3 -l reads_len -o output_file')
parser.add_option('-f','--for',
                dest='forward',
                help='Reads end of forward. This information can get from bedtools genomecov -scale 1 [-5,-3].which depend on the stratagy of library construction (e.g. 3-ligation: -5; dUTP: -3). Take care that when you use hisat2 --rna-strandness, After mapping you can get "XS:A:+" from sam or bam files. XS attribute tag: + means a read belongs to a transcript on + strand of genome. - means a read belongs to a transcript on - strand of genome. which is very different from the strand information in bedtools genomecov (forward/reverse), specifically, + strand in bedtools is equal to the - strand in XS attribute tag. Therefore, we suggest you split sam file by XS attribute tag and use the split file as the input of genomecov instead of using strand parameter in genomecov.')

parser.add_option('-r','--rev',
                dest='reverse',
                help='Reads end of reverse. Same to forward.')

parser.add_option('-o','--out',
                dest='out',
                help='Output annotation file. output the next position of truncation site.')

parser.add_option('-n','--num',
                dest='num',
				default = 3,
                help='Min number of R2 start (number of truncation). default [0]. only positions with reads end number greater than this number will take into account.')

parser.add_option('-s','--score',
                dest='score',
				default = 2,
                help='Min score of (R2 start number)/(average flanking number), which means that the reads end number of the current position is greater than the average of flanking. default [2]. only positions with scores greater than this number will take into account.')

parser.add_option('-l','--flank',
                dest='flank',
                default = 15,
                help='Flanking region size (background for the arrest position. the first choice is reads length). default [15]')


(options,args) = parser.parse_args()
if len(sys.argv) == 1:
	parser.print_help()
	sys.exit(1)

print("INFO {} Thread 1: Start analysis......".format(time.strftime("%Y-%m-%d %H:%M:%S",time.localtime(time.time()))))

forward = open(options.forward,"r")
reverse = open(options.reverse,"r")
number = int(options.num)
print("INFO: Min number of Truncation is: {}".format(number))
# min number of R2 start. if not detected, then use the default [3]
score = float(options.score)
print("INFO: Min score of Truncation is: {}".format(score))
# min score of (number of R2 start)/(average number of flanking region). if not detected, then use the default [2]
flank = int(options.flank)
out = open(options.out,"w")

out.write("#1) Create a summary that shows potential m6A sites of reads end.\n")
out.write("#2) number of Reads start or end must be [score]-fold greater than the average number of [flank] flanking region \n")
out.write("#3) if there is more than 5 continuouse positon have zero reads of start/end in either flanking orientation. the continuouse positon and the left position in the flanking region will not take into account.\n")
out.write("#4) if the number of reads start/end is zero, but the position is not located in the continuouse zero positon.The reads start or end will be assigned to one.\n")
out.write("#5) the number of reads end of position i (m6A) must be no less than (args) in the samples\n")
out.write("#Chrom\tStart\tStop\tR2_start\tStrand\tScore(left|right)\tPvalue\n")
forward_end_dict = {}
reverse_end_dict = {}
print("INFO {}\n Expand the result of bedtools genomecov".format(time.strftime("%Y-%m-%d %H:%M:%S",time.localtime(time.time()))))
for i in forward:
	i = i.strip().split("\t")
	chrom = i[0]
	start = int(i[1])
	stop = int(i[2])
	cov = mean([int(x) for x in i[3::]])
	while stop - start >= 1:
		key = chrom+"\t"+str(start-1)+"\t"+str(start)
		forward_end_dict[key] = cov
		start += 1
forward.close()
#first step: expand the result of bedtools genomecov and construct dict
for i in reverse:
	i = i.strip().split("\t")
	chrom = i[0]
	start = int(i[1])
	stop = int(i[2])
	cov = mean([int(x) for x in i[3::]])
	while stop - start >= 1:
		key = chrom+"\t"+str(start-1)+"\t"+str(start)
		forward_end_dict[key] = cov
		start += 1	
reverse.close()
#same to step1

def flank_continue(chrom,start,stop,tag):
	if tag == "forward":
		end_dict = forward_end_dict
	else:
		end_dict = reverse_end_dict
	l = 1
	left_continuouse = 0
	r = 1
	right_continuouse = 0
	coverage_left = []
	coverage_right = []
	# if there is five continuouse position don't have reads end, then pass
	while l <= flank and left_continuouse < 5: 
		#left and current position
		left = chrom + "\t" + str(start-l) + "\t" + str(start-l+1)
		current = chrom + "\t" + str(start-l+1) + "\t" + str(start-l+2)
		#if left position in the reads end dict, append this position to the list.
		if left in end_dict.keys():
			coverage_left.append(end_dict[left])
		else:
			#if left position not in the reads end dict but the current position in the reads end dict, append left position to the list too.
			if current in end_dict.keys():
				left_continuouse = 0
				coverage_left.append(1)
			# if left and current position neither in the the reads end dict, pass
			else:
				left_2 = chrom + "\t" + str(start-l-1) + "\t" + str(start-l)
				left_3 = chrom + "\t" + str(start-l-2) + "\t" + str(start-l-1)
				if left_2 or left_3 in end_dict.keys():
					coverage_left.append(1)
				left_continuouse += 1
		l += 1
		#return(coverage_left)
	#echo the average of the left coverage; if there are five continuouse position dont have reads end. then stop flanking.
	while r <= flank and right_continuouse < 5:
		right = chrom + "\t" + str(stop+r-1) + "\t" + str(stop+r)
		current = chrom + "\t" + str(start-r+1) + "\t" + str(start-r+2)
		if right in end_dict.keys():
			coverage_right.append(end_dict[right])
		else:
			if current in end_dict.keys():
				right_continuouse = 0
				coverage_right.append(1)
			else:
				right_2 = chrom + "\t" + str(start+r-2) + "\t" + str(start+r-1)
				right_3 = chrom + "\t" + str(start+r-3) + "\t" + str(start+r-2)
				if right_2 or right_3 in end_dict.keys():
					coverage_right.append(1)
				right_continuouse += 1
		r += 1
		return([coverage_left,coverage_right])


forward = open(options.forward,"r")
reverse = open(options.reverse,"r")
print("INFO {} Thread 1: Working on forward......".format(time.strftime("%Y-%m-%d %H:%M:%S",time.localtime(time.time()))))
for i in forward:
	tag = "forward"
	i = i.strip().split("\t")
	chrom = i[0]
	start = int(i[1])
	stop = int(i[2])
	coverage = [int(x) for x in i[3::]]
	end = [str(x) for x in i[3::]]
	left_coverage = flank_continue(chrom,start,stop,tag)[0]
	#print(left_coverage)
	if len(left_coverage) > 1:
		left_coverage_mean = mean(left_coverage)
		left_coverage_sum = sum(left_coverage)
	else:
		left_coverage_mean = 1
		left_coverage_sum = 0
	right_coverage = flank_continue(chrom,start,stop,tag)[1]
	if len(right_coverage) > 1:
		right_coverage_mean = mean(right_coverage)
		right_coverage_sum = sum(right_coverage)
	else:
		right_coverage_mean = 1
		right_coverage_sum = 0
	score_left = mean(coverage)/left_coverage_mean
	score_right = mean(coverage)/right_coverage_mean
	total_reads_end = sum(coverage) + left_coverage_sum + right_coverage_sum
	expect = total_reads_end/(len(left_coverage)+1+len(right_coverage))
	#possion distribution pmf: coverage is True; pdf: less than coverage is True
	pvalue = stats.poisson.pmf(mean(coverage),expect)
	#Beta distribution
	if score_left >= score and score_right >= score and all([c >= number for c in coverage]):
		out.write(chrom + "\t" + str(stop) + "\t" + str(stop+1)  + "\t" + '|'.join(end)+ "\t.\t" + "+\t" + str(round(score_left,2))+"|"+str(round(score_right,2)) + "\t" + str(round(pvalue,4)) +"\n")
forward.close()

print("INFO {} Thread 1: Working on reverse......".format(time.strftime("%Y-%m-%d %H:%M:%S",time.localtime(time.time()))))
for i in reverse:
	tag = "reverse"
	i = i.strip().split("\t")
	chrom = i[0]
	start = int(i[1])
	stop = int(i[2])
	coverage = [int(x) for x in i[3::]]
	end = [str(x) for x in i[3::]]
	left_coverage = flank_continue(chrom,start,stop,tag)[0]
	if len(left_coverage) > 1:
		left_coverage_mean = mean(left_coverage)
		left_coverage_sum = sum(left_coverage)
	else:
		left_coverage_mean = 1
		left_coverage_sum = 0
	right_coverage = flank_continue(chrom,start,stop,tag)[1]
	if len(right_coverage) > 1:
		right_coverage_mean = mean(right_coverage)
		right_coverage_sum = sum(right_coverage)
	else:
		right_coverage_mean = 1
		right_coverage_sum = 0
	score_left = mean(coverage)/left_coverage_mean
	score_right = mean(coverage)/right_coverage_mean
	total_reads_end = sum(coverage) + left_coverage_sum + right_coverage_sum
	#possion distribution pmf: coverage is True; pdf: less than coverage is True
	expect = total_reads_end/(len(left_coverage)+1+len(right_coverage))
        #possion distribution pmf: coverage is True; pdf: less than coverage is True	
	pvalue = stats.poisson.pmf(mean(coverage),expect)
	if score_left >= score and score_right >= score and all([c >= number for c in coverage]):
		out.write(chrom + "\t" + str(start-1) + "\t" + str(start)  + "\t" + '|'.join(end) + "\t.\t" + "-" + "\t" + str(round(score_left,2))+"|"+str(round(score_right,2)) + "\t" + str(round(pvalue,4))+ "\n")
reverse.close()
out.close()
print("INFO {} DONE.".format(time.strftime("%Y-%m-%d %H:%M:%S",time.localtime(time.time()))))




















