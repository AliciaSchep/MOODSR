#!/usr/bin/python

import MOODS.scan
import MOODS.tools

import sys
import time
import os
import random

def load_matrix(filename):
    with open(filename, "r") as file_handle:
        mat = [[float(i) for i in line.split()] for line in file_handle]
        return mat

if len(sys.argv) < 2:
	print "usage: python ex-scanner.py [matrix directory]"
	sys.exit(1)

# read the matrix files
# should be in JASPAR .pfm format
matrix_directory = sys.argv[1]
matrix_names = [filename for filename in os.listdir(matrix_directory) if filename[-4:] == '.pfm']

matrices = [load_matrix(matrix_directory + filename) for filename in matrix_names]
bg = MOODS.tools.flat_bg(4)
# log-odds transformation for matrices
matrices = [MOODS.tools.log_odds(m, bg, 1) for m in matrices]
# thresholds computed from p-value
thresholds = [MOODS.tools.threshold_from_p(m, bg, 0.0001) for m in matrices]


# instead of just calling the scan function, we'll build a scanner object
# to get a persistent scanner that can be used to scan multiple sequences
# without having to repeat the preprocessing each time
scanner = MOODS.scan.Scanner(7) # parameter is the window size

scanner.set_motifs(matrices, bg, thresholds)

# we'll generate 100 random sequences
for i in xrange(100):
	seq = ''.join([random.choice(['A', 'C', 'G', 'T']) for j in range(10000)])
	results = scanner.scan(seq)

	# you can specify a limit for the number of hits:
	#
	# results = scanner.scan(seq, 10)
	#
	# this is mostly intended to prevent things from slowing down too much
	# when the threshold is too loose for some reason

	print "Sequence:", i, "Hits:", sum([len(r) for r in results])