#!/usr/bin/env python3

#Example Input: ./automate_weight_extraction.py -s /Volumes/T7/prive_sf_211_PLR_9ancestries/ -c /Users/adrianharris/Desktop/Gini-PGS/scripts/ -w /Users/adrianharris/Desktop/

import subprocess as sp
import argparse as ap
import os

parser = ap.ArgumentParser()
parser = ap.ArgumentParser(description="automation for the creation of SNP.txt files for each trait")
parser.add_argument("-c", "--code", help="directory containing the weight_extraction.py script", required=True)
parser.add_argument("-s", "--summary", help="directory containing PGS summary files", required=True)
parser.add_argument("-w", "--working", help="working directory where SNP.txt files will be stored", required=True)
args = parser.parse_args()

directory = args.summary
f = []
for filename in os.listdir(directory):
    file = os.path.join(directory, filename)
    f.append(file)

os.chdir(args.code)

for file in f:
	print(f"Running file {file}")
	autorun = sp.call(["./weight_extraction.py", "-s", file, "-w", args.working])
