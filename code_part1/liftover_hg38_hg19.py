#!/usr/bin/env python3

# Example input: ./liftover_hg38_hg19.py -m <genetic_map_path>

# 4 - liftover.py

# Converts the genetic map in `../input_data/aau1043_datas3` from GRCh38 to GRCh37 to match existing data
# Must be unzipped from the supplemental files (aau1043_datas3.gz) of:
# https://doi.org/10.1126/science.aau1043
# Script outputs to the same directory as the input file

import os
import argparse as ap
import sys 
import liftover # alternatively, the 'pyliftover' package also works, although you must switch which
				# "lo = " line is commented out below	

# Parses arguments to script
parser = ap.ArgumentParser()
parser = ap.ArgumentParser(description='Liftover Script')
parser.add_argument("-m", "--loc_map", help="input path to genetic recombination map", required=True, default="../other_data/aau1043_datas3")
args = parser.parse_args()

lo = liftover.get_lifter('hg38', 'hg19')	# Uncomment this line of code if using 'liftover' package (default)
#lo = pyliftover.LiftOver('hg38', 'hg19')	# Uncomment this line of code instead if using 'pyliftover' package

# helper function for later code
def convert(x):
	return str(x)

# reads input genetic recombination map file
loc_map = args.loc_map
inputFile = open(loc_map, 'r')
contents = inputFile.readlines()

# sets system output to be printed to output file
loc_map_out = args.loc_map + "_hg19"
outputFile = open(loc_map_out, 'w')
original_stdout = sys.stdout
sys.stdout = outputFile

# loops through lines of genetic map and updates chrom_positions accordingly
for x in range(0, len(contents)):
	line = contents[x]

	# prints header of table and identifies indices of columns
	if 'Chr	Begin	End	cMperMb	cM' in line:
		line = line.strip().split()

		chrom_index = line.index('Chr')
		begin_index = line.index('Begin')
		end_index = line.index('End')

		line = '\t'.join(line)
		print(line)
	# prints converted coordinates of map bins
	elif "#" not in line:
		line = line.strip().split()

		# extracts chromosome, beginning bp, and ending bp of line
		chrom = str(line[chrom_index])
		begin = float(line[begin_index])
		end = float(line[end_index])

		# gets converted begin and end bp coordinates
		new_begin_tuple = lo.convert_coordinate(chrom, begin)
		new_end_tuple = lo.convert_coordinate(chrom, end)

		if new_begin_tuple and new_end_tuple: #if the list contains elements (i.e. has hg19 coordinates)
			# extracts begin and end coordinates from tuple
			new_begin = new_begin_tuple[0][1]
			new_end = new_end_tuple[0][1]

			# replaces old coordinates with new coordinates
			line[begin_index] = int(new_begin)
			line[end_index] = int(new_end)

			# prints updated line
			lineList = list(map(convert, line))
			line = '\t'.join(lineList)
			print(line) 
		else:
			pass
# closes files
inputFile.close()
outputFile.close()



		



