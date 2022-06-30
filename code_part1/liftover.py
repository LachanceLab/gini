#!/usr/bin/env python3

#Example input: ./liftover.py -s /Volumes/T7/prive_sf_211_PLR_9ancestries/153-1KG_PLR.txt -o /Users/adrianharris/Desktop/
import argparse as ap
import sys 
import pyliftover
import os

parser = ap.ArgumentParser()
parser = ap.ArgumentParser(description='Liftover Script')
parser.add_argument("-s", "--summary", help="input pgs summary file here", required=True)
parser.add_argument("-o", "--out_dir", help="output directory", required=True)
args = parser.parse_args()

def convert(x):
	return str(x)

chrom_index = ''
chrom_pos_index = ''
lo = pyliftover.LiftOver('hg19', 'hg38')
file = args.summary
inputFile = open(file, 'r')
filename = file.split('/')
filename = filename[-1]
filename = filename[:-4]
filename += '_hg38.txt'
filename = args.out_dir + filename
outputFile = open(filename, 'w')
original_stdout = sys.stdout
sys.stdout = outputFile
contents = inputFile.readlines()
for x in range(0, len(contents)):
	line = contents[x]
	if 'chr_position' in line:
		line = line.strip().split()
		chrom_index = line.index('chrom')
		chrom_pos_index = line.index('chr_position')
		line = '\t'.join(line)
		print(line)
	else:
		line = line.strip().split()
		chrom = str(line[chrom_index])
		if 'chr' not in chrom:
			chrom = 'chr' + chrom
		chrom_position = float(line[chrom_pos_index])
		new_position_tuple = lo.convert_coordinate(chrom, chrom_position)
		if new_position_tuple: #if the list contains elements
			new_position = new_position_tuple[0][1]
			line[chrom_pos_index] = int(new_position)
			lineList = list(map(convert, line))
			line = '\t'.join(lineList)
			print(line)
			#need a way to write to this specific position of 
		else:
			pass
#sys.stdout = original_stdout
inputFile.close()
outputFile.close()



		



