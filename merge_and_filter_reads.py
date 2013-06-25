#!/usr/bin/env python

import sys
import os
import subprocess

# This script takes one or more input files containing mappings in SAM,
# BAM, or ELAND format. The files are converted to ELAND format, filtered
# to include only uniquely mapped reads, and written to a single ELAND
# output file.
#
# TODO: This script could be eliminated if the pipeline app took Mappings
# tables as input. Internally these could be dumped to ELAND format for
# compatibility with the existing pipeline scripts.

from eland import ElandExtendedFile, ElandMultiFile, ElandFile, BwaSamFile, BowtieSamFile, ElandSamFile, IlluminaSamFile

def convert_bwasam(eland_output, sam_input, mismatches):
	'''Converts SAM file from BWA to Eland, filters out all but unique reads with no more than the specified number of mismatches'''
	input = BwaSamFile(sam_input, 'r')
	total_passed = 0
	for i, line in enumerate(input):
		if line.num_best_hits != 1:
			continue
		elif line.num_mismatches > mismatches:
			continue
		else:
			total_passed += 1
			eland_output.write(line.convert_to_eland())
	print "bwasam: total lines", i, "total passed", total_passed
	input.close()

def convert_bowtiesam(eland_output, sam_input, mismatches):
	'''Converts SAM file from Bowtie to Eland, filters out all but unique reads with no more than the specified number of mismatches.'''
	input = BowtieSamFile(sam_input, 'r')
	total_passed = 0
	for i, line in enumerate(input):
		if line.num_hits != 1:
			continue
		elif line.edit_distance > mismatches:
			continue
		else:
			total_passed += 1
			eland_output.write(line.convert_to_eland())
	print "bowtiesam: total lines", i, "total passed", total_passed
	input.close()

def convert_illuminasam(eland_output, sam_input, mismatches):
	'''Converts SAM file from Bowtie to Eland, filters out reads with no more than the specified number of mismatches. NOTE: Due to limitations of Illumina format, does not filter out unique reads. Assumes prefiltered.'''
	input = IlluminaSamFile(sam_input, 'r')
	total_passed = 0
	for i, line in enumerate(input):
		if line.mismatching_positions > mismatches:
			continue
		else:
			total_passed += 1
			eland_output.write(line.convert_to_eland())
	print "illuminasam: total lines", i, "total passed", total_passed
	input.close()
	eland_output.close()

def convert_elandsam(eland_output, sam_input, mismatches):
	input = ElandSamFile(sam_input, 'r')
	total_passed = 0
	for i, line in enumerate(input):
		if line.num_hits != 1:
			continue
		elif line.edit_distance > mismatches:
			continue
		else:
			total_passed += 1
			eland_output.write(line.convert_to_eland())
	print "elandsam: total lines", i, "total passed", total_passed
	input.close()

def convert_sam(eland_output, sam_input, mismatches):
	'''Determines SAM format and converts to eland'''
	input = open(sam_input, 'r')
	line = input.readline()
	# Look for proper headers first
	while line.startswith('@'):
		if line.startswith('@PG'):
			fields = line.split('\t')
			program_name = None
			for f in fields:
				if f.startswith('PN:'):
					program_name = f[3:]
					continue
				if f.startswith('ID:'):
					program_name = f[3:]
					continue
			if program_name == 'bwa':
				input.close()
				convert_bwasam(eland_output, sam_input, mismatches)
				return
			elif program_name == 'bowtie':
				input.close()
				convert_bowtiesam(eland_output, sam_input, mismatches)
			elif program_name == 'illumina_export2sam.pl':
				input.close()
				convert_illuminasam(eland_output, sam_input, mismatches)
				return
		line = input.readline()

	# If no proper headers, try to guess appropriate SAM format
	fields = line.split('\t')
	for f in fields[11:]:
		if f.startswith('X0'):
			input.close()
			convert_bwasam(eland_output, sam_input, mismatches)
			return
		if f.startswith('NM'):
			input.close()
			convert_elandsam(eland_output, sam_input, mismatches)
			return
	input.close()
	if len(fields) == 11:  # No Custom flags
		raise Exception("Cannot convert regular SAM file")
	else:
		convert_bowtiesam(eland_output, sam_input, mismatches)

def convert_bam(eland_output, bam_input, mismatches):
        # Figure out how to put this in a uniquely-named temp file
	sam_input = os.path.basename(bam_input) + '.sam'
	bam2sam_cmd = 'samtools view -h %s > %s' % (bam_input, sam_input)
	subprocess.call(bam2sam_cmd, shell=True)
	convert_sam(eland_output, sam_input, mismatches)
	os.remove(sam_input)

def merge_unique_eland(output, mapped_reads_files, mismatches=2):
	eland_out = ElandFile(output, 'w')
	for i in mapped_reads_files:
		if not os.path.exists(i):
			raise Exception("File %s does not exist" % i)
		if i.endswith('.bam'):
                        print "Parsing file as BAM..."
			convert_bam(eland_out, i, mismatches)
			continue
		if i.endswith('.sam'):
                        print "Parsing file as SAM..."
			convert_sam(eland_out, i, mismatches)
			continue
		if 'multi' in i:
                        print "Parsing file as ELAND multi..."
			eland_in = ElandMultiFile(i, 'r')
		elif 'extended' in i:
                        print "Parsing file as ELAND extended..."
			eland_in = ElandExtendedFile(i, 'r')
		else:
                        print "Parsing file as plain ELAND..."
			eland_in = ElandFile(i, 'r')
		total_passed = 0
		for i, line in enumerate(eland_in):
			best_hits = line.best_matches()
			if len(best_hits) == 0:
				continue
			elif len(best_hits) > 1:
				continue  # Only merge unique hits
			elif best_hits[0].number_of_mismatches() > mismatches:
				continue
			else:
				total_passed += 1
				eland_out.write(line.convert_to_eland())
		print "unique eland: total lines", i, "total passed", total_passed
		eland_in.close()
	eland_out.close()

if __name__ == '__main__':
	if len(sys.argv) < 3:
		print "Usage: merge_and_filter_reads.py output_file reads_file [reads_file ...]"
		raise SystemExit(1)
	merge_unique_eland(sys.argv[1], sys.argv[2:])
