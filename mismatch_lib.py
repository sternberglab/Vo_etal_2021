from Bio import SeqIO
from Bio.Seq import Seq
import multiprocessing
import os
import random
import csv
from pathlib import Path

def get_1bp_mms(target):
	seq = Seq(target['Spacer'])
	spacers = []
	for bp in range(len(seq)):
		if bp == 0:
			spacer = seq.complement()[0] + seq[1:]
		else:
			spacer = seq[0:bp] + seq.complement()[bp] + seq[bp+1:]
		spacers.append(spacer)
	return spacers

def get_2bp_mms(target):
	seq = Seq(target['Spacer'])
	spacers = []
	for bp in range(len(seq)-1):
		if bp == 0:
			spacer = seq.complement()[0:2] + seq[2:]
		else:
			spacer = seq[0:bp] + seq.complement()[bp:bp+2] + seq[bp+2:]
		spacers.append(spacer)
	return spacers

def get_4bp_mms(target):
	seq = Seq(target['Spacer'])
	spacers = []
	for bp in range(len(seq)-3):
		if bp == 0:
			spacer = seq.complement()[0:4] + seq[4:]
		else:
			spacer = seq[0:bp] + seq.complement()[bp:bp+4] + seq[bp+4:]
		spacers.append(spacer)
	return spacers

def hamming_dist(s1, s2):
    assert len(s1) == len(s2)
    return sum(ch1 != ch2 for ch1, ch2 in zip(s1, s2))

def get_barcode(barcodes):
	# returns a 12bp barcode at least 5 edit distance from all existing once
	valid = False
	while not valid:
		valid = True
		potential = ''.join(random.choices('CGAT', k=12))
		for barcode in barcodes:
			if hamming_dist(barcode, potential) < 5:
				valid = False
				continue
	barcodes.append(potential)
	return potential

def main():
	# get targets from input csv
	targets = []
	with open('mm_library_input.csv', 'r', encoding='utf-8-sig') as in_file:
		reader = csv.DictReader(in_file)
		targets = [r for r in reader]

	barcodes = []
	spacer_id = 1
	# open output file for writing
	with open('mm_library_output.csv', 'w', newline='') as out_file:
		writer = csv.writer(out_file)
		writer.writerow(['ID', 'TargetID', 'TargetDescription', 'SpacerDescription', 'Spacer', 'Barcode (12bp)'])
		# generate mm sets for each target
		for target in targets:
			writer.writerow([spacer_id, target['TargetID'], target['Description'], 'Wild-type', target['Spacer'], get_barcode(barcodes)])
			spacer_id += 1
			for (i, spacer) in enumerate(get_1bp_mms(target)):
				writer.writerow([spacer_id, target['TargetID'], target['Description'], f"1mm({i+1})", spacer,  get_barcode(barcodes)])
				spacer_id += 1
			for (i, spacer) in enumerate(get_2bp_mms(target)):
				writer.writerow([spacer_id, target['TargetID'], target['Description'], f"2mm({i+1}-{i+2})", spacer,  get_barcode(barcodes)])
				spacer_id += 1
			for (i, spacer) in enumerate(get_4bp_mms(target)):
				writer.writerow([spacer_id, target['TargetID'], target['Description'], f"4mm({i+1}-{i+4})", spacer,  get_barcode(barcodes)])
				spacer_id += 1

main()