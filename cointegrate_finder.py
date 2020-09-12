from Bio import SeqIO, SearchIO
from Bio.SeqRecord import SeqRecord
from Bio.Blast.Applications import NcbiblastnCommandline
import gzip
from datetime import date
import shutil
import multiprocessing
import subprocess
import os
import csv
import boto3
import botocore
from pathlib import Path
from simplesam import Reader as samReader

today = date.today()
output_name = f'all-{today.year}-{today.month}-{today.day}'
# install blast locally
# curl -O ftp://ftp.ncbi.nlm.nih.gov/blast/executables/LATEST/ncbi-blast-2.10.1+-x64-linux.tar.gz && tar xzf ncbi-blast-2.10.1+-x64-linux.tar.gz && sudo cp ncbi-blast-2.10.1+/bin/* /usr/local/bin
run_local = False

def main():
	os.makedirs(os.path.join(f"./outputs"), exist_ok=True)
	os.makedirs(os.path.join(f"./bt2index"), exist_ok=True)
	os.makedirs(os.path.join(f"./tmp"), exist_ok=True)
	today = date.today()
	with open(f'outputs/{output_name}.csv', 'w', newline='') as outfile:
		writer = csv.writer(outfile)
		writer.writerow(['Read_file', 'total_tn_reads', 'cointegrates', 'genomic_insertions', 'plasmids', 'insufficient', 'unknown', 'Sample Description'])
	
	with open('input2.csv', 'r', encoding='utf-8-sig') as infile:
		reader = csv.DictReader(infile)
		for row in reader:
			reads_file = row['Reads File']
			tn_file = row['Tn File']
			plasmid_file = row['Plasmid File']
			genome_file = row['Genome File']
			sample = row['Sample Code']
			sample_desc = row['Sample Desc']
			print(f"Processing {sample}...")
			print("-----")
			process_sample(reads_file, tn_file, plasmid_file, genome_file, sample, sample_desc)
			print("-----")

	if not run_local:
		s3 = boto3.client('s3')
		for filename in os.listdir('outputs'):
			s3key = f"cointegrate_outputs/{filename}"
			s3.upload_file(f"outputs/{filename}", "sternberg-sequencing-data", s3key)
	print("done")

def download_s3(filename, isNotSample=True):
	isGzip = False
	if isNotSample:
		local_path = filename
		s3key = f'mssm_inputs/{filename}'
	else:
		s3key = f'mssm/{filename}'
		if filename.endswith('.gz'):
			isGzip = True
			local_path = "sample.fasta.gz"
		else:
			local_path = "sample.fasta"
	print(f"Downloading... {s3key}")
	s3 = boto3.client('s3')
	s3.download_file('sternberg-sequencing-data', s3key, local_path)
	if isGzip:
		with gzip.open(local_path, 'rb') as f_in:
		    with open(local_path[:-3], 'wb') as f_out:
		        shutil.copyfileobj(f_in, f_out)
		local_path = local_path[:-3]
	return local_path

def process_sample(reads_file, tn_file, plasmid_file, genome_file, sample, sample_desc):
	reads_file = download_s3(reads_file, False)
	# reads_file = 'sample.fasta'
	tn_file = download_s3(tn_file)
	plasmid_file = download_s3(plasmid_file)
	genome_file = download_s3(genome_file)
	plasmid = SeqIO.read(plasmid_file, 'fasta')
	tn = SeqIO.read(tn_file, 'fasta')
	genome = SeqIO.read(genome_file, 'fasta')
	tn_length = len(tn)

	tn_start = plasmid.seq.find(tn.seq)
	plasmid_l = plasmid.seq[tn_start-20:tn_start]
	plasmid_r = plasmid.seq[tn_start+len(tn.seq):tn_start+len(tn.seq)+20]

	basename = "tmp/" + sample
	blast_filename =  f'cointegrate_outputs/{sample}_blastresults.xml'
	do_blast(tn_file, reads_file, blast_filename)
	
	# We are querying the transposon against all the reads, so only one result with many hits is output
	res = SearchIO.read(blast_filename, 'blast-xml')
	tn_read_ids = [hit.id for hit in res.hits]

	# filter the reads to just those with the transposon, write to file in case we want them later
	all_reads = SeqIO.parse(reads_file, 'fasta')
	tn_reads = [r for r in all_reads if r.id in tn_read_ids]
	SeqIO.write(tn_reads, f'{basename}_tnreads.fasta', 'fasta')
	tn_reads = list([r for r in SeqIO.parse(f'{basename}_tnreads.fasta', 'fasta')])

	all_results = []
	# each hit represents one or more matches of the query against a read sequence
	for hit in res.hits:
		tn_read = next(r for r in tn_reads if r.id == hit.id)
		read_obj = get_read_obj(hit, tn_read, tn_length)
		if read_obj:
			all_results.append(read_obj)
	attach_alignments(all_results, basename, plasmid_file, genome_file)
	with open(f'outputs/output_{sample}.csv', 'w', newline='') as outfile:
		writer = csv.writer(outfile)
		writer.writerow(['read_id', 'type', 'hsps', 'ends', 'types', 'read_length'])
		for read in all_results:
			hsps_formatted = ''
			for pair in read['hsps']:
				if len(hsps_formatted):
					hsps_formatted += ' _ '
				hsps_formatted += '-'.join([str(i) for i in pair])
			hsps = hsps_formatted
			ends = [e['id'][-19:] for e in read['ends']]
			types = [e['type'] for e in read['ends']]
			read['end_types'] = types
			read_length = len(read['read_seqrec'].seq)
			writer.writerow([read['id'], read['type'], hsps, ends, types, read_length])
	with open(f'outputs/{output_name}.csv', 'a', newline='') as outfile:
		writer = csv.writer(outfile)
		cointegrates = [r for r in all_results if r['type'] is 'COINTEGRATE']
		# only use full cointegrate sequences for the sample fasta output
		selected_cointegrates = [r for r in cointegrates if len(r['end_types']) is 4]
		if len(selected_cointegrates) < 2 and len(cointegrates) > 9:
			selected_cointegrates = cointegrates
		if len(selected_cointegrates):
			SeqIO.write([r['read_seqrec'] for r in selected_cointegrates[:10]], f"outputs/{sample}_cointegrates.fasta", "fasta")
		genome_insertions = [r for r in all_results if r['type'] is 'GENOME']
		if len(genome_insertions):
			SeqIO.write([r['read_seqrec'] for r in genome_insertions[:10]], f"outputs/{sample}_genomic.fasta", "fasta")
		plasmids = [r for r in all_results if r['type'] is 'pl']
		if len(plasmids):
			SeqIO.write([r['read_seqrec'] for r in plasmids[:10]], f"outputs/{sample}_plasmids.fasta", "fasta")
		insufficients = [r for r in all_results if r['type'] is 'partialRead']
		if len(insufficients):
			SeqIO.write([r['read_seqrec'] for r in insufficients[:10]], f"outputs/{sample}_insufficient.fasta", "fasta")
		unknown = [r for r in all_results if r['type'] is 'unknown']
		if len(unknown):
			SeqIO.write([r['read_seqrec'] for r in unknown[:10]], f"outputs/{sample}_unknowns.fasta", "fasta")
		writer.writerow([sample, len(all_results), len(cointegrates), len(genome_insertions), len(plasmids), len(insufficients), len(unknown), sample_desc])
	print("done")

def attach_alignments(results, basename, plasmid_file, genome_file):
	all_fp_seqs = [end['seqrec'] for r in results for end in r['ends']]
	tmp_fp_fasta_name = f'{basename}_fps.fasta'
	SeqIO.write(all_fp_seqs, tmp_fp_fasta_name, 'fasta')
	do_bowtie(tmp_fp_fasta_name, plasmid_file, f"{basename}_plasmid.sam")
	do_bowtie(tmp_fp_fasta_name, genome_file, f"{basename}_genome.sam")
	with open(f"{basename}_plasmid.sam") as sam:
		plasmid_reads = [r for r in samReader(sam)]
		plasmid_read_ids = [r.qname for r in plasmid_reads]
	with open(f"{basename}_genome.sam") as sam:
		genome_reads = [r for r in samReader(sam)]
		genome_read_ids = [r.qname for r in genome_reads]

	for read in results:
		for end in read['ends']:
			pl_score = next((p.tags['NM'] for p in plasmid_reads if p.qname == end['id']), 1000)
			gn_score = next((g.tags['NM'] for g in genome_reads if g.qname == end['id']), 1000)
			if pl_score < gn_score:
				end['type'] = 'pl'
			elif pl_score > gn_score:
				end['type'] = 'gn'
			else:
				for i in range(len(read['hsps'])):
					gaps = []
					if i>=1:
						hit_gap = read['hsps'][i][1] - read['hsps'][i-1][1]
						gaps.append(hit_gap)
				if gaps and all([g<3450 and g>3420 for g in gaps]):
					end['type'] = 'pl'
				else:
					end['type'] = 'unknown'
		types = [e['type'] for e in read['ends']]
		if len(types) < 2:
			read['type'] = 'partialRead'
		elif len(set(types)) == 1 and types[0] == 'pl':
			read['type'] = 'pl'
		elif len(set(types)) == 1 and types[0] == 'gn':
			read['type'] = 'GENOME'
		elif 'unknown' in types:
			read['type'] = 'unknown'
		elif 'gn' in types and 'pl' in types:
			read['type'] = 'COINTEGRATE'
	return results

def get_read_obj(hit, tn_read, tn_length):
	read_result = {'id': hit.id, 'hsps': [], 'ends': [], 'result': None, 'len': hit.seq_len, 'read_seqrec': tn_read}

	prev_end = 0
	# each hsps represents a unique alignment in the read
	# should be ~equal to the transposon in length, or at start or end of the sequence
	for (i, hsp) in enumerate(sorted(hit.hsps, key=lambda x: x.hit_start)):
		hit_length = hsp.hit_end - hsp.hit_start
		if hit_length < (tn_length - 30) or hsp.evalue > 0.000001:
			continue
		read_result['hsps'].append((hsp.hit_start, hsp.hit_end))
		is_start = hsp.hit_start < 5
		is_end = hit.seq_len - hsp.hit_end < 5

		if not is_start:
			seqid = f'{hit.id}___{i}l'
			left_fp = tn_read.seq[max(prev_end, hsp.hit_start-40): hsp.hit_start]
			if len(left_fp) > 14:
				read_result['ends'].append({
					'id': seqid,
					'eval': hsp.pos_num,
					'seqrec': SeqRecord(left_fp, id=seqid, description=seqid, name=seqid)
				})

		if not is_end:
			seqid = f'{hit.id}___{i}r'
			next_start = hit.seq_len
			if i+1 < len(hit.hsps):
				next_start = sorted(hit.hsps, key=lambda x: x.hit_start)[i+1].hit_start
			right_fp = tn_read.seq[hsp.hit_end:min(next_start, hsp.hit_end+40)]
			if len(right_fp) > 14:
				read_result['ends'].append({
					'id': seqid,
					'eval': hsp.evalue,
					'seqrec': SeqRecord(right_fp, id=seqid, description=seqid, name=seqid)
				})
		prev_end = hsp.hit_end
	if len(read_result["hsps"]):
		return read_result

def do_blast(query_file, subject_file, output_name):
	#if Path(f'./{output_name}').exists():
	#	return output_name
	cline = NcbiblastnCommandline(query=query_file, subject=subject_file, num_alignments=10000, out=Path(f'./{output_name}'), outfmt=5)
	subprocess.run(str(cline), shell=True)
	return output_name

def do_bowtie(query_file, subject_file, output_name):
	cores = multiprocessing.cpu_count()
	name = subject_file.split('.')[0]
	build_cmd = f'bowtie2-build {subject_file} bt2index/{name} -q'
	subprocess.run(build_cmd, shell=True)
	align_command = f'bowtie2 -x bt2index/{name} --no-unal --omit-sec-seq -f {query_file} -p {cores-1} -S {output_name}'
	subprocess.run(align_command, shell=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)

if __name__ == "__main__":
	main()