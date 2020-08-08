from Bio import SeqIO, SearchIO
from Bio.SeqRecord import SeqRecord
from Bio.Blast.Applications import NcbiblastnCommandline
import multiprocessing
import subprocess
import os
import csv
import boto3
import botocore
from pathlib import Path
from simplesam import Reader as samReader

def main():
	os.makedirs(os.path.join(f"./outputs"), exist_ok=True)
	os.makedirs(os.path.join(f"./bt2index"), exist_ok=True)
	os.makedirs(os.path.join(f"./tmp"), exist_ok=True)
	with open('outputs/all.csv', 'w', newline='') as outfile:
		writer = csv.writer(outfile)
		writer.writerow(['Read_file', 'total_tn_reads', 'cointegrates', 'genomic_insertions', 'plasmids', 'insufficient', 'unknown'])
	'''
	reads_file = 'CCS_SLMS_23.fasta'
	tn_file = 'transposon_SLMS_23.fasta'
	plasmid_file = 'transposon_plasmid_SLMS23.fasta'
	genome_file = 'bl21de3-genome.fasta'
	'''
	print("1")
	process_sample('CCS_SLMS_22.fasta', 'SLMS_22_tn.fasta', 'SLMS_22_plasmid.fasta', 'bl21de3-genome.fasta')
	print("2")
	process_sample('CCS_SLMS_23.fasta', 'transposon_SLMS_23.fasta', 'transposon_plasmid_SLMS23.fasta', 'bw25113-reca-pb-genome.fasta')
	print("3")
	process_sample('CCS_SLMS_24.fasta', 'SLMS_24_tn.fasta', 'SLMS_24_plasmid.fasta', 'bw25113-reca-pb-genome.fasta')
	print("4")
	process_sample('CCS_SLMS_3.fasta', 'SLMS_3_tn.fasta', 'SLMS_3_plasmid.fasta', 'bl21de3-genome.fasta')
	print("uploading outputs")
	s3 = boto3.client('s3')
	for filename in os.listdir('outputs'):
		s3key = f"pilot_samples/outputs/{filename}"
		s3.upload_file(f"outputs/{filename}", "sternberg-sequencing-data", s3key)
	print("done")

def download_s3(filename):
	if Path(filename).exists():
		return
	s3 = boto3.client('s3')
	s3.download_file('sternberg-sequencing-data', f'pilot_samples/{filename}', filename)

def process_sample(reads_file, tn_file, plasmid_file, genome_file):
	download_s3(reads_file)
	download_s3(tn_file)
	download_s3(plasmid_file)
	download_s3(genome_file)
	plasmid = SeqIO.read(plasmid_file, 'fasta')
	tn = SeqIO.read(tn_file, 'fasta')
	genome = SeqIO.read(genome_file, 'fasta')

	tn_start = plasmid.seq.find(tn.seq)
	plasmid_l = plasmid.seq[tn_start-20:tn_start]
	plasmid_r = plasmid.seq[tn_start+len(tn.seq):tn_start+len(tn.seq)+20]

	basename = "tmp/" + reads_file.split('.')[0]
	blast_filename =  f'{basename}_blastresults.xml'
	do_blast(tn_file, reads_file, blast_filename)
	
	# We are querying the transposon against all the reads, so only one result with many hits is output
	res = SearchIO.read(blast_filename, 'blast-xml')
	tn_read_ids = [hit.id for hit in res.hits]

	# filter the reads to just those with the transposon, write to file in case we want them later
	all_reads = SeqIO.parse(reads_file, 'fasta')
	tn_reads = [r for r in all_reads if r.id in tn_read_ids]
	SeqIO.write(tn_reads, f'{basename}_tnreads.fasta', 'fasta')
	#tn_reads = list([r for r in SeqIO.parse(f'{basename}_tnreads.fasta', 'fasta')])

	all_results = []
	# each hit represents one or more matches of the query against a read sequence
	for hit in res.hits:
		tn_read = next(r for r in tn_reads if r.id == hit.id)
		all_results.append(get_read_obj(hit, tn_read))
	attach_alignments(all_results, basename, plasmid_file, genome_file)
	with open(f'outputs/output_{reads_file}.csv', 'w', newline='') as outfile:
		writer = csv.writer(outfile)
		writer.writerow(['read_id', 'type', 'hsps', 'ends', 'types'])
		for read in all_results:
			hsps = ['-'.join([str(i) for h in read['hsps'] for i in h])]
			ends = [e['id'][-19:] for e in read['ends']]
			types = [e['type'] for e in read['ends']]
			writer.writerow([read['id'], read['type'], hsps, ends, types])
	with open(f'outputs/all.csv', 'a', newline='') as outfile:
		writer = csv.writer(outfile)
		cointegrates = [r for r in all_results if r['type'] is 'COINTEGRATE']
		SeqIO.write([r['read_seqrec'] for r in cointegrates[:10]], f"outputs/{reads_file}_cointegrates.fasta", "fasta")
		genome_insertions = [r for r in all_results if r['type'] is 'GENOME']
		SeqIO.write([r['read_seqrec'] for r in genome_insertions[:10]], f"outputs/{reads_file}_genomic.fasta", "fasta")
		plasmids = [r for r in all_results if r['type'] is 'pl']
		SeqIO.write([r['read_seqrec'] for r in plasmids[:10]], f"outputs/{reads_file}_plasmids.fasta", "fasta")
		insufficients = [r for r in all_results if r['type'] is 'partialRead']
		print(set([r['type'] for r in all_results]))
		SeqIO.write([r['read_seqrec'] for r in insufficients[:10]], f"outputs/{reads_file}_insufficient.fasta", "fasta")
		unknown = [r for r in all_results if r['type'] is 'unknown']
		SeqIO.write([r['read_seqrec'] for r in unknown[:10]], f"outputs/{reads_file}_unknowns.fasta", "fasta")
		writer.writerow([reads_file, len(all_results), len(cointegrates), len(genome_insertions), len(plasmids), len(insufficients), len(unknown)])
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

def get_read_obj(hit, tn_read):
	read_result = {'id': hit.id, 'hsps': [], 'ends': [], 'result': None, 'len': hit.seq_len, 'read_seqrec': tn_read}

	prev_end = 0
	# each hsps represents a unique alignment in the read
	# should be ~equal to the transposon in length, or at start or end of the sequence
	for (i, hsp) in enumerate(sorted(hit.hsps, key=lambda x: x.hit_start)):
		read_result['hsps'].append((hsp.hit_start, hsp.hit_end))
		hit_length = hsp.hit_end - hsp.hit_start
		if hit_length < 100:
			continue
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
	return read_result

def do_blast(query_file, subject_file, output_name):
	if Path(f'./{output_name}').exists():
		return output_name
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