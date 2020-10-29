from Bio import SeqIO, SearchIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
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

def hamming_dist(s1, s2):
    assert len(s1) == len(s2)
    return sum(ch1 != ch2 for ch1, ch2 in zip(s1, s2))

def main():
	os.makedirs(os.path.join(f"./outputs"), exist_ok=True)
	os.makedirs(os.path.join(f"./bt2index"), exist_ok=True)
	os.makedirs(os.path.join(f"./tmp"), exist_ok=True)
	today = date.today()
	with open(f'outputs/{output_name}.csv', 'w', newline='') as outfile:
		writer = csv.writer(outfile)
		writer.writerow(['Read_file', 'total_tn_reads', 'cointegrates', 'genomic_insertions', 'plasmids', 'insufficient', 'unknown', 'Sample Description', 'Uninterrupted insertion site reads', 'Normal site reads', 'Approx. Efficiency %', 'On-target %', 'multi_cointegrate_ct'])
	
	with open('input.csv', 'r', encoding='utf-8-sig') as infile:
		reader = csv.DictReader(infile)
		for row in reader:
			reads_file = row['Reads File']
			tn_file = row['Tn File']
			plasmid_file = row['Plasmid File']
			genome_file = row['Genome File']
			sample = row['Sample Code']
			sample_desc = row['Sample Desc']
			target = row['Target Sequence'].upper() if len(row['Target Sequence']) else None
			print(f"Processing {sample}...")
			print("-----")
			process_sample(reads_file, tn_file, plasmid_file, genome_file, sample, sample_desc, target)
			print("-----")
	
	if not run_local:
		s3 = boto3.client('s3')
		for filename in os.listdir('outputs'):
			s3key = f"cointegrate_outputs/{filename}"
			s3.upload_file(f"outputs/{filename}", "sternberg-sequencing-data", s3key)
	print("done")

def download_s3(filename, isNotSample=True):
	if Path(filename).exists() and run_local:
		return filename

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

def get_efficiency(reads_file, genome, target):
	# bowtie for the uninterrupted insertion site to get approximate integration efficiency
	uninterrupted_reads, other_site_reads, efficiency = None, None, None
	if target:
		site = genome.seq.upper().find(target)
		# add length of target, and the span of +35 to +65
		if site > 0:
			uninterrupted_insertion_site = genome.seq[site+len(target)+35:site+len(target)+65]
		else:
			site = genome.seq.upper().reverse_complement().find(target)
			uninterrupted_insertion_site = genome.seq.reverse_complement()[site+len(target)+35:site+len(target)+65]
		uninterrupted_fasta = f'tmp/uninterrupted.fasta'
		SeqIO.write(SeqRecord(uninterrupted_insertion_site, id="target"), uninterrupted_fasta, 'fasta')
		blast_filename =  f'tmp/uninterrupted_blastresults.xml'
		do_blast(uninterrupted_fasta, reads_file, blast_filename)
		# We are querying the 30bp sequence around the insertion site for WT sequence
		res = SearchIO.read(blast_filename, 'blast-xml')
		hits = [hit for hit in res.hits if len(hit.hsps) == 1 and hit.hsps[0].ident_num > 28]
		uninterrupted_reads = len(hits)

		# also do a blast for another random site in the genome to compare
		other_site = site - 5000 if site > 20000 else site + 5000
		SeqIO.write(SeqRecord(genome.seq[other_site:other_site+30], id='site'), 'tmp/other_site.fasta', 'fasta')
		do_blast('tmp/other_site.fasta', reads_file, 'tmp/other_site_blast.xml')
		res = SearchIO.read('tmp/other_site_blast.xml', 'blast-xml')
		hits = [hit for hit in res.hits if len(hit.hsps) == 1 and hit.hsps[0].ident_num > 28]
		other_site_reads = len(hits)

		efficiency = (other_site_reads - uninterrupted_reads) / other_site_reads

	return (uninterrupted_reads, other_site_reads, efficiency)

def process_sample(reads_file, tn_file, plasmid_file, genome_file, sample, sample_desc, target):
	reads_file = download_s3(reads_file, False)
	# reads_file = 'sample.fasta'
	tn_file = download_s3(tn_file)
	plasmid_file = download_s3(plasmid_file)
	genome_file = download_s3(genome_file)
	plasmid = SeqIO.read(plasmid_file, 'fasta')
	tn = SeqIO.read(tn_file, 'fasta')
	genome = SeqIO.read(genome_file, 'fasta')
	tn_length = len(tn)

	target_insertion_site = None
	if target:
		target_insertion_site = genome.seq.upper().find(target)
		if target_insertion_site >= 0:
			target_insertion_site += len(target) + 49
		else:
			target_insertion_site = genome.seq.upper().reverse_complement().find(target) + len(target) + 49
			target_insertion_site = len(genome.seq) - target_insertion_site

	tn_start = plasmid.seq.find(tn.seq)
	plasmid_l = plasmid.seq[tn_start-20:tn_start].upper()
	plasmid_r = plasmid.seq[tn_start+len(tn.seq):tn_start+len(tn.seq)+20].upper()
	plasmid_ends = [plasmid_l, plasmid_r]

	basename = "tmp/" + sample
	blast_filename =  f'{basename}_blastresults.xml'
	do_blast(tn_file, reads_file, blast_filename)
	
	# We are querying the transposon against all the reads, so only one result with many hits is output
	res = SearchIO.read(blast_filename, 'blast-xml')
	tn_read_ids = [hit.id for hit in res.hits]

	# filter the reads to just those with the transposon, write to file in case we want them later
	all_reads = SeqIO.parse(reads_file, 'fasta')
	tn_reads = [r for r in all_reads if r.id in tn_read_ids]
	#SeqIO.write(tn_reads, f'{basename}_tnreads.fasta', 'fasta')
	#tn_reads = list([r for r in SeqIO.parse(f'{basename}_tnreads.fasta', 'fasta')])

	efficiency_results = get_efficiency(reads_file, genome, target)

	all_results = []
	# each hit represents one or more matches of the query against a read sequence
	for hit in res.hits:
		tn_read = next(r for r in tn_reads if r.id == hit.id)
		read_obj = get_read_obj(hit, tn_read, tn_length)
		if read_obj:
			all_results.append(read_obj)

	attach_alignments(all_results, basename, plasmid_file, genome_file, plasmid_ends, target_insertion_site)
	
	with open(f'outputs/output_{sample}.csv', 'w', newline='') as outfile:
		writer = csv.writer(outfile)
		writer.writerow(['read_id', 'type', 'hsps', 'ends', 'types', 'read_length', 'genome_location', 'is_multi_coint'])
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
			location = read['genome_location'] if 'genome_location' in read else None

			writer.writerow([read['id'], read['type'], hsps, ends, types, read['len'], location, read['is_multi_coint'], read['strange_pl']])
	
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
		
		if target_insertion_site:
			reads_w_location = [r for r in all_results if 'genome_location' in r and r['type'] != 'unknown']
			ontarget_reads = [r for r in reads_w_location if abs(target_insertion_site - r['genome_location']) < 100]
		ontarget_perc = (len(ontarget_reads) / len(reads_w_location)) if (target_insertion_site and reads_w_location) else None
		multi_coint_ct = len([r for r in all_results if r['is_multi_coint']])
		strange_pl_ct = len([r for r in all_results if r['strange_pl'] is "TRUE"])
		writer.writerow([sample, len(all_results), len(cointegrates), len(genome_insertions), len(plasmids), len(insufficients), len(unknown), sample_desc, efficiency_results[0], efficiency_results[1], efficiency_results[2], ontarget_perc, multi_coint_ct, strange_pl_ct])
	print("done")

def get_short_end_type(end, read, plasmid_ends):
	end_id = end['id']
	other_end_id = end['id'][:-1]
	if end_id[-1] == 'l':
		other_end_id += 'r'
	else:
		other_end_id += 'l'
	other_end = next((e for e in read['ends'] if e['id'] == other_end_id), None)
	if not other_end:
		raise Exception("Can't find the other end for a short end")
	other_end_seq = other_end['seqrec'].seq.upper()

	end_seq = end['seqrec'].seq.upper()
	if end['rv']:
		end_seq = end_seq.reverse_complement()

	end_length = len(end_seq)
	dupe_length = min(end_length, 5)
	if other_end_id == 'l':
		end_dupe = other_end_seq[-dupe_length:]
	else:
		end_dupe = other_end_seq[0:dupe_length]
	
	plasmid_l = plasmid_ends[0][-end_length:]
	plasmid_r = plasmid_ends[1][:end_length]

	if hamming_dist(end_dupe, plasmid_r[:len(end_dupe)]) < 2 or hamming_dist(end_dupe, plasmid_l[-len(end_dupe):]) < 2:
		return 'unknown'
	if end_seq == plasmid_r or end_seq == plasmid_l:
		return 'pl'
	elif end_seq == end_dupe and other_end['type'] == 'gn':
		return 'gn'
	else:
		return 'unknown'

def attach_alignments(results, basename, plasmid_file, genome_file, plasmid_ends, target_insertion_site):
	long_fp_seqs = [end['seqrec'] for r in results for end in r['ends'] if len(end['seqrec']) > 14]
	short_fp_seqs = [end['seqrec'] for r in results for end in r['ends'] if len(end['seqrec']) <= 14]
	
	genome = SeqIO.read(genome_file, 'fasta')
	plasmid = SeqIO.read(plasmid_file, 'fasta')
	plasmid_length = len(plasmid)

	tmp_fp_fasta_name = f'{basename}_fps.fasta'
	SeqIO.write(long_fp_seqs, tmp_fp_fasta_name, 'fasta')
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
			if len(end['seqrec'].seq) > 14:
				pl_score = next((p.tags['NM'] for p in plasmid_reads if p.qname == end['id']), 1000)
				gn_reads = [g for g in genome_reads if g.qname == end['id']]
				
				min_nm = min([g.tags['NM'] for g in gn_reads]) if len(gn_reads) else 1000
				gn_read = [g for g in gn_reads if g.tags['NM'] == min_nm][0] if gn_reads else None
				gn_score = min_nm
				if pl_score < gn_score and int(pl_score) < 3:
					end['type'] = 'pl'
				elif pl_score > gn_score and int(gn_score) < 3:
					end['type'] = 'gn'
					# get the location of the read based on if it's to the left or right of the transposon, and the orientation of the read
					# the blast result flips the read, not the tn, so left is always BEFORE the raw tn sequence, not necessarily lower ref in the genome
					# so if it was a reverse read, we actually want the 'right' end (which is lower ref in the genome)
					is_left = end['id'][-1] is 'l'
					if (is_left and not gn_read.reverse) or (gn_read.reverse and not is_left):
						# no -1 because the +1 location site and bowtie 1-indexing cancel
						location = gn_read.coords[-1]
					else:
						# -1 because this read already starts AFTER the transposon and bowtie is 1-indexed
						location = gn_read.coords[0] - 1 + 5
					end['genome_location'] = location
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
		for end in read['ends']:
			if 'type' not in end:
				end['type'] = get_short_end_type(end, read, plasmid_ends)
		types = [e['type'] for e in read['ends']]

		locations = [e for e in read['ends'] if 'genome_location' in e]
		if len(locations) > 1:
			if len(locations[1]['seqrec'].seq) < 25 or (len(locations[1]['seqrec'].seq) > 25 and not locations[0]['rv']):
				read['genome_location'] = locations[0]['genome_location']
			else:
				read['genome_location'] = locations[1]['genome_location']
		if len(locations) == 1:
			read['genome_location'] = locations[0]['genome_location']

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

		if read['type'] is 'COINTEGRATE' and len([t for t in types if t == 'pl']) > 2:
			read['is_multi_coint'] = True
		else:
			read['is_multi_coint'] = False

		if read['type'] is 'pl' and len(read['seqrec']) > (len(ends)+2)/2 *plasmid_length:
			read['strange_pl'] = "TRUE"
		else:
			read['strange_pl'] = "F"
	return results

def get_read_obj(hit, tn_read, tn_length):
	read_result = {'id': hit.id, 'hit': hit, 'hsps': [], 'ends': [], 'result': None, 'len': hit.seq_len, 'read_seqrec': tn_read}
	prev_end = 0
	valid_hits = [hsp for hsp in hit.hsps if (hsp.hit_end - hsp.hit_start) > (tn_length - 5)]
	valid_hits = sorted(valid_hits, key=lambda x: x.hit_start)
	# each hsps represents a unique alignment in the read
	# should be ~equal to the transposon in length, or at start or end of the sequence
	for (i, hsp) in enumerate(valid_hits):
		hit_length = hsp.hit_end - hsp.hit_start
		if hsp.evalue > 0.000001:
			continue
		read_result['hsps'].append((hsp.hit_start, hsp.hit_end))
		is_start = hsp.hit_start < 4
		is_end = hit.seq_len - hsp.hit_end < 4

		if not is_start:
			seqid = f'{hit.id}___{i}l'
			left_fp = tn_read.seq[max(prev_end, hsp.hit_start-40): hsp.hit_start]
			read_result['ends'].append({
				'id': seqid,
				'rv': hsp.hit_strand == -1,
				'eval': hsp.evalue,
				'seqrec': SeqRecord(left_fp, id=seqid, description=seqid, name=seqid)
			})


		if not is_end:
			seqid = f'{hit.id}___{i}r'
			next_start = hit.seq_len
			if i+1 < len(valid_hits):
				next_start = sorted(valid_hits, key=lambda x: x.hit_start)[i+1].hit_start
			
			right_fp = tn_read.seq[hsp.hit_end:min(next_start, hsp.hit_end+40)]
			read_result['ends'].append({
				'id': seqid,
				'eval': hsp.evalue,
				'rv': hsp.hit_strand == -1,
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