#!/usr/bin/env python3

## parse_ARCHER.py - a tool for manipulating ARCHER-DX summarised vcfs
## Usage: 
## Author: @neil.alistair.robertson@hotmail.co.uk

import getopt, sys, os
import csv
import resource

DELIMITER = "\t"
MISSING_VALUE = "NA"
GIAB_FILES = ("GenomeInABottle1","GenomeInABottle2","GenomeInABottle3","GenomeInABottle4","GenomeInABottle5","GenomeInABottle6","GenomeInABottle7","GenomeInABottle8","Genome-in-a-Bottle-A","Genome-in-a-Bottle-B","GenomeInABottle_C","GenomeInABottle_D","Genome_in_a_Bottle_E","Genome_in_a_Bottle_F")
BASE_DIRECTORY = "/mnt/tchandra-lab/Neil/data/DNA-seq/LBC_ARCHER_Panel-seq/data"

CONSEQUENCE_FILTER = ("5_prime_UTR_variant", "coding_sequence_variant", "feature_elongation", "feature_truncation", "frameshift_variant", "incomplete_terminal_codon_variant", "inframe_deletion", "inframe_insertion", "missense_variant", "protein_altering_variant", "start_lost", "stop_gained", "stop_lost", "transcript_ablation", "transcript_amplification", "splice_acceptor_variant", "splice_donor_variant", "splice_region_variant")
AO_FILTER = 5
UAO_FILTER = 3
GNOMAD_AF_FILTER = 0.05
SYNONYMOUS_AF_GERMLINE_LIMIT = 0.45
AF_FILTER = 0.005
SYNONYMOUS_AF_FILTER = 0.005
SEQ_BIAS_FILTER = 'Yes'

from itertools import groupby as g
def most_common_consequence(consequences_list):
	HIGH_IMPACT = ("transcript_ablation","splice_acceptor_variant","splice_donor_variant","stop_gained","frameshift_variant","stop_lost","start_lost","transcript_amplification")
	MODERATE_IMPACT = ("inframe_insertion","inframe_deletion","missense_variant","protein_altering_variant","regulatory_region_ablation")
	LOW_IMPACT = ("splice_region_variant","incomplete_terminal_codon_variant","start_retained_variant","stop_retained_variant","synonymous_variant")
	MODIFIER_IMPACT = ("coding_sequence_variant","mature_miRNA_variant","5_prime_UTR_variant","3_prime_UTR_variant","non_coding_transcript_exon_variant","intron_variant","NMD_transcript_variant","non_coding_transcript_variant","upstream_gene_variant","downstream_gene_variant","TFBS_ablation","TFBS_amplification","TF_binding_site_variant","regulatory_region_amplification","feature_elongation","regulatory_region_variant","feature_truncation","intergenic_variant")
	if len(list(set(consequences_list).intersection(HIGH_IMPACT))) > 0:
		L = intersection(HIGH_IMPACT, consequences_list)
		return max(g(sorted(L)), key=lambda(x, v):(len(list(v)),-L.index(x)))[0]
	elif len(list(set(consequences_list).intersection(MODERATE_IMPACT))) > 0:
		L = intersection(MODERATE_IMPACT, consequences_list)
		return max(g(sorted(L)), key=lambda(x, v):(len(list(v)),-L.index(x)))[0]
	elif len(list(set(consequences_list).intersection(MODIFIER_IMPACT))) > 0:
		L = intersection(MODIFIER_IMPACT, consequences_list)
		return max(g(sorted(L)), key=lambda(x, v):(len(list(v)),-L.index(x)))[0]
	elif len(list(set(consequences_list).intersection(LOW_IMPACT))) > 0:
		L = intersection(LOW_IMPACT, consequences_list)
		return max(g(sorted(L)), key=lambda(x, v):(len(list(v)),-L.index(x)))[0]
	return max(g(sorted(L)), key=lambda(x, v):(len(list(v)),-L.index(x)))[0]

def intersection(lst1, lst2): 
	return [value for value in lst1 if value in lst2] 

def parse_manifest(manifest_filename):
	print("-> Beginning to parse mainifest file: {0}".format(manifest_filename))
	manifest_keys = []
	sample_wave_dict = {}
	count = 0
	with open(manifest_filename, mode='r') as csv_file:
		csv_reader = csv.DictReader(csv_file, delimiter = DELIMITER)
		for i, row in enumerate(csv_reader):
			if i == 0: 
				manifest_keys = row.keys()
			else:
				participant_id = row['Particiapant_ID']
				if participant_id not in GIAB_FILES:
					if participant_id not in sample_wave_dict:
						sample_wave_dict[participant_id] = {}
					wave = row['Wave']
					sample_wave_dict[participant_id][wave] = row
	return manifest_keys, sample_wave_dict
	#OrderedDict([('Original_Filename', 'LBC361225_03159A3000121_S48_R1_001.fastq.gz'), ('Sequencing_Batch', 'FC04'), ('Particiapant_ID', 'LBC361225'), ('Sample_ID', '03159A3000121'), ('Batch_Sample_Number', 'S48'), ('Paired', 'Paired'), ('Platform', 'Illumina'), ('Run_Folder', 'FC04_4708'), ('Wave', '1'), ('Cohort', 'LBC1936'), ('Raw_Variant_Count', '18970')])


def build_vcf_filename(sample_metadata):
	file_prefix = sample_metadata['Original_Filename'].split(".fastq.gz")[0]
	vcf_filename = "{0}/{1}/{2}.vcf.summary.tsv".format(BASE_DIRECTORY, sample_metadata['Run_Folder'], file_prefix)
	return vcf_filename


def write_file(output_directory, filename, header, rows):
	file_path = "{0}/{1}".format(output_directory, filename)
	with open(file_path, 'w') as output_file:
		csv_writer = csv.writer(output_file, delimiter = DELIMITER)
		csv_writer.writerow(header)
		for row in rows:
			csv_writer.writerow(row)


def make_directory(dir_path):
	if os.path.isdir(dir_path) is True:
		print("-> Directory {0} exists...".format(dir_path))
	else:
		os.mkdir(dir_path)
		print("-> Created directory: {0} \n".format(dir_path))


def parse_vcf_consequences(consequences_string):
	consequences_list = consequences_string.strip().split("|")
	consequences_list = [x.strip().split("&") for x in consequences_list]
	consequences_list = list(set([x[0] for x in consequences_list]))
	return consequences_list


def parse_transcript_id(transcript_string):
	transcript_parts = transcript_string.strip().split("|")
	for trans in transcript_parts:
		if trans.find("NM_") != -1:
			if len(trans) == 11:
				return trans


def parse_HGVSc(hgvc_string):
	hgvc_parts = hgvc_string.strip().split("|")
	for part in hgvc_parts:
		if part.find("NM_") != -1:
			if len(part.split(".")[0]) == 9:
				return part


def parse_HGVSp(hgvp_string):
	hgvp_parts = hgvp_string.strip().split("|")
	for part in hgvp_parts:
		if part.find("NP_") != -1:
			if len(part.split(".")[0]) == 9:
				return part


def parse_substitution(hgvs_string):
	substitution = ''
	if hgvs_string not in (None, ''):
		substitution = hgvs_string.strip().split(":")[-1]
	return substitution


def get_simplified_header():
	return ["PreferredSymbol", "type", "HGVSp", "protein_substitution", "HGVSc", "base_substitution", "DP", "AO", "UAO", "AF", "participant_id", "wave", "chromosome", "position", "reference", "mutation", "gnomAD_AF", "consequence", "annotation", "transcript_id", "quality", "FunctionalStatus", "COSMICID", "SIFT", "PolyPhen", "HasSeqDirBias", "PosPrimerRefCount", "NegPrimerRefCount", "PosPrimerAltCount", "NegPrimerAltCount", "UDP", "DAO", "URO", "UAO", "SampleStrandBiasRatio", "SampleStrandBiasProb", "Min_Outlier_PValue", "AF_Outlier_Pvalue", "95MDAF"]


def simplify_vcf_row(row, participant_id, wave):
	preferred_symbol = row['PreferredSymbol']
	mtype = row['type']
	HGVSp = row['HGVSp']
	HGVSp = parse_HGVSp(HGVSp)
	protein_substitution = parse_substitution(HGVSp)
	HGVSc = row['HGVSc']
	HGVSc = parse_HGVSc(HGVSc)
	base_substitution = parse_substitution(HGVSc)

	DP = row['DP']
	AO = row['AO']
	UAO = row['UAO']
	AF = row['AF']

	chromosome = row['chromosome']
	position = row['position']
	reference = row['reference']
	mutation = row['mutation']

	gnomAD_AF = row['gnomAD_AF']
	consequence = row['consequence']
	consequence = parse_vcf_consequences(consequence)
	consequence = most_common_consequence(consequence)

	annotation = row['annotation']
	transcript_id = row['transcript_id']
	transcript_id = parse_transcript_id(transcript_id)

	quality = row['quality']

	FunctionalStatus = row['FunctionalStatus']
	COSMICID = row['COSMICID']
	SIFT = row['SIFT'].split("|")[0]
	PolyPhen = row['PolyPhen'].split("|")[0]

	seq_bias = row['HasSeqDirBias']

	PosPrimerRefCount = row['PosPrimerRefCount']
	NegPrimerRefCount = row['NegPrimerRefCount']
	PosPrimerAltCount = row['PosPrimerAltCount']
	NegPrimerAltCount = row['NegPrimerAltCount']
	UDP = row['UDP']
	DAO = row['DAO']
	URO = row['URO']
	UAO = row['UAO']
	SampleStrandBiasRatio = row['SampleStrandBiasRatio']
	SampleStrandBiasProb = row['SampleStrandBiasProb']
	
	Min_Outlier_PValue = row['Min_Outlier_PValue']
	AF_Outlier_Pvalue = row['AF_Outlier_Pvalue']
	MDAF95 = row['95MDAF']

	output_row = [preferred_symbol, mtype, HGVSp, protein_substitution, HGVSc, base_substitution, DP, AO, UAO, AF, participant_id, wave, chromosome, position, reference, mutation, gnomAD_AF, consequence, annotation, transcript_id, quality, FunctionalStatus, COSMICID, SIFT, PolyPhen, seq_bias, PosPrimerRefCount, NegPrimerRefCount, PosPrimerAltCount, NegPrimerAltCount, UDP, DAO, URO, UAO, SampleStrandBiasRatio, SampleStrandBiasProb, Min_Outlier_PValue, AF_Outlier_Pvalue, MDAF95]
	return output_row


def simplify_vcf_summary(sample_vcf_filename, participant_id, wave):
	vcf_keys = None

	if os.path.isfile(sample_vcf_filename) != True:
		print("Cannot find file: {0}".format(sample_vcf_filename))
		exit()
	else:
		print("-> Opening file: {0}".format(sample_vcf_filename))
		with open(sample_vcf_filename, 'r') as vcf_file:
			vcf_reader = csv.DictReader(vcf_file, delimiter = DELIMITER)
			for i, row in enumerate(vcf_reader):
				if i == 0: 
					vcf_keys = row.keys()
					#print(vcf_keys)
				else: 
					output_row = simplify_vcf_row(row, participant_id, wave)
					yield output_row


def filter_vcf_summary(sample_vcf_filename, participant_id, wave):
	vcf_keys = None

	if os.path.isfile(sample_vcf_filename) != True:
		print("Cannot find file: {0}".format(sample_vcf_filename))
		exit()
	else:
		print("-> Opening file: {0}".format(sample_vcf_filename))
		with open(sample_vcf_filename, 'r') as vcf_file:
			vcf_reader = csv.DictReader(vcf_file, delimiter = DELIMITER)
			for i, row in enumerate(vcf_reader):
				if i == 0: 
					vcf_keys = row.keys()
					#print(vcf_keys)
				else: 
					AO = row['AO']
					UAO = row['UAO']
					AF = row['AF']
					if int(AO) >= AO_FILTER:
						if int(UAO) >= UAO_FILTER:
							if float(AF) >= AF_FILTER:
								gnomAD_AF = row['gnomAD_AF']
								if (gnomAD_AF == '') or (float(gnomAD_AF) <= GNOMAD_AF_FILTER):

									seq_bias = row['HasSeqDirBias']
									if seq_bias.strip() != SEQ_BIAS_FILTER:

										consequence = row['consequence']
										consequences_list = parse_vcf_consequences(consequence)
										if list(set(consequences_list).intersection(CONSEQUENCE_FILTER)):
											output_row = simplify_vcf_row(row, participant_id, wave)
											yield output_row


def filter_for_synonymous(sample_vcf_filename, participant_id, wave):
	vcf_keys = None

	if os.path.isfile(sample_vcf_filename) != True:
		print("Cannot find file: {0}".format(sample_vcf_filename))
		exit()
	else:
		print("-> Opening file: {0}".format(sample_vcf_filename))
		with open(sample_vcf_filename, 'r') as vcf_file:
			vcf_reader = csv.DictReader(vcf_file, delimiter = DELIMITER)
			for i, row in enumerate(vcf_reader):
				if i == 0: 
					vcf_keys = row.keys()
					#print(vcf_keys)
				else: 
					AO = row['AO']
					UAO = row['UAO']
					AF = row['AF']
					if int(AO) >= AO_FILTER:
						if int(UAO) >= UAO_FILTER:
							if float(AF) >= SYNONYMOUS_AF_FILTER:
								if float(AF) <= SYNONYMOUS_AF_GERMLINE_LIMIT:
									#gnomAD_AF = row['gnomAD_AF']
									#if (gnomAD_AF == '') or (float(gnomAD_AF) <= GNOMAD_AF_FILTER):
									seq_bias = row['HasSeqDirBias']
									if seq_bias.strip() != SEQ_BIAS_FILTER:

										consequence = row['consequence']
										consequences_list = parse_vcf_consequences(consequence)
										if list(set(consequences_list).intersection(['synonymous_variant'])):
											output_row = simplify_vcf_row(row, participant_id, wave)
											yield output_row


def build_mutation_dictionary(sample_wave_dict, participant_id, file_directory):
	summarized_variant_dict = {}
	waves = []
	for wave in sorted(sample_wave_dict[participant_id].keys()):
		sample_metadata = sample_wave_dict[participant_id][wave]
		filename = "{0}/{1}.wave_{2}.variants.summary.tsv".format(file_directory, sample_metadata['Original_Filename'].split(".fastq.gz")[0], wave)
		with open(filename, 'r') as summarised_file:
			vcf_reader = csv.DictReader(summarised_file, delimiter = DELIMITER)
			for i, row in enumerate(vcf_reader):
				if i == 0: pass
				else:
					key = row['HGVSc']
					if key not in (None, ''):
						if key not in summarized_variant_dict:
							summarized_variant_dict[key] = {}
						wave = row['wave']
						waves.append(wave)
						summarized_variant_dict[key][wave] = row
	return summarized_variant_dict, sorted(set(waves))


def build_sample_tables(sample_metadata, participant_id):
	pass


options, remainder = getopt.getopt(sys.argv[1:], 'm:o', ['manifest=','output_directory='])
manifest_filename = None
output_directory = None

for opt, arg in options:
	if opt in ("-m", "--manifest"):
		manifest_filename = arg.strip()
	elif opt in ("-o", "--output_directory"):
		output_directory = arg.strip()

assert manifest_filename != None, ("Please provide participant metadata manifest: --manifest")
assert output_directory != None, ("Please provide output directory for processed files: --output_directory")

print("MANIFEST: {0}".format(manifest_filename))
print("OUTPUT DIRECTORY: {0} \n".format(output_directory))
print("VAF LIMIT: {0} \n".format(AF_FILTER))

manifest_keys, sample_wave_dict = parse_manifest(manifest_filename)

output_directory = output_directory.rstrip("/")
make_directory(output_directory)

run_start = True
if run_start:
	simplified_directory = "{0}/simplified_variants".format(output_directory)
	make_directory(simplified_directory)
	print("BEGINNING TO SIMPLIFY AND WRITE ALL VARIANTS...")
	for counter, participant_id in enumerate(sample_wave_dict.keys()):
	#for participant_id in ["LBC360021","LBC360072","LBC360508","LBC360549","LBC360702","LBC360725"]:
		for wave in sorted(sample_wave_dict[participant_id].keys()):
			sample_metadata = sample_wave_dict[participant_id][wave]
			sample_vcf_filename = build_vcf_filename(sample_metadata)

			output_rows = []
			filename = sample_metadata['Original_Filename'].split(".fastq.gz")[0] + ".wave_{0}.variants.summary.tsv".format(wave)
			for summarised_row in simplify_vcf_summary(sample_vcf_filename, participant_id, wave):
				output_rows.append(summarised_row)
			print("--> Writing simplified file for {0} at wave {1}".format(participant_id, wave))
			write_file(simplified_directory, filename, get_simplified_header(), output_rows)

if run_start:
	filter_directory = "{0}/filtered_variants".format(output_directory)
	make_directory(filter_directory)
	print("BEGINNING TO FILTER AND WRITE VARIANTS...")
	for counter, participant_id in enumerate(sample_wave_dict.keys()):
	#for participant_id in ["LBC360021","LBC360072","LBC360508","LBC360549","LBC360702","LBC360725"]:
		for wave in sorted(sample_wave_dict[participant_id].keys()):
			sample_metadata = sample_wave_dict[participant_id][wave]
			sample_vcf_filename = build_vcf_filename(sample_metadata)

			found = 0
			output_rows = []
			filename = sample_metadata['Original_Filename'].split(".fastq.gz")[0] + ".wave_{0}.variants.summary.tsv".format(wave)
			for summarised_row in filter_vcf_summary(sample_vcf_filename, participant_id, wave):
				output_rows.append(summarised_row)
				found += 1
			print("--> Writing summarized file  for {0} at wave {1}. Variants found: {2}".format(participant_id, wave, str(found)))
			write_file(filter_directory, filename, get_simplified_header(), output_rows)


syn_run = True
if syn_run:
	synonymous_directory = "{0}/synonymous_variants".format(output_directory)
	make_directory(synonymous_directory)
	print("BEGINNING TO FILTER AND WRITE SYNONYMOUS VARIANTS...")
	for counter, participant_id in enumerate(sample_wave_dict.keys()):
	#for participant_id in ["LBC360021","LBC360072","LBC360508","LBC360549","LBC360702","LBC360725"]:
		for wave in sorted(sample_wave_dict[participant_id].keys()):
			sample_metadata = sample_wave_dict[participant_id][wave]
			sample_vcf_filename = build_vcf_filename(sample_metadata)

			found = 0
			output_rows = []
			filename = sample_metadata['Original_Filename'].split(".fastq.gz")[0] + ".wave_{0}.variants.summary.tsv".format(wave)
			for summarised_row in filter_for_synonymous(sample_vcf_filename, participant_id, wave):
				output_rows.append(summarised_row)
				found += 1
			print("--> Writing summarized file  for {0} at wave {1}. Variants found: {2}".format(participant_id, wave, str(found)))
			write_file(synonymous_directory, filename, get_simplified_header(), output_rows)


run_processing = True
if run_processing:
	table_directory = "{0}/processed_tables".format(output_directory)
	make_directory(table_directory)
	processed_directory = "{0}/processed_variants".format(output_directory)
	make_directory(processed_directory)
	print("BUILDING TABLES FROM DATA...")
	for counter, participant_id in enumerate(sample_wave_dict.keys()):
		output_rows = []
		table_rows = []
		accepted_variant_dict, waves = build_mutation_dictionary(sample_wave_dict, participant_id, "{0}/filtered_variants".format(output_directory))
		all_variant_dict, waves = build_mutation_dictionary(sample_wave_dict, participant_id, "{0}/simplified_variants".format(output_directory))
		for mutation_key in accepted_variant_dict.keys():
			allele_fractions = []
			is_wave = None
			for wave in waves:
				if wave in accepted_variant_dict[mutation_key].keys():
					is_wave = wave
					AF = accepted_variant_dict[mutation_key][wave]['AF']
					allele_fractions.append(AF)
					output_rows.append(simplify_vcf_row(accepted_variant_dict[mutation_key][wave], participant_id, wave))
				elif (mutation_key in all_variant_dict.keys()) and (wave in all_variant_dict[mutation_key].keys()):
					AF = all_variant_dict[mutation_key][wave]['AF']
					allele_fractions.append(AF)
					output_rows.append(simplify_vcf_row(all_variant_dict[mutation_key][wave], participant_id, wave))
					#print ("Missing value found at wave: {0} mutation: {1}".format(wave, mutation_key))
				else:
					allele_fractions.append(MISSING_VALUE)
			#print("Mutation: {0}    AFs: {1}".format(mutation_key, allele_fractions))
			row = accepted_variant_dict[mutation_key][is_wave]
			table_row = [row["PreferredSymbol"]] + allele_fractions + [row["HGVSp"], row["protein_substitution"], row["HGVSc"], row["base_substitution"], row["participant_id"], row["chromosome"], row["position"], row["reference"], row["mutation"], row["gnomAD_AF"], row["consequence"], row["annotation"], row["transcript_id"], row["FunctionalStatus"], row["COSMICID"], row["SIFT"], row["PolyPhen"]]
			table_rows.append(table_row)
		print("--> Writing processed file and table for {0} at all waves".format(participant_id))
		filename = "{0}.processed_variants.summary.tsv".format(participant_id)
		write_file(processed_directory, filename, get_simplified_header(), output_rows)

		filename = "{0}.processed_variants.summary_table.tsv".format(participant_id)
		header = ["PreferredSymbol"] + ["AF_wave{0}".format(x) for x in waves] + ["HGVSp", "protein_substitution", "HGVSc", "base_substitution", "participant_id", "chromosome", "position", "reference", "mutation", "gnomAD_AF", "consequence", "annotation", "transcript_id", "FunctionalStatus", "COSMICID", "SIFT", "PolyPhen"]
		write_file(table_directory, filename, header, table_rows)


	synonymous_table_directory = "{0}/synonymous_tables".format(output_directory)
	make_directory(synonymous_table_directory)
	processed_syn_directory = "{0}/processed_synonymous_variants".format(output_directory)
	make_directory(processed_syn_directory)
	print("BUILDING TABLES FROM DATA...")
	for counter, participant_id in enumerate(sample_wave_dict.keys()):
		output_rows = []
		table_rows = []
		accepted_variant_dict, waves = build_mutation_dictionary(sample_wave_dict, participant_id, "{0}/synonymous_variants".format(output_directory))
		all_variant_dict, waves = build_mutation_dictionary(sample_wave_dict, participant_id, "{0}/simplified_variants".format(output_directory))
		for mutation_key in accepted_variant_dict.keys():
			allele_fractions = []
			is_wave = None
			for wave in waves:
				if wave in accepted_variant_dict[mutation_key].keys():
					is_wave = wave
					AF = accepted_variant_dict[mutation_key][wave]['AF']
					allele_fractions.append(AF)
					output_rows.append(simplify_vcf_row(accepted_variant_dict[mutation_key][wave], participant_id, wave))
				elif (mutation_key in all_variant_dict.keys()) and (wave in all_variant_dict[mutation_key].keys()):
					AF = all_variant_dict[mutation_key][wave]['AF']
					allele_fractions.append(AF)
					output_rows.append(simplify_vcf_row(all_variant_dict[mutation_key][wave], participant_id, wave))
					#print ("Missing value found at wave: {0} mutation: {1}".format(wave, mutation_key))
				else:
					allele_fractions.append(MISSING_VALUE)
			#print("Mutation: {0}    AFs: {1}".format(mutation_key, allele_fractions))
			row = accepted_variant_dict[mutation_key][is_wave]
			table_row = [row["PreferredSymbol"]] + allele_fractions + [row["HGVSp"], row["protein_substitution"], row["HGVSc"], row["base_substitution"], row["participant_id"], row["chromosome"], row["position"], row["reference"], row["mutation"], row["gnomAD_AF"], row["consequence"], row["annotation"], row["transcript_id"], row["FunctionalStatus"], row["COSMICID"], row["SIFT"], row["PolyPhen"]]
			table_rows.append(table_row)
		print("--> Writing processed file and table for {0} at all waves".format(participant_id))
		filename = "{0}.processed_synonymous_variants.summary.tsv".format(participant_id)
		write_file(processed_syn_directory, filename, get_simplified_header(), output_rows)

		filename = "{0}.synonymous_variants.summary_table.tsv".format(participant_id)
		header = ["PreferredSymbol"] + ["AF_wave{0}".format(x) for x in waves] + ["HGVSp", "protein_substitution", "HGVSc", "base_substitution", "participant_id", "chromosome", "position", "reference", "mutation", "gnomAD_AF", "consequence", "annotation", "transcript_id", "FunctionalStatus", "COSMICID", "SIFT", "PolyPhen"]
		write_file(synonymous_table_directory, filename, header, table_rows)








