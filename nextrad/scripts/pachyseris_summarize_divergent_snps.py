#!/usr/bin/env python
"""

"""
__author__ = 'Pim Bongaerts'
__copyright__ = 'Copyright (C) 2021 Pim Bongaerts'
__license__ = 'GPL'

import os
import sys
from glob import glob


def create_pachy2uniprot_dict(pachy2uniprot_filename, uniprot_filename):
	"""" Creates look-up dict with pachyseris gene ID -> uniprot gene/protein """ 
	uniprots = {}
	uniprot_file = open(uniprot_filename, 'r')
	for line in uniprot_file:
		cols = line.rstrip().split('\t')
		uniprots[cols[0]] = '{0},{1},{2},https://www.uniprot.org/uniprot/{0}'.format(cols[0], cols[2], cols[3])
	uniprot_file.close()

	pachy2uniprot = {}
	prev_gene = ''
	pachy2uniprot_file = open(pachy2uniprot_filename, 'r') # sorted by best match first
	for line in pachy2uniprot_file:
		cols = line.rstrip().split('\t')
		pachy_gene = cols[0]
		uniprot_id = cols[1]
		if pachy_gene != prev_gene:
			if uniprot_id in uniprots:
				pachy2uniprot[pachy_gene] = uniprots[uniprot_id]
				prev_gene = pachy_gene
			else:
				pachy2uniprot[pachy_gene] = 'NA,NA,NA,NA'

	pachy2uniprot_file.close()
	return pachy2uniprot

def create_pachy_snps_dict(pachy_snps_filename, pachy2uniprot):
	"""" Creates dict with pachyseris annotated pachyseris snps """ 
	pachy_snps = {}
	pachy_snps_file = open(pachy_snps_filename, 'r')
	for line in pachy_snps_file:
		cols = line.rstrip().split('\t')
		chrom_snp = '{0}_{1}'.format(cols[0], cols[1])
		pachy_snps[chrom_snp] = ','.join(cols[2:])
		# Include uniprot gene and protein name(s)
		pachy_gene1 = 'NA,NA,NA,NA'
		pachy_gene2 = 'NA,NA,NA,NA'
		if len(cols) > 5: 
			pachy_genes = cols[5].split('-')
			if pachy_genes[0] in pachy2uniprot:
				pachy_gene1 = pachy2uniprot[pachy_genes[0]]
			if len(pachy_genes) > 1 and pachy_genes[1] in pachy2uniprot:
				pachy_gene2 = pachy2uniprot[pachy_genes[1]]
	
			pachy_snps[chrom_snp] = '{0},{1},{2}'.format(pachy_snps[chrom_snp],
														   pachy_gene1,
														   pachy_gene2)
		else:
			pachy_gene1 = ''
			pachy_gene2 = ''
	pachy_snps_file.close()
	return pachy_snps

def create_outlier_dict(foldername):
	""" Append outlier classification to list of pachyseris SNPs """
	outlier_snps = {}
	outlier_filenames = glob('{0}/*.txt'.format(foldername))
	for outlier_filename in outlier_filenames:
		outlier_file = open(outlier_filename, 'r')
		outlier_type = 'out_{0}'.format(outlier_filename.split('_')[2])
		for line in outlier_file:
			cols = line.rstrip().split(' ') 
			chrom_snp = '{0}_{1}'.format(cols[0], cols[1])
			if chrom_snp in outlier_snps:
				outlier_snps[chrom_snp].append(outlier_type)
			else:
				outlier_snps[chrom_snp] = [outlier_type]
	return outlier_snps

def create_fixed_dict(foldername):
	""" Append fixed variant classification to list of pachyseris SNPs """
	fixed_snps = {}
	fixed_filenames = glob('{0}/*.csv'.format(foldername))
	for fixed_filename in fixed_filenames:
		fixed_file = open(fixed_filename, 'r')
		fixed_type = 'fix_{0}'.format(fixed_filename.split('_')[2].split('.')[0])
		for line in fixed_file:
			cols = line.rstrip().split(',') 
			chrom_snp = '{0}_{1}'.format(cols[0], cols[1])
			if 'AFD_OUTLIER' in cols[2]:
				if chrom_snp in fixed_snps:
					fixed_snps[chrom_snp].append(fixed_type)
				else:
					fixed_snps[chrom_snp] = [fixed_type]
	return fixed_snps

def output_table(pachy_snps, outlier_snps, fixed_snps):
	""" Output summary table """
	full_output_csv = open('pachy_divergent_snps.csv', 'w')
	headers = ['CHROM', 'POS', 'fixed', 'pcadapt', 'Allele', 'Annotation', 'putative_impact', 'gene(s)', 
			   'uniprot1_id', 'uniprot1_gene', 'uniprot1_protein', 'uniprot1_url', 
			   'uniprot2_id', 'uniprot2_gene', 'uniprot2_protein', 'uniprot2_url']
	full_output_csv.write(','.join(headers) + '\n')
	for chrom_pos in sorted(pachy_snps.keys()):
		if chrom_pos in outlier_snps:
			outlier_status = (';').join(outlier_snps[chrom_pos])
		else:
			outlier_status = ''
		if chrom_pos in fixed_snps:
			fixed_status = (';').join(fixed_snps[chrom_pos])
		else:
			fixed_status = ''
		if outlier_status or fixed_status:
			chrom_pos_formatted = chrom_pos.replace('_', ',')
			output_line = '{0},{1},{2},{3}\n'.format(chrom_pos_formatted,
													  fixed_status,
													  outlier_status,
													  pachy_snps[chrom_pos])
			full_output_csv.write(output_line)
	full_output_csv.close()

def main():
	# Create look-up table for pachy gene ID -> uniprot
	pachy2uniprot = create_pachy2uniprot_dict('pachy_gene2uniprot_b4c.txt', 'uniprot_names.tab')
	# Iterate over all annotated SNPs and uniprot identification
	pachy_snps = create_pachy_snps_dict('pachy_annotations_b4b.txt', pachy2uniprot)
	# Iterate over all alternatively fixed SNPs
	outlier_snps = create_outlier_dict('pcadapt')
	# Iterate over all alternatively fixed SNPs
	fixed_snps = create_fixed_dict('fixed')
	# Output summary table
	output_table(pachy_snps, outlier_snps, fixed_snps)

if __name__ == '__main__':
	main()