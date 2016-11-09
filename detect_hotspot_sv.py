import os
import argparse
import shlex
import subprocess
import vcf, vcf.model
import pandas as pd
from iCallSV.dellyVcf2Tab import vcf2tab

'''Filter a VCF for PRECISE deletions with PE > 10

Args:
	vcf_reader: vcf.Reader object

Yields:
	vcf._Record: the calls that pass filter

'''
def filter_vcf(vcf_reader):
	for record in vcf_reader:
		# assert type(record) == vcf.model._Record

		info = record.INFO

		if (info['SVTYPE'] == 'DEL') and ('PRECISE' in info) and (info['PE'] > 10):
			yield record

'''Find pseudogenes herustically.

Counts the number of events that aren't in exons, "in frame" or "out of frame",
and returns the genes with more at least ``threshold`` of these events.
Adapted from check_cDNA_contamination.py by Ronak Shah.

Args:
	dataDF (pandas.DataFrame): Output of iAnnotateSV
	THRESHOLD (int): Threshold number of qualifying events

Returns:
	List of pseudogenes
'''
def find_pseudogenes(dataDF, THRESHOLD):
	# Group the data by gene1 name
	gDF = dataDF.groupby('gene1').groups
	
	geneList = []
	# traverse through gDF dictionary 
	for gene1, value in gDF.iteritems():
		# check how many entries are there of sv type
		entries = len(value)
		# run only if the event is deletion and has more then 2 entries for the same gene
		if entries >= THRESHOLD:
			
			# number of cDNA events
			count = 0
			for idx in value:

				record = dataDF.loc[idx]
				site1 = str(record.loc['site1'])
				site2 = str(record.loc['site2'])
				fusion = str(record.loc['fusion'])
				gene2 = record.loc['gene2']
				
				# Skip entries that are:
				# - within exon
				# - in-frame or out-of-frame
				# - involve IGRs.
				# - involve multiple genes
				if ("Exon" in site1 and "Exon" in site2) or \
				   ("in frame" in fusion or "out of frame" in fusion) or \
				   ('IGR' in site1 or 'IGR' in site2) or \
				   (gene1 != gene2):
					continue
				
				else:
					count += 1

			# count the entries,genes and fill the 
			if count >= THRESHOLD:
				geneList.append(gene1)

	return geneList


'''Main method

Passing a string will parse command line args from the string instead of stdin
'''
def main(command=None):
	
	# 0: Parse aand interpret args

	parser = argparse.ArgumentParser()
	# Required args
	parser.add_argument('in_vcf', 
		action='store', 
		metavar='structural_variants.vcf',
		help='VCF with deletion event calls'
	)
	parser.add_argument('wd', 
		action='store', 
		metavar='/path/to/output_directory/',
		help='Directory to put output in'
	)
	# Optional args
	parser.add_argument('-s', '--sample-name',
		action='store',
		metavar='sample_name',
		help='Prefix [default: filename of in_vcf]'
	)
	parser.add_argument('-t', '--threshold',
		action='store',
		metavar='n',
		type=int,
		help='Output pseudogenes must have at least this many qualifying deletion events. [default: 2]',
		default=2
	)
	parser.add_argument('--iAnnotateSV', 
		action='store', 
		metavar='iAnnotateSV.py',
		default='/home/shahr2/git/iAnnotateSV/iAnnotateSV/iAnnotateSV.py',
		help='''path to iAnnotateSV.py script
		[default: /home/shahr2/git/iAnnotateSV/iAnnotateSV/iAnnotateSV.py]'''
	)
	parser.add_argument('--genome',
		action='store',
		default='hg19',
		choices=['hg18', 'hg19', 'hg38'],
		help='Reference genome version to use during annotation [default: hg19]'
	)
	
	args = parser.parse_args(command.split()) if command else parser.parse_args()

	WD = args.wd
	IN_VCF = args.in_vcf
	PREFIX = args.sample_name or os.path.splitext(os.path.basename(IN_VCF))[0]
	iAnnotateSV = args.iAnnotateSV
	REF_GENOME = args.genome
	THRESHOLD = args.threshold
	scratch_dir = 'scratch/'

	os.chdir(WD)
	
	try:
		os.mkdir(scratch_dir)
	except OSError:
		pass # Assume directory exists already

	#############################################################################

	# 1: Filter VCF
	print 'Filtering vcf...'

	filtered_vcf = os.path.join(scratch_dir, PREFIX + '.filtered.vcf')

	with open(IN_VCF, 'rU') as vcf_in_fp, open(filtered_vcf, 'w') as filtered_vcf_fp:
		reader = vcf.Reader(vcf_in_fp)
		writer = vcf.Writer(filtered_vcf_fp, template=reader)

		for record in filter_vcf(reader):
			writer.write_record(record)


	# 2: Convert to tab delimited
	print 'Converting to tab delimited...'
	vcf2tab(filtered_vcf, scratch_dir, verbose=False)


	# 3: Annotate
	print 'Annotating...'

	annotate_command = 'python {iAnnotateSV} -r {genome} -ofp {prefix} -o {wd} -i {tab_file} -d 3000'.format(
					iAnnotateSV = iAnnotateSV,
					genome = REF_GENOME,
					prefix = PREFIX,
					wd = scratch_dir,
					tab_file = os.path.splitext(filtered_vcf)[0] + ".tab"
		)
	print annotate_command
	subprocess.check_call(shlex.split(annotate_command))


	# 4: find pseudogenes
	print 'Looking for pseudogenes...'

	# The file that iAnnotateSV created that we want
	annotations_file = os.path.join(scratch_dir, PREFIX + '_Annotated.txt')

	annotations_df = pd.read_table(annotations_file)
	pseudogenes = find_pseudogenes(annotations_df, THRESHOLD=THRESHOLD)
	pseudogenes.sort() # A to Z

	out_file = PREFIX + '.pseudogenes.txt'

	with open(out_file, 'w') as out_fp:
		for gene in pseudogenes:
			out_fp.write('%s\t%s\n' % (PREFIX, gene))


	print 'Done!'


if __name__ == '__main__':
	main()
