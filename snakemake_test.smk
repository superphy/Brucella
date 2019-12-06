import glob 
import os
cwd = os.getcwd()

rule all:
	input:
		'quast_files/quast_files.txt'

#Aquiring the data from the ncbi database
rule ncbi_data_retrieval:
	output:
		"refseq/"
	shell:
		'ncbi-genome-download --parallel 55 --genus "brucella" bacteria --format fasta'

mds = glob.glob('refseq/bacteria/*/MD5SUMS')
unzips = glob.glob('refseq/bacteria/*/*.gz')

#Removes extra files downloaded from NCBI
rule rmv_MD5SUMS:
	input:
		refseq = "refseq/",
		md = '{md}'
	output:
		'flags/md5sums/{md}.txt'
	run:
		shell('rm {input.md}')
		shell('touch {output}')

#Unzips NCBI sequences
rule unzip:
	input:
		refseq = 'refseq/',
		unzip = '{unzip}'
	output:
		'flags/unzip/{unzip}.txt'
	run:
		shell('gunzip {input.unzip}')
		shell('touch {output}')	

#Obtains a list of file names for refseq masher
rule r_mash_files:
	input:
		md = expand('flags/md5sums/{md}.txt', md=mds),
		unzip = expand('flags/unzip/{unzip}.txt', unzip = unzips)
	output:
		'rmash_files/rmash_files.txt'
	run:
		files = glob.glob('refseq/bacteria/*/*')
		files = " ".join(files)
		r_mash_file = open('rmash_files/rmash_files.txt', 'w')
		r_mash_file.write(files)
		r_mash_file.close()

#Runs refseq masher 
rule r_mash:
	conda:
		'envs/mash.yaml'
	input:
		'rmash_files/rmash_files.txt'
	output:
	 	'rmash_files/rmash_output.csv' 
	params:
		files = open('rmash_files/rmash_files.txt').read()
	shell:
		'refseq_masher matches {params.files} -o rmash_files/rmash_output.csv '+cwd+' --output-type csv -n 5'

rule r_mash_filtering:
	input:
		rmash_out = 'rmash_files/rmash_output.csv',
		filtering = 'source/r_mash_filtering.py'
	output:
		'rmash_files/rmash_filtered_output.csv'
	run:
		shell('python {input.filtering}')	

rule species_validation:
	input:
		filtered_out = 'rmash_files/rmash_filtered_output.csv',
		species_val = 'source/species_validation.py'
	output:
		'rmash_files/Validated_Species.csv'
	run:
		shell('python {input.species_val}')

rule quast_input:
	input:
		'rmash_files/Validated_Species.csv'
	output:
		'quast_files/quast_files.txt'
	run:
		files = glob.glob('refseq/bacteria/*/*')
		quast_file = open('quast_files/quast_files.txt', 'w')
		quast_file.write(quast_file)
'''
rule quast:	
	input:
		'quast_files/quast_files.txt'
	output:
		fnke
	run:

'''
