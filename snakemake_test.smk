import glob 
import os
import pandas as pd
from multiprocessing import cpu_count
kmer_size = 31
cpu = cpu_count() - 6
cwd = os.getcwd()

rule all:
	input:
		'ksnp_files/kSNP3_Output/Logfile.txt'
'''
#Aquiring the data from the ncbi database
rule ncbi_data_retrieval:
	output:
		"refseq/"
	shell:
		'ncbi-genome-download --parallel 55 --genus "brucella" bacteria --format fasta'
'''
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
	shell:
		'refseq_masher matches $(cat {input}) -o rmash_files/rmash_output.csv '+cwd+' --output-type csv -n 5'

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

rule quast:
	input:
		'rmash_files/Validated_Species.csv'
	output:
		'quast_files/quast_files.txt'
	run:
		Validated_Species = pd.read_csv('rmash_files/Validated_Species.csv')
		samples = list(Validated_Species['Sample'])
		for sample in samples:
			shell('python quast-5.0.2/quast.py --silent --fast -o quast_files/quast_'+sample[4:13]+' refseq/bacteria/'+sample[0:15]+'/'+sample+'.fna')
		shell('touch {output}')

rule quast_pp:
	input:
		quast_complete = 'quast_files/quast_files.txt', 
		species_val = 'rmash_files/Validated_Species.csv',
		quast_pp = 'source/quast_post_processing.py'
	output:
		metadata = 'Metadata.csv',
		data_vis = 'Data Quality Visualization.png',
		approved_seq = directory('Approved_Sequences/')
	run:
		shell('python {input.quast_pp}')

rule ksnp_input:
	input:
		'Metadata.csv'
	output:
		'ksnp_files/kSNP3_input.txt'
	run:
		files = os.listdir('Approved_Sequences/')
		kfiles = [cwd+'/Approved_Sequences/'+file+"	"+file[4:13]+'\n' for file in files]
		with open('ksnp_files/kSNP3_input.txt', 'w') as ksnp_input_file:
			for row in kfiles:
				ksnp_input_file.write(row)
		ksnp_input_file.close()		

rule ksnp:
	input:
		'ksnp_files/kSNP3_input.txt'
	output:
		'ksnp_files/kSNP3_Output/Logfile.txt'
	run:
		#shell("PATH=$PATH:~/Desktop/Fall_2019/Brucella/kSNP3.1_Linux_package/kSNP3")
		shell('kSNP3 -in {input} -outdir ksnp_files/kSNP3_Output/ -k 31 -CPU {cpu} | tee {output}')