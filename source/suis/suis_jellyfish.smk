import pandas as pd 
suis_biovar = pd.read_csv('Suis_Biovar.csv', index_col = 0)
suis_biovar = suis_biovar[suis_biovar['Biovar Classification'] != 'No Matched Biovar']
samples = list(suis_biovar.index.values)
#samples = [sample+'.fna' for sample in samples]
kmer_size = 31

rule all:
	input:
		all_dump = 'j_all_suis.fna',
		comp = 'Suis_Ranks.csv'

#Performs the jellyfish count opperation on all_brucella.fna
rule all_jellyfish_count: 
	input: 
		'suis_master.fna'
	output:
		temp('j_all_suis.jf')
	run:
		shell('jellyfish count -C -m {kmer_size} -s 100M -t 2 {input} -o {output} -L 190 -U 501 --out-counter-len 1')

#Performs the jellyfish dump opperation on all_brucella.fna
rule all_jellyfish_dump: 
	input: 
		'j_all_suis.jf'
	output:
		'j_all_suis.fna'
	run:
		shell('jellyfish dump {input} > {output}')

#Performs the jellyfish count opperation on everything in the dataset
rule count:
	input:
		"../../Approved_Sequences/{samples}.fna"
	output:
		temp("suis_jellyfish_output/{samples}.jf")
	shell:
		"jellyfish count -C -m {kmer_size} -s 100M -t 2 {input} -o {output} --out-counter-len 1"

#Performs the jellyfish dump opperation on everything in the dataset
rule dump:
	input:
		"suis_jellyfish_output/{samples}.jf"
	output:
		"suis_jellyfish_output/{samples}.fna", 
	shell:
		"jellyfish dump {input} > {output}"

#Makes count/dump actually work
rule jellyfish_complete:
	input:
		expand("suis_jellyfish_output/{samples}.fna", samples = samples)
	output:
		'flags/suis_jellyfish_flag.txt'
	run:
		shell('touch flags/suis_jellyfish_flag.txt')


#Creates kmer_counts.csv 
rule kmer_dataframe:
	input:
		jf = 'flags/suis_jellyfish_flag.txt',
		c_df = 'suis_count_to_df.py'
	output:
		kc = 'suis_kmer_counts.csv'
	run:
		shell('python {input.c_df}')

#Analyzes the occurances of a specific kmer in each species
rule ranks:
	input:
		ra = 'suis_ranks.py', 
		kc = 'suis_kmer_counts.csv'
	output:
		'Suis_Ranks.csv'
	run:
		shell('python {input.ra}')