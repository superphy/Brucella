import glob 
import os
import pandas as pd
from multiprocessing import cpu_count
kmer_size = 31
cpu = cpu_count() - 6
cwd = os.getcwd()

rule all:
	input:
		kover_in = 'Kover_Data/Kmer_Matrix.tsv', 
		tree = 'phylogenetic_tree.pdf', 
		rank_hist = 'Visualizations/Rank_Histograms/', 
		suis_k_loc = 'suis_biovar_files/Suis_Kmer_Locations.csv', 
		Top_200 = 'Top_200_bp.csv'

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
		kfiles = [cwd+'/Approved_Sequences/'+file+"	"+file[4:13]+'\n' for file in files if file[0] != '.']
		with open('ksnp_files/kSNP3_input.txt', 'w') as ksnp_input_file:
			for row in kfiles:
				ksnp_input_file.write(row)
		ksnp_input_file.close()		

rule ksnp:
	input:
		'ksnp_files/kSNP3_input.txt'
	output:
		'ksnp_files/Logfile.txt'
	run:
		#shell("PATH=$PATH:~/Desktop/Fall_2019/Brucella_Snakemake_Test_Evn/kSNP3.1_Linux_package/kSNP3")
		shell('kSNP3 -in {input} -outdir ksnp_files/ -k 31 -CPU {cpu} | tee {output}')

rule phylogenetic_tree:
	input:
		ksnp_flag = 'ksnp_files/Logfile.txt', 
		tree_manipulation = 'source/tree_manipulation.py'
	output:
		'phylogenetic_tree.pdf'
	run:
		shell('python {input.tree_manipulation}')

rule suis_master_fasta:
	input:
		suis_meta = 'source/suis/suis_meta.py',
		approved_seq = 'Approved_Sequences/'
	output:
		'suis_biovar_files/suis_master.fna'
	run:
		shell('python {input.suis_meta}')

rule suis_blast_db:
	input:
		master = 'suis_biovar_files/suis_master.fna'
	output:
		blast_db = 'suis_biovar_files/blast/suis_db.ndb'
	run:
		shell('makeblastdb -in {input} -parse_seqids -blastdb_version 5 -title "Brucella Suis" -dbtype nucl -out suis_biovar_files/blast/suis_db')

rule suis_blast_search:
	input:
		query_file = 'source/suis/query_file.fna',
		blast_db = 'suis_biovar_files/blast/suis_db.ndb'
	output:
		serch_out = 'suis_biovar_files/blast/blast_search_output.tsv'
	run:
		shell('blastn -db suis_biovar_files/blast/suis_db -query {input.query_file} -dust no -word_size 7 -evalue 100 -outfmt 6 -out {output}')

rule suis_blast_output:
	input:
		blast_search_output = 'suis_biovar_files/blast/blast_search_output.tsv',
		bo = 'source/suis/blast_output.py'
	output:
		'suis_biovar_files/Blast_Primer_Matching.csv'
	run:
		shell('python {input.bo}')	

rule suis_biovar_assignment:
	input:
		ba = 'source/suis/biovar_assignment.py',
		bpm = 'suis_biovar_files/Blast_Primer_Matching.csv'

	output:
		'suis_biovar_files/Suis_Biovar.csv'
	run:
		shell('python {input.ba}')

#Generates a fasta file that is the sum of all fasta files in the dataset
rule jellyfish_all_brucella:
	input:
		'Approved_Sequences/'
	output:
		'jellyfish/all_brucella.fna'
	run:
		shell('cat {input}*.fna > {output}')
#Performs the jellyfish count opperation on all_brucella.fna
rule all_jellyfish_count: 
	input: 
		'jellyfish/all_brucella.fna'
	output:
		temp('jellyfish/j_all_brucella.jf')
	run:
		shell('jellyfish count -C -m {kmer_size} -s 100M -t {cpu} {input} -o {output} -L 190 -U 501 --out-counter-len 1')

#Performs the jellyfish dump opperation on all_brucella.fna
rule all_jellyfish_dump: 
	input: 
		'jellyfish/j_all_brucella.jf'
	output:
		'jellyfish/j_all_brucella.fna'
	run:
		shell('jellyfish dump {input} > {output}')

rule suis_all_jellyfish_count:
	input:
		'suis_biovar_files/suis_master.fna'
	output:
		temp('suis_biovar_files/j_all_suis.jf')
	run:
		shell('jellyfish count -C -m {kmer_size} -s 100M -t {cpu} {input} -o {output} -L 190 -U 501 --out-counter-len 1')

rule suis_all_jellyfish_dump:
	input:
		'suis_biovar_files/j_all_suis.jf'
	output:
		'suis_biovar_files/j_all_suis.fna'
	run:
		shell('jellyfish dump {input} > {output}')

#Performs the jellyfish count opperation on everything in the dataset
rule jellyfish:
	input:
		"Approved_Sequences/"
	output:
		flag = 'jellyfish/jellyfish_complete.txt', 
		output_direct = directory('jellyfish/output/')
	run:
		samples = glob.glob('Approved_Sequences/*.fna')
		for sample in samples:	
			shell('jellyfish count -C -m {kmer_size} -s 100M -t 2 '+sample+' -o jellyfish/output/'+sample[19:len(sample)-4]+'.jf --out-counter-len 1') # runs the jellyfish count opperation 
			shell('jellyfish dump jellyfish/output/'+sample[19:len(sample)-4]+'.jf > jellyfish/output/'+sample[19:len(sample)-4]+'.fna') # runs the jellyfish dump opperation 
			shell('rm jellyfish/output/'+sample[19:len(sample)-4]+'.jf')
		shell('touch {output}')

rule suis_jellyfish:
	input:
		approved_seq = "Approved_Sequences/", 
		suis_biovar = 'suis_biovar_files/Suis_Biovar.csv'
	output:
		flag = "suis_biovar_files/jellyfish_complete.txt", 
		output_direct = directory('suis_biovar_files/jellyfish/')
	run:
		suis_biovar = pd.read_csv('suis_biovar_files/Suis_Biovar.csv', index_col = 0)
		suis_biovar = suis_biovar[suis_biovar['Biovar Classification'] != 'No Matched Biovar']
		samples = suis_biovar.index.values
		for sample in samples:	
			shell('jellyfish count -C -m {kmer_size} -s 100M -t 2 Approved_Sequences/'+sample+'.fna -o suis_biovar_files/jellyfish/'+sample+'.jf --out-counter-len 1') # runs the jellyfish count opperation 
			shell('jellyfish dump suis_biovar_files/jellyfish/'+sample+'.jf > suis_biovar_files/jellyfish/'+sample+'.fna') # runs the jellyfish dump opperation 
			shell('rm suis_biovar_files/jellyfish/'+sample+'.jf')
		shell('touch {output}')

rule kmer_counts:
	input:
		flag = 'jellyfish/jellyfish_complete.txt',
		k_counts = 'source/kmer_counts.py'
	output:
		'Kmer_Counts.csv'
	run:
		shell('python {input.k_counts}')

rule suis_kmer_counts:
	input:
		flag = 'suis_biovar_files/jellyfish_complete.txt',
		k_counts = 'source/suis/suis_kmer_counts.py'
	output:
		'suis_biovar_files/Suis_Kmer_Counts.csv'
	run:
		shell('python {input.k_counts}')	

rule species_occ:
	input:
		'Metadata.csv'
	output:
		'Species_Occurrence.csv'
	run:	
		metadata = pd.read_csv('Metadata.csv', index_col = 0 )
		species = list(metadata['Strain'].unique())
		species_occ = {}
		for sp in species:
			species_occ[sp] = 0
		for i, row in metadata.iterrows():
			species_occ[row['Strain']]+=1
		df = pd.DataFrame()
		for i in species_occ.keys():
			df[i] = pd.Series(species_occ[i])
		df = df.T.rename(columns = {0:'Species Occurrence'})
		df.to_csv('Species_Occurrence.csv')

rule ranks:
	input:
		species_occ = 'Species_Occurrence.csv', 
		k_counts = 'Kmer_Counts.csv',
		ranks = 'source/ranks.py'
	output:
		'Ranks.csv'
	run:
		shell('python {input.ranks}')

rule suis_ranks:
	input:
		suis_biovar = 'suis_biovar_files/Suis_Biovar.csv',
		sk_counts = 'suis_biovar_files/Suis_Kmer_Counts.csv',
		s_ranks = 'source/suis/suis_ranks.py'
	output:
		'suis_biovar_files/Suis_Ranks.csv'
	run:	
		shell('python {input.s_ranks}')

rule rank_histograms:
	input:
		ranks = 'Ranks.csv',
		sp_occ = 'Species_Occurrence.csv',
		s_ranks = 'suis_biovar_files/Suis_Ranks.csv',
		s_biovar = 'suis_biovar_files/Suis_Biovar.csv',
		rank_hist = 'source/rank_histograms.py',
		s_rank_hist = 'source/suis/suis_rank_histograms.py'
	output:
		directory('Visualizations/Rank_Histograms/')
	run:
		shell('python {input.rank_hist}')
		shell('python {input.s_rank_hist}')

rule chrom_mapping:
	input:
		metadata = 'Metadata.csv',
		ranks = 'Ranks.csv',
		sp_occ = 'Species_Occurrence.csv',
		chrom_mapping = 'source/chrom_mapping.py'
	output:
		visualizations = directory('Visualizations/Kmer_Locations/'), 
		ref_genomes = 'Refrence_Genomes.csv',
		k_loc = 'Kmer_Locations.csv'
	run:
		shell('python {input.chrom_mapping}')

rule suis_chrom_mapping:
	input:
		metadata = 'Metadata.csv',
		ranks = 'suis_biovar_files/Suis_Ranks.csv',
		s_biovar = 'suis_biovar_files/Suis_Biovar.csv',
		chrom_mapping = 'source/suis/suis_chrom_mapping.py'
	output:
		ref_genomes = 'suis_biovar_files/Suis_Refrence_Genomes.csv',
		suis_k_loc = 'suis_biovar_files/Suis_Kmer_Locations.csv'
	run:
		shell('python {input.chrom_mapping}')

rule top_200_bp:
	input:
		ranks = 'Ranks.csv',
		ref_genomes = 'Refrence_Genomes.csv',
		k_loc = 'Kmer_Locations.csv',
		strain_occ = 'Species_Occurrence.csv',
		top_200 = 'source/top_200.py'
	output:
		vis = directory('Visualizations/Top_200_bp/'), 
		Top_200 = 'Top_200_bp.csv'
	run:
		shell('python {input.top_200}')

#Not implemented on suis files - ran out of time...Shouldnt be too dificult to map the existing function to the suis data
rule suis_top_200_bp:
	input:
		refrence_genomes = 'suis_biovar_files/Suis_Refrence_Genomes.csv'
		k_loc = 'suis_biovar_files/Suis_Kmer_Locations.csv'
		suis_biovar ='suis_biovar_files/Suis_Biovar.csv'
		ranks_df = 'suis_biovar_files/Suis_Ranks.csv'
		top_200 = 'source/suis/suis_top_200.py'
	output:
		'suis_biovar_files/Suis_Top_200_bp.csv'
	run:
		shell('python {input.top_200}')


rule kover_input:
	input:
		metadata = 'Metadata.csv',
		k_counts = 'Kmer_Counts.csv',
		strain_occ = 'Species_Occurrence.csv',
		k_input = 'source/kover_inputs.py'
	output:
		'Kover_Data/Kmer_Matrix.tsv'
	run:
		shell('python {input.k_input}')


#If you would like to run Kover on this dataset, run the seperate kover smk. 