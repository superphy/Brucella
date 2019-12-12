import os, sys
import pandas as pd
from Bio import SeqIO
import numpy as np
from concurrent.futures import ProcessPoolExecutor
from multiprocessing import cpu_count

if __name__ == '__main__':
	i = 0
	genome_counter = 0
	times_printed = 0
	num_start = cpu_count()-15
	num_stop = 0 
	total = 0
	kmer_dict = {}

	genome_list = os.listdir("jellyfish/output/")
	total = len(genome_list)

	# Creating a dict of all the kmers in the master file
	for record in SeqIO.parse("jellyfish/j_all_brucella.fna", "fasta"):
		kmer_seq_all = record.seq
		kmer_seq_all = kmer_seq_all._get_seq_str_and_check_alphabet(kmer_seq_all)
		kmer_dict[kmer_seq_all]=i
		i+=1

	df = pd.DataFrame(index = kmer_dict, columns = genome_list, data = np.zeros((len(kmer_dict), len(genome_list)), dtype = "uint8"))

	# For each genome a df column is created to representing the amount of times each kmer is represented 
	def kmer_check(genome):
		location = "jellyfish/output/"+genome
		col_list = np.zeros(len(kmer_dict), dtype = "uint8")
		new_fasta = []
		for record in SeqIO.parse(location, "fasta"):
			kmer_seq = record.seq
			kmer_seq = kmer_seq._get_seq_str_and_check_alphabet(kmer_seq)
			if kmer_seq in kmer_dict:
				col_list[kmer_dict[kmer_seq]] = np.uint8(1)
				new_fasta.append(record)
		SeqIO.write(new_fasta, "filtered_mer_counts/"+genome,"fasta")
		return(genome, col_list)

	# Tracks the progress of the script
	def progress():
		sys.stdout.write('\r')
		sys.stdout.write("Loading Genomes: {} started, {} finished, {} total".format(num_start,num_stop,total))
		sys.stdout.flush()
		if(num_stop==total):
			print("\nAll Genomes Loaded!\n")

	progress()
	# Runs kmer_check as concurent features
	with ProcessPoolExecutor(cpu_count()-15) as ppe:
		for genome, column in ppe.map(kmer_check, genome_list):
			df[genome] = column
			num_stop+=1
			if(num_start<total):
				num_start+=1
			progress()

	df.to_csv("kmer_counts.csv")
