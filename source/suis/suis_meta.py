import pandas as pd 
from Bio import SeqIO
import os

metadata = pd.read_csv('Metadata.csv', index_col = 0)
suis_meta = metadata[metadata['Strain']=='Brucella suis']
suis_files = ['Approved_Sequences/'+sample+".fna" for sample in suis_meta.index.values]

def new_contig(fasta_path):
	for record in SeqIO.parse(fasta_path, "fasta"):
		contig_seq = record.seq
		contig_seq = contig_seq._get_seq_str_and_check_alphabet(contig_seq)
		contig_header = record.id
		file_name = fasta_path.split('GCF_')[-1].split('_')[0]
		contig_header = ">{}_{}".format(file_name, contig_header)
		with open("suis_biovar_files/suis_master.fna",'a') as master:
			master.write(contig_header)
			master.write("\n")
			master.write(contig_seq)
			master.write("\n")
	return 0

for fasta in suis_files:
	new_contig(fasta)
