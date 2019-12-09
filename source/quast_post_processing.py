import os 
from Bio import SeqIO
import pandas as pd 
import matplotlib.pyplot as plt

contigs_list = []
contigs_under_1000_list =[]

metadata = pd.read_csv('rmash_files/Validated_Species.csv', converters={'Record ID': lambda x: str(x)})

output_names = []
for index, row in metadata.iterrows():
	record = row['Record ID']
	sample = row['Sample']
	
	quast_report = pd.read_csv ("quast_files/quast_"+record+"/report.tsv", sep = '\t', index_col = 0).T
	total_contigs = int(quast_report['# contigs'])
	contigs_under_1000 = int(quast_report['# contigs (>= 0 bp)'])
	contigs_list.append(total_contigs)
	contigs_under_1000_list.append(contigs_under_1000)

	location = "refseq/bacteria/"+sample[0:15]+"/"+sample+".fna"
	output_name = sample.replace('.', '_')
	output_names.append(output_name)
	output_location = 'Approved_Sequences/'+output_name+'.fna'
	
	if contigs_under_1000 > 0:
		over_500_bp = []
		for record in SeqIO.parse(location,'fasta'):
			if len(record.seq) > 500:
				over_500_bp.append(record)
		SeqIO.write(over_500_bp,output_location,'fasta')
	else:
		all_records = []
		for record in SeqIO.parse(location,'fasta'):
			all_records.append(record)
		SeqIO.write(all_records, output_location,'fasta')

metadata['Sample'] = output_names
metadata['Number of Contigs']= contigs_list
metadata.to_csv('Metadata.csv', index = False)

scatterplot = plt.scatter(contigs_list, contigs_under_1000_list, s=10, c='k')
plt.xlabel ('Number of Contigs')
plt.ylabel('Number of Contigs <1000bp')
 
plt.title('Data Quality Visualization')
plt.savefig("Data Quality Visualization.png")