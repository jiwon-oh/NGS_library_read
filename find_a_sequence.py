import os
from Bio import SeqIO
from pathlib import Path
import csv

#read sequences in fastq file 
work_dir = Path('/Users/jiwonoh/Desktop/Python').absolute().parent
file_path = work_dir / '/Users/jiwonoh/Desktop/Python/Ben1_merged.fastqjoin'

with open('Ben1.csv', 'w', newline='') as csvfile:
        writer = csv.writer(csvfile)
        for seq_record in SeqIO.parse(file_path, "fastq"):
             print(seq_record.id)
             print(repr(seq_record.seq))
             print(len(seq_record))
             writer.writerow([seq_record.id, len(seq_record.seq), str(seq_record.seq)])

#split the sequences based on barcode
import subprocess
import shutil 
output_dir = '/Users/jiwonoh/Desktop/output/' 
input_dir = '/Users/jiwonoh/Desktop/' 

os.chdir(input_dir)
fastx_barcode_splitter = shutil.which('fastx_barcode_splitter.pl')

command = [
    fastx_barcode_splitter,
    '--bcfile', 'Azprimers.txt',
    '--bol',
    '--prefix', output_dir,
    '--suffix', '.fastq',
    '--mismatches', '0',
    '<', 'Ben1_merged.fastqjoin'
]

#Translate nucleotide into amino acid sequence at the start codon
from pathlib import Path
from Bio import SeqIO
from Bio.Seq import Seq


work_dir = Path('/Users/jiwonoh/Desktop/output/')
Ben1_PD = work_dir / 'BC1.fastq'
Ben1_FT = work_dir / 'unmatched.fastq'

def translate_sequences(file_paths, csv_writer):
    for file_path in file_paths:
        for record in SeqIO.parse(file_path, "fastq"):
            start_codon = 'TACATG'
            start_index = record.seq.find(start_codon)
            trimmed_seq = None  
            trimmed_seq = record.seq[start_index:len(record.seq) - (len(record.seq) % 3)]
            translated_seq = trimmed_seq.translate()
            print(translated_seq)
            csv_writer.writerow([record.id, len(str(translated_seq)), str(translated_seq)])

#save the translated sequences into csvfile
with open('Ben1_PD.csv', 'w', newline='') as csvfile:
    writer = csv.writer(csvfile)
    translate_sequences([Ben1_PD], writer)

with open('Ben1_FT.csv', 'w', newline='') as csvfile:
    writer = csv.writer(csvfile)
    translate_sequences([Ben1_FT], writer)

from collections import Counter

PD_file_path = '/Users/jiwonoh/Desktop/Ben1_PD.csv'
FT_file_path = '/Users/jiwonoh/Desktop/Ben1_FT.csv'

matching_sequences = []
target_sequence = 'PPLYGDNLDQHFRLLAQKQSLPYLQAANLLLQAQLPPKPPAWAWAEGWTR'

with open(PD_file_path, 'r') as csvfile:
    csv_reader = csv.reader(csvfile)

    for row in csv_reader:
        sequence = row[2]

        if target_sequence in sequence:  
            matching_sequences.append(sequence)  

frequency_data = Counter(matching_sequences)

for seq in matching_sequences:
    print (seq)
for count in frequency_data.items():  
    print(count)
with open('matching_sequences.csv', 'w', newline='') as csvfile:
    writer = csv.writer(csvfile)   
    writer.writerow(['Sequence', 'Length', 'Frequency']) 

    for seq, count in frequency_data.items():  
        writer.writerow([seq, len(seq), count])  


matching_sequences1 = []
with open(FT_file_path, 'r') as csvfile:
    csv_reader = csv.reader(csvfile)

    for row in csv_reader:
        sequence = row[2]

        if target_sequence in sequence:  
            matching_sequences1.append(sequence)  

frequency_data = Counter(matching_sequences1)

for seq in matching_sequences1:
    print (seq)
for count in frequency_data.items():  
    print(count)
with open('matching_sequences1.csv', 'w', newline='') as csvfile:
    writer = csv.writer(csvfile)   
    writer.writerow(['Sequence', 'Length', 'Frequency']) 

    for seq, count in frequency_data.items():  
        writer.writerow([seq, len(seq), count])  