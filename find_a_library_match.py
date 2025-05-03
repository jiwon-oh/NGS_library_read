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

#Split the sequences based on barcode
import subprocess
import shutil 
output_dir = '/Users/jiwonoh/Desktop/NGS_python/output' 
input_dir = '/Users/jiwonoh/Desktop/NGS_python/Ben1_Screen_test' 

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

#Translate nucleotides into amino acid sequence at the start codon and stop at the stop codon
from pathlib import Path
from Bio import SeqIO
from Bio.Seq import Seq


work_dir = Path('/Users/jiwonoh/Desktop/NGS_python/output')
Ben1_PD = work_dir / 'BC1.fastq'
Ben1_FT = work_dir / 'unmatched.fastq'

def translate_sequences(file_paths, csv_writer):
    for file_path in file_paths:
        for record in SeqIO.parse(file_path, "fastq"):
            start_codon = 'TACATG'
            stop_codon = 'TGTACAAA'
            start_index = record.seq.find(start_codon)
            stop_index = record.seq.find(stop_codon)
            trimmed_seq = None  
            trimmed_seq = record.seq[start_index:stop_index - (len(record.seq) % 3)]
            translated_seq = trimmed_seq.translate()
            print(translated_seq[2:])
            csv_writer.writerow([record.id, len(str(translated_seq)), str(translated_seq[2:])])

#save the translated sequences into csvfile
with open('Ben1_PD.csv', 'w', newline='') as csvfile:
    writer = csv.writer(csvfile)
    translate_sequences([Ben1_PD], writer)

with open('Ben1_FT.csv', 'w', newline='') as csvfile:
    writer = csv.writer(csvfile)
    translate_sequences([Ben1_FT], writer)

#Find sequences in a library and count the frequency 
#target_sequence = 'PPLYGDNLDQHFRLLAQKQSLPYLQAANLLLQAQLPPKPPAWAWAEGWTR'
from collections import Counter
import pandas as pd
frequency_data = Counter()

PD_file_path = '/Users/jiwonoh/Desktop/NGS_python/Ben1_Screen_test/Ben1_PD.csv'
FT_file_path = '/Users/jiwonoh/Desktop/NGS_python/Ben1_Screen_test/Ben1_FT.csv'
target_library = '/Users/jiwonoh/Desktop/NGS_python/Ben1_Screen_test/out2_WES_pt1.csv'

target_lib_df = pd.read_csv(r'/Users/jiwonoh/Desktop/NGS_python/Ben1_Screen_test/out2_WES_pt1.csv', usecols=[4])  
pd_df = pd.read_csv(r'/Users/jiwonoh/Desktop/NGS_python/Ben1_Screen_test/Ben1_PD.csv')  
ft_df = pd.read_csv(r'/Users/jiwonoh/Desktop/NGS_python/Ben1_Screen_test/Ben1_FT.csv')
target_lib_sequences = target_lib_df.iloc[:, 0].dropna().tolist()  
pd_sequences = pd_df.iloc[:, 2].dropna().tolist()  
ft_sequences = ft_df.iloc[:, 2].dropna().tolist()     

matching_sequences = []
for pd_seq in pd_sequences:
    pd_seq_str = str(pd_seq) 
    for target_seq in target_lib_sequences:
        target_seq_str = str(target_seq)  
        if pd_seq_str in target_seq_str:  
            matching_sequences.append(pd_seq_str)
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
for ft_seq in ft_sequences:
    ft_seq_str = str(ft_seq) 
    for target_seq in target_lib_sequences:
        target_seq_str = str(target_seq)  
        if ft_seq_str in target_seq_str:  
            matching_sequences1.append(ft_seq_str)
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

