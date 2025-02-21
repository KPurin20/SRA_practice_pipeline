#!/usr/bin/env python
# coding: utf-8

# In[ ]:


import subprocess
import os
from Bio import SeqIO
import statistics
from multiprocessing import Pool
import matplotlib.pyplot as plt
from datetime import datetime
import logging

logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')

def set_output_directory(base_directory, sra_id):
    """
    Set the output directory for the SRA ID.
    """
    directories = {
        "fastq": os.path.join(base_directory, sra_id, 'fastq'),
        'reports': os.path.join(base_directory, sra_id, 'reports'),
        'processed': os.path.join(base_directory, sra_id, 'processed')
    }
    for directory_path in directories.values():
        os.makedirs(directory_path, exist_ok=True)
    return directories

def download_and_convert_sra(sra_id, output_directory):
    """
    Download and convert SRA to FASTQ.
    """
    try:
        logging.info(f"Downloading {sra_id}")
        subprocess.run(['prefetch', sra_id], check=True)

        logging.info(f"Converting {sra_id} to FASTQ")
        subprocess.run(['fastq-dump', '--split-files', sra_id, '-O', output_directory], check=True)

        return [f for f in os.listdir(output_directory) if f.endswith('.fastq')]
    
    except subprocess.CalledProcessError as e:
        logging.error(f"Error processing {sra_id}: {str(e)}")
        return None

def quality_control(records, min_length=50, min_quality=20):
    """
    Filter sequences based on length and quality.
    """
    filtered = []
    for record in records:
        qual = record.letter_annotations['phred_quality']
        if (len(record.seq) >= min_length and
            statistics.mean(qual) >= min_quality):
            filtered.append(record)
    return filtered

def trim_adapters(records, adapter_seq="AGATCGGAAGAGC"):
    """
    Trim adapters from sequences.
    """
    trimmed = []
    for record in records:
        seq_str = str(record.seq)
        if adapter_seq in seq_str:
            adapter_position = seq_str.find(adapter_seq)
            seq_str = seq_str[:adapter_position]

        quality = record.letter_annotations['phred_quality']
        start = 0
        end = len(quality)

        while start < end and quality[start] < min_quality:
            start += 1
        while end > start and quality[end - 1] < min_quality:
            end -= 1

        trimmed_record = record[start:end]
        if len(trimmed_record) > 0:
            trimmed.append(trimmed_record)
    return trimmed

def find_motif(records, motif_sequence):
    """
    Find motif in sequences.
    """
    motif_positions = {}
    for i, record in enumerate(records):
        position = []
        pos = 0
        while True:
            pos = str(record.seq).find(motif_sequence, pos)
            if pos == -1:
                break
            position.append(pos)
            pos += 1
        if position:
            motif_positions[record.id] = position
    return motif_positions

def generate_sequence_stats(records):
    """
    Generate sequence statistics.
    """
    stats = {
        'total_sequences': len(records),
        'length_stats': {'min': float('inf'), 'max': 0, 'mean': 0},
        'gc_content': [],
        'quality_scores': [],
        'base_composition': {base: 0 for base in 'ATGCN'}
    }
    
    lengths = []
    for record in records:
        # Length stats
        length = len(record.seq)
        lengths.append(length)
        stats['length_stats']['min'] = min(stats['length_stats']['min'], length)
        stats['length_stats']['max'] = max(stats['length_stats']['max'], length)

        # GC content
        gc_content = (record.seq.count('G') + record.seq.count('C')) / length
        stats['gc_content'].append(gc_content)

        # Quality scores
        stats['quality_scores'].append(record.letter_annotations['phred_quality'])

        # Base composition
        for base in 'ATGCN':
            stats['base_composition'][base] += record.seq.count(base)

    stats['length_stats']['mean'] = statistics.mean(lengths)
    stats['mean_gc_content'] = statistics.mean(stats['gc_content'])
    stats['mean_quality_scores'] = statistics.mean(stats['quality_scores'])

    return stats

def generate_quality_plots(records, output_directory, sra_id):
    """
    Generate quality plots.
    """
    plt.figure(figsize=(10, 6))

    # Position-based quality scores
    positions = range(len(records[0].letter_annotations['phred_quality']))
    qualities = [[] for _ in positions]

    for record in records:
        for i, quality in enumerate(record.letter_annotations['phred_quality']):
            qualities[i].append(quality)

    avg_qualities = [statistics.mean(pos_qualities) for pos_qualities in qualities]

    plt.plot(positions, avg_qualities)
    plt.xlabel('Position in read')
    plt.ylabel('Average Quality Score')
    plt.title(f'Quality Score Distribution for {sra_id}')
    plt.savefig(os.path.join(output_directory, f'{sra_id}_quality_distribution.png'))
    plt.close()

def process_single_sra(sra_id, base_output_directory):
    """
    Process a single SRA file.
    """
    try:
        directory = set_output_directory(base_output_directory, sra_id)
        # Download and convert SRA
        fastq_files = download_and_convert_sra(sra_id, directory['fastq'])
        if not fastq_files:
            return False
        for fastq_file in fastq_files:
            file_path = os.path.join(directory['fastq'], fastq_file)
            logging.info(f"Processing {fastq_file}")
            # Read sequences
            records = list(SeqIO.parse(file_path, 'fastq'))
            # Initial stats
            initial_stats = generate_sequence_stats(records)
            # Quality control
            filtered_records = quality_control(records)
            # Adapter trimming
            trimmed_records = trim_adapters(filtered_records)
            # Motif finding
            motif_positions = find_motif(trimmed_records, 'AGATCGGAAGAGC')
            # Generate stats
            final_stats = generate_sequence_stats(trimmed_records)
            # Generate plots
            generate_quality_plots(trimmed_records, directory['reports'], sra_id)
            # Save processed sequences
            output_file = os.path.join(directory['processed'], f'{sra_id}_processed.fastq')
            SeqIO.write(trimmed_records, output_file, 'fastq')

            # Save stats
            with open(os.path.join(directory['reports'], f'{sra_id}_stats.txt'), 'w') as f:
                f.write(f"Processing report for {sra_id}\n")
                f.write(f"Date: {datetime.now()}\n")
                f.write(f"Initial Stats:\n")
                f.write(str(initial_stats))
                f.write(f"Final Stats:\n")
                f.write(str(final_stats))
                f.write(f"Motif Positions:\n")
                f.write(str(motif_positions))

        return True
    except Exception as e:
        logging.error(f"Error processing {sra_id}: {str(e)}")
        return False

def process_multiple_sra(sra_list, base_output_directory, num_processes=4):
    """
    Process multiple SRA files.
    """
    with Pool(processes=num_processes) as pool:
        from functools import partial
        process_func = partial(process_single_sra, base_output_directory=base_output_directory)
        results = pool.map(process_func, sra_list)
    return results

# Test case
if __name__ == "__main__":
    sra_numbers = ['SRR031507', 'SRR031508', 'SRR031509']
    output_directory = './sra_processed_data'

    results = process_multiple_sra(sra_numbers, output_directory)

    for sra, success in zip(sra_numbers, results):
        status = "completed" if success else "failed"
        logging.info(f"Processing of {sra} {status}")

