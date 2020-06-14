#! usr/bin/python3.6

# 2018-01-03, Carolina Monzo

# Modified for exomes-imprinting 2018-03-13

import os
import pandas as pd
import re
import argparse
import json
import datetime
import glob
import subprocess

def read_input_fof(fof_file):
    '''
    Function to get merged fastq files fof
    Input: path to the project of interest
    Output: fof file to do the mapping
    '''
    # Read fof file

    with open(fof_file, "r") as fi:
        fof = fi.read().splitlines()
    return(fof)

def merged_fastq_dataframe(fastq_files):
    '''
    Function to take lines from fof file to create a sorted pandas dataframe with info from merged fastqs
    Input: list of fastq files read from the fof file
    Ouptut: sorted dataframe with the fastq name metadata
    '''
    #<TODO> join this function with the fastq_dataframe function from merge_fastq.py
        #Ideally we would only create the rex variable outside and pass it as input to the function
        #Then we create an if statement and chose the corresponding columns for the dataframe

    # Regular expression to create an annonymous dictionary
    dict_list = []

    rex = re.compile(r"{}(?P<sample>ATA_\d+)_R(?P<read>\d+)_merged.fastq.gz".format("fastq/"))

    for fastq in fastq_files:
        m = rex.match(fastq)
        dicc = m.groupdict()
        dicc['fastq_file'] = fastq.split('/')[-1]
        dicc['fastq_path'] = str(fastq)
        dict_list.append(dicc)

    # Read fastq info array into a pandas dataframe
    
    df_fastq = pd.DataFrame(dict_list, columns = ['sample', 'read', 'fastq_file', 'fastq_path'])

    # Sort by sample and read

    df_fastq_sorted = df_fastq.sort_values(by=['fastq_path'], ascending=[True]).reset_index(drop=True)

    return(df_fastq_sorted)

def cmd_map_fastq(df_fastq_sorted):
    '''
    Function to create commands to map merged fastq files, we have
        a R1 and R2 fastq corresponding to each sample
    Input: dataframe with metadata from merged fastq files
    Output: files with commands for mapping with BWA
    '''

    cmd_sh = datetime.datetime.now().strftime("cmd_bwa_mem_%Y%m%d_%H-%M-%S.sh")

    # Create list by sample

    samples = list(df_fastq_sorted['sample'].unique())

    # Set templates

    for i in range(len(samples)):
        sample_fastqs = list(df_fastq_sorted[df_fastq_sorted['sample']==samples[i]]['fastq_path'])
        R1_list = sample_fastqs[0]
        R2_list = sample_fastqs[1]
        output_bam_name = df_fastq_sorted['sample'].unique()[i]

        # Create template for bwa
        
        cmd_bwa = """time bwa mem -M -L 5 -t 10 {}.gz -R "@RG\\tID:{}\\tPL:ILLUMINA\\tSM:{}\\tDS:ref=37d5\\tCN:UGDG\\tDT:{}\\tPU:{}" {} {} 2> {}{}_bwa_mem.err | samtools sort -@ 10 -O bam -o {}{}_sorted.bam; time samtools index {}{}_sorted.bam""".format("/nfs/qnapugdg8tb3a/nfs_references/human/GRCh38_GCA-000001405-15/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna", output_bam_name, output_bam_name, datetime.datetime.now().strftime("%Y-%m-%d"), output_bam_name, R1_list, R2_list, "./general_bams/", output_bam_name, "./general_bams/", output_bam_name, "./general_bams/", output_bam_name)

        # Write the cmd.sh file

        with open('{}{}'.format("./", cmd_sh), 'a') as cmd_file:
            cmd_file.write(cmd_bwa + '\n')
    print("[INFO]: CMD_FILE - {}{}".format("./", cmd_sh))

    return(cmd_sh)


def run_parallel(cmd_sh):
    '''
    Function to run the cmd file in parallel
    Input: cmd file and project path to find it
    Output: merged fastq files
    '''
    log_str = datetime.datetime.now().strftime("bwa_mem_%Y%m%d_%H-%M-%S.log")

    cmd = "parallel --joblog {}{} -j 5 :::: {}{}".format("./general_bams/", log_str, "./", cmd_sh)

    print("[CMD]: " + cmd)

    subprocess.call(cmd + " 2> /dev/null", shell = True)

def main():
    '''
    Function to sort the order of functions to use
    '''
    fastq_files = read_input_fof("./fastq_files.fof")

    # Create fastq files dataframe

    df_fastq_sorted = merged_fastq_dataframe(fastq_files)

    # Create BWA MEM commands for mapping

    cmd_sh = cmd_map_fastq(df_fastq_sorted)

    #run_parallel(cmd_sh)

if __name__ == '__main__':
    main()
