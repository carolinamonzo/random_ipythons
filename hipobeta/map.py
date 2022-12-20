#! usr/bin/python3.6

# 2018-01-03, Carolina Monzo

import os
import pandas as pd
import re


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

    rex = re.compile(r"\./fastq_merged/fastq_merged_(seqtk|cutadapt)/(?P<sample>DHB-\d+\.\d+)_R(?P<read>\d+).+")

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



def map_fastq(df_fastq_sorted):
    '''
    Function to create commands to map merged fastq files, we have
        a R1 and R2 fastq corresponding to each sample
    Input: dataframe with metadata from merged fastq files
    Output: files with commands for mapping with BWA and SNAP
    '''
    # Create list by sample

    samples = list(df_fastq_sorted['sample'].unique())

    # Get trimming tool

    trimming_tool = list(df_fastq_sorted['fastq_path'])[0].split('/')[2].split('_')[-1]

    # Set templates

    for i in range(len(samples)):
        sample_fastqs = list(df_fastq_sorted[df_fastq_sorted['sample']==samples[i]]['fastq_path'])
        R1_list = sample_fastqs[0]
        R2_list = sample_fastqs[1]
        output_bam_name = df_fastq_sorted['sample'].unique()[i]

        # Create template for bwa
        
        cmd_bwa = "time bwa mem -M -L 5 -t 5 /storage/ethernus_miscellanea/References/human/1kg/hs37d5/hs37d5.fa.gz \
                {} {} 2> ./mapping_bwa_mem/bwa_{}/{}_bwa_mem.err | samtools sort -O bam -o \
                ./mapping_bwa_mem/bwa_{}/{}.{}.sorted.bam; time samtools index ./mapping_bwa_mem/bwa_{}/{}.{}.sorted.bam".format(R1_list, R2_list, trimming_tool, output_bam_name, trimming_tool, output_bam_name, trimming_tool, trimming_tool, output_bam_name, trimming_tool)

        # Write the cmd.sh file

        with open('cmd_bwa_mem_{}.sh'.format(trimming_tool), 'a') as cmd_file:
            cmd_file.write(cmd_bwa + '\n')

        # Create template for SNAP

        cmd_snap = "time /home/cmc/software/snap/snap-aligner paired /storage/ethernus_miscellanea/scratch_local/workspace/cmc_projects_tmp/hipobeta_work/snap_reference_genome/ {} {} -o -bam ./mapping_snap/snap_{}/{}.{}.unsorted.bam".format(R1_list, R2_list, trimming_tool, output_bam_name, trimming_tool)
        # Create template for samtools sort and remove unsorted bam file

        cmd_samtools_sort = "time samtools sort ./mapping_snap/snap_{}/{}.{}.unsorted.bam -O bam -o ./mapping_snap/snap_{}/{}.{}.sorted.bam; rm -f ./mapping_snap/snap_{}/{}.{}.unsorted.bam; time samtools index ./mapping_snap/snap_{}/{}.{}.sorted.bam".format(trimming_tool, output_bam_name, trimming_tool, trimming_tool, output_bam_name, trimming_tool, trimming_tool, output_bam_name, trimming_tool, trimming_tool, output_bam_name, trimming_tool)

        # Write the cmd.sh file

        with open('cmd_snap_{}.sh'.format(trimming_tool), 'a') as cmd_file:
            cmd_file.write(cmd_snap + '; ' + cmd_samtools_sort + '\n')



def main():
    '''
    Function to sort the order of functions to use
    '''

    fastq_fof = 'merged_fastq_seqtk.fof'

    with open(fastq_fof) as fof:
        fastq_files = fof.read().splitlines()

    # Create fastq files dataframe

    df_fastq_sorted = merged_fastq_dataframe(fastq_files)

    # Create BWA, SNAP and STAR commands for mapping

    map_fastq(df_fastq_sorted)



if __name__ == '__main__':
    main()
