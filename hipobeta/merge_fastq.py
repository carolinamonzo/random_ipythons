#! usr/bin/python3.6

# 2017-12-31, Carolina Monzo

import os
import re   
import pandas as pd

def fastq_dataframe(fastq_files):
    '''
    Function to take lines from fof file and create a sorted pandas dataframe with info from fastqs
    Input: list of fastq files read from the fof file
    Output: sorted dataframe with the fastq name metadata
    '''

    # Regular expression to create an annonymous dictionary
    dict_list = []
    rex = re.compile(r"\./fastq_trimmed/(seqtk|cutadapt)_trimming/(?P<sample>DHB-\d+\.\d+)_L(?P<lane>\d+)_R(?P<read>\d+).+")
    for fastq in fastq_files:
        m = rex.match(fastq)
        dicc = m.groupdict()
        dicc['fastq_file'] = fastq.split('/')[-1]
        dicc['fastq_path'] = str(fastq)
        dict_list.append(dicc)

    # Read fastq info array into a pandas dataframe
    df_fastq = pd.DataFrame(dict_list, columns = ['sample', 'lane', 'read', 'fastq_file', 'fastq_path'])
    # Sort by sample, lane and read
    df_fastq_sorted = df_fastq.sort_values(by=["fastq_path"], ascending=[True]).reset_index(drop=True)

    return(df_fastq_sorted)


def zcat_fastq(df_fastq_sorted):
    '''
    Function to generate the zcat commands, it selects R1 and R2 from sorted samples
    Input: dataframe of fastq sorted names
    Output: cmd file with the zcat commands of ordered fastqs to merge
    '''

    # Create list by sample
    samples = list(df_fastq_sorted['sample'].unique())

    # Get trimming tool
    
    trimming_tool = list(df_fastq_sorted['fastq_file'])[0].split('-')[2]
    
    # Set zcat template

    for i in range(len(samples)):
            
            sample_fastqs = list(df_fastq_sorted[df_fastq_sorted['sample']==samples[i]]['fastq_path'])
            R1_list = sample_fastqs[0::2]
            R2_list = sample_fastqs[1::2]
            merged_fastq_name = df_fastq_sorted['sample'].unique()[i]
            R1_str = "zcat {} | gzip > ./fastq_merged/fastq_merged_{}/{}_R1_merged.fastq.gz".format(' '.join(R1_list), trimming_tool, merged_fastq_name)
            R2_str = "zcat {} | gzip > ./fastq_merged/fastq_merged_{}/{}_R2_merged.fastq.gz".format(' '.join(R2_list), trimming_tool, merged_fastq_name)
                                    
            # Write the cmd.sh file
            with open('cmd_zcat_{}_fastq.sh'.format(trimming_tool), 'a') as cmd_file:
                cmd_file.write(R1_str + '\n')
                cmd_file.write(R2_str + '\n')


def main():
    '''
    Function to sort the orders
    '''
    fastq_fof = 'fastq_seqtk_trimmed.fof'

    with open(fastq_fof) as fof:
            fastq_files = fof.read().splitlines()


    df_fastq_sorted = fastq_dataframe(fastq_files)

    zcat_fastq(df_fastq_sorted)



if __name__ == '__main__':
    main()
