#! /usr/bin/python3.5

import argparse

def parseArguments():
    '''

    '''

    # Create argument parser class
    parser = argparse.ArgumentParser(description = "Get gene names from positions in biomart")

    # Define arguments to parse
    parser.add_argument("--input", "-i", required=True, type=str, help="Argument to set the name of the file with possitions")

    parser.add_argument("--genes", "-g", required=True, type=str, help="Argument to set the name of the file with the genes")

    # Call for arguments
    args = parser.parse_args()

    return(args)


def read_files(ireads, igenes):

    new_reads = []

    with open(ireads, 'r') as inreads:
        reads = inreads.read().splitlines()
        for e in reads:

            if "GL0" not in e:
                if "hs37d5" not in e:
                    if "NC" not in e:

                        new_reads.append(e.split(" "))

    new_genes = []

    with open(igenes, "r") as ingenes:
        genes = ingenes.read().splitlines()
        for en in genes:
            new_genes.append(en.split("\t"))

    return(new_reads, new_genes)


def count_reads(reads, genes, filename):

    with open(filename.replace("pos", "counts"), "w") as fileout:

        count = 0
        i = 0
        for j in range(0, len(reads)):

            if reads[j][1] == genes[i][0]:

                if reads[j][2] >= genes[i][1] and reads[j][2] <= genes[i][2]:

                    count += 1
                    fileout.write("\t".join(genes[i]) + "\t" + str(count) + "\t" + str(reads[j][0]) + "\n" )

                elif reads[j][2] > genes[i][2]:

                #print("dentro_next")
                    if i < len(genes):
                        i += 1
                        count = 0
                    elif i > len(genes):
                        break

            elif reads[j][1] != genes[i][0]:

                if i < len(genes):
                    i +=1

                elif i > len(genes):
                    break


def main():
    
    args = parseArguments()

    reads, genes = read_files(args.input, args.genes)

    count_reads(sorted(reads, key=lambda x: x[1:]), sorted(genes), args.input)
   


if __name__ == '__main__':
    main()
