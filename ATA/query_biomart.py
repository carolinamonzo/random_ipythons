#! /usr/bin/python3

from biomart import BiomartServer
import argparse
import time

def parseArguments():
    '''

    '''

    # Create argument parser class
    parser = argparse.ArgumentParser(description = "Get gene names from positions in biomart")

    # Define arguments to parse
    parser.add_argument("--input", "-i", required=True, type=str, help="Argument to set the name of the file with possitions")

    # Call for arguments
    args = parser.parse_args()

    return(args)

def get_positions(input_file):
    """

    """

    with open(input_file, 'r') as infile:
        inf = infile.read().splitlines()

    new_inf = []

    # Create list with positions to ask biomart

    for e in inf:

        if "hs37d5" and "GL00" not in e:
            plus = int(str(e).split(" ")[2]) + 1
            ei = ":".join(str(e).split(" ")[1:])        

            new_inf.append(ei + ":" + str(plus))

        else:
            continue

    new2_inf = sorted(set(new_inf))

    return(new2_inf)

def ask_biomart(new_inf, infile):
    """
    """
    server = BiomartServer("http://grch37.ensembl.org/biomart")

    # set verbose to true
    #server.verbose = True
    server.verbose=False

    # show server databases
    #server.show_databases()

    # Show server datasets
    #server.show_datasets()

    # User the ensemble genes dataset
    genes = server.datasets['hsapiens_gene_ensembl']

    # Show all available filters and atributes of the hsapiens gene dataset
    #genes.show_filters()
    #genes.show_attributes()

    with open(infile.replace("pos", "genes"), "a") as outfile:

        for i in range(0, len(new_inf), 100):
 
            response = genes.search({
                'filters':{
                    'chromosomal_region':new_inf[i:i+100]
                },
                'attributes':["chromosome_name", "start_position", "end_position", "ensembl_gene_id", "external_gene_name"]
            })


            # response format in TSV
            for line in response.iter_lines():
                line = line.decode('utf-8')

            # If the position we have asked for has information, print it
            # in the same line

                outfile.write(line + "\n")
                #print(line.split("\t"))
            


def main():
    args = parseArguments()

    new_inf = get_positions(args.input)

    #print(new_inf)

    ask_biomart(new_inf, args.input)


if __name__ == "__main__":

    main()
