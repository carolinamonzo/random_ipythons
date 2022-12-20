#! /usr/bin/python

def read_fasta(filename):
    sequences = []
    identifiers = []
    genome = ''

    with open(filename, 'r') as f:
        for line in f:
            if line.startswith(">"):
                identifiers.append(line.rstrip().replace(">", ""))
                sequences.append(genome)
                genome = ''
            else:
                genome += line.rstrip()
        sequences.append(genome)
        sequences.pop(0)                    
    if len(sequences) == len(identifiers):
        #dictionary = dict(zip(identifiers, sequences))
        return(identifiers, sequences)

    else:
        print("ERROR: ID list length = {}/n Genome list length = {}".format(len(identifiers), len(sequences)))

def parse_identifiers(identifiers):
    formated = [">Viruses", "ssDNAviruses", "Anelloviridae"]
    new_seq = []
    for s in identifiers:
        seq = s.replace(" ", "_")
        if "Torque_teno_mini_virus" in seq:
            new_seq.append(";".join(formated + ["Betatorquevirus", seq + ';']))
            
        elif "Torque_teno_midi_virus" in seq:
            new_seq.append(";".join(formated + ["Gammatorquevirus", seq + ';']))
            
        elif "Torque_teno_virus" in seq:
            new_seq.append(";".join(formated + ["Alphatorquevirus", seq + ";"]))

        elif "like" in seq:
            new_seq.append(";".join(formated + ["Revise", seq + ";"]))        
 
        else:
            new_seq.append(";".join(formated + ["Unclassified", seq + ";"]))

    print(len(identifiers))
    print(len(new_seq))

    return(new_seq)


def main():

    identifiers, sequences = read_fasta("all_torques.fa")
    identifiers = parse_identifiers(identifiers)

    with open("formated_ttv-genomes.fa", 'w') as fh:
        for i in range(len(identifiers)):

            fh.write(identifiers[i] + '\n')
            fh.write(sequences[i] + '\n')


if __name__ == "__main__":
    main()
