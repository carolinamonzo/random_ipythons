#! /usr/bin/python3.5

from biomart import BiomartServer

server = BiomartServer("http://grch37.ensembl.org/biomart")

    # set verbose to true
server.verbose = True

    # show server databases
    #server.show_databases()

    # Show server datasets
#server.show_datasets()

    # User the ensemble genes dataset
genes = server.datasets['hsapiens_gene_ensembl']

    # Show all available filters and atributes of the hsapiens gene dataset
#genes.show_filters()
#genes.show_attributes()

response = genes.search({
    'filters':{
        'chromosomal_region':["10:100148059:100148060", "10:100456043:100456044", "10:100385714:100385715"]
    },
    'attributes':["chromosome_name", "start_position", "end_position", "ensembl_gene_id", "external_gene_name"]
})


            # response format in TSV
for line in response.iter_lines():
    line = line.decode('utf-8')

    print(line.split("\t"))


