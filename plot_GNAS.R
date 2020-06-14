library(data.table)
library(Gviz)
library(GenomicFeatures)

options(Gviz.ucscUrl="http://genome-euro.ucsc.edu/cgi-bin/")

## PLOTTING GNAS
#atrack = UcscTrack(track="CpG Islands", chromosome="chr20", genome="hg19", start="chromStart", end="chromEnd", name="CpG Islands", from=57225763, to=57486250)
# Good range

#itrack = IdeogramTrack(genome = "hg19", chromosome = "chr20")
#gtrack = GenomeAxisTrack()
#grtrack = UcscTrack(genome = "hg19", chromosome = "chr20", track = "ensGene", from = 57225763, to = 57486250, 
#                    trackType = "GeneRegionTrack", rstarts = "exonStarts", rends = "exonEnds", gene = "name", 
#                    symbol = "name2", transcript = "name", strand = "strand", name = "GNAS complex locus", stacking = "pack")

plotTracks(list(itrack, gtrack, atrack, grtrack), from = 57225763, to = 57486250, showBandId = TRUE,
           sizes = c(0.1, 0.1, 0.15, 0.65))



## PLOTTING SNRPN
atrack = UcscTrack(track="CpG Islands", chromosome="chr15", genome="hg19", start="chromStart", end="chromEnd", 
                   name="CpG Islands", from=25068780, to=25223870)
# Good range

itrack = IdeogramTrack(genome = "hg19", chromosome = "chr15")
gtrack = GenomeAxisTrack()
grtrack = UcscTrack(genome = "hg19", chromosome = "chr15", track = "ensGene", from = 25068780, to = 25223870, 
                    trackType = "GeneRegionTrack", rstarts = "exonStarts", rends = "exonEnds", gene = "name", 
                    symbol = "name2", transcript = "name", strand = "strand", name = "SNRPN small nuclear ribonucleoprotein polypeptide N", stacking = "pack")

plotTracks(list(itrack, gtrack, atrack, grtrack), from = 25068780, to = 25223870, showBandId = TRUE,
           sizes = c(0.1, 0.1, 0.15, 0.65))