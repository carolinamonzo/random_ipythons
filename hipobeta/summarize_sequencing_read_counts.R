##
###
##


library ('ggplot2')

## # obtain the read counts
## pmg@Illuina-pipeline[~/illumina_runs/120719_H114_0214_AC0RJ6ACXX/Unaligned-1mismatch]36
## $  for x in $(ls | grep TARIN2-Pool); do echo "folder $x"; for sample in $(ls $x | grep -P 'Sample'); do echo "$x/$sample"; zcat $x/$sample/*gz | wc -l | perl -lne 'print "reads '$x'\t'$sample'\t". $_/4' ;done ;done | grep reads > reads_counts_tarin2-1_mismatch_AC0RJ6ACXX.tab

# no mismatch
#flowcells
fc <- c('AC0RJ6ACXX','Ad0hruacxx')
setwd("/home/anab/javier_labo/tarin/tarin2")
run1 <- read.delim("/home/anab/javier_labo/tarin/tarin2/reads_counts_tarin2-1_mismatch_Ad0hruacxx.tab", header=FALSE)
run2 <- read.delim("/home/anab/javier_labo/tarin/tarin2/reads_counts_tarin2-1_mismatch_AC0RJ6ACXX.tab", header=FALSE)

colum_names <- c('flowcell', 'pool-lane','sample','counts')
names(run1) <- colum_names
names(run2) <- colum_names

# sum by factor
library('plyr')
sum_run1 <- ddply(run1, 'sample', summarize, sum_run1=sum(counts)) 
sum_run2 <- ddply(run2, 'sample', summarize, sum_run2=sum(counts))

all_counts <- merge(sum_run1, sum_run2, by='sample')
all_counts$total <- rowSums(all_counts[2:3])
all_counts <- all_counts[rev(order(all_counts$total)),]

# add the pool name
all_counts_with_pool <- merge(all_counts, run2[,c('sample', 'pool-lane')])



write.table(all_counts_with_pool, file='tarin2_summary_counts.tab',  row.names=FALSE, sep="\t", dec=',', quote=FALSE)
#write.csv2(all_counts_with_pool, file='tarin2_summary_counts.csv',  row.names=FALSE)




ggplot(dat.nm, aes(x = lane, y = read_counts)) +
  geom_boxplot(aes(fill = lane), alpha = 0.3, outlier.color =
               NA) +
  geom_point(position = position_jitter(width = 0.05),
             colour = "blue", fill = "blue") +
  labs(x = "lane", y = "read_counts", fill = "lane"

ggplot(mtcars, aes(x = factor(vs), y = mpg)) +
  geom_boxplot(aes(fill = factor(vs)), alpha = 0.3, outlier.color =
               NA) +
  geom_dotplot(binaxis = "y", stackdir = "center",
               position = "dodge", colour = "blue", fill = "blue") +
  labs(x = "vs", y = "mpg", fill = "vs")

