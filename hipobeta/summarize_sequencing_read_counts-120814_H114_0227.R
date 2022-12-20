##
###
##


library ('ggplot2')

date <-Sys.Date()

## # obtain the read counts
## pmg@Illuina-pipeline[~/illumina_runs/120719_H114_0214_AC0RJ6ACXX/Unaligned-1mismatch]36
## $  for x in $(ls | grep TARIN2-Pool); do echo "folder $x"; for sample in $(ls $x | grep -P 'Sample'); do echo "$x/$sample"; zcat $x/$sample/*gz | wc -l | perl -lne 'print "120814_H114_0227\t'$x'\t'$sample'\t". $_/4' ;done ;done | grep reads > reads_counts_tarin2-1_mismatch_120814_H114_0227.tab


############################################
##   TARIN 2 
############################################

tarin2.last_run <- "120928_H114_0230"

# no mismatch
#flowcells
setwd("/home/anab/javier_labo/tarin/tarin2")
run1 <- read.delim("/home/anab/javier_labo/tarin/tarin2/reads_counts_tarin2-1_mismatch_Ad0hruacxx.tab", header=FALSE)
run2 <- read.delim("/home/anab/javier_labo/tarin/tarin2/reads_counts_tarin2-1_mismatch_AC0RJ6ACXX.tab", header=FALSE)
run3 <- read.delim("/home/anab/javier_labo/tarin/tarin2/120814_H114_0227_Ac127uacxx/count_reads-tarin2_by_sample-120814_H114_0227.tab", header=FALSE)
run4 <- read.delim("/home/anab/javier_labo/tarin/tarin2/reads_counts-Tarin_2-1_mismatch-120828_H114_0228.tab", header=FALSE)
run5 <- read.delim("/home/anab/javier_labo/tarin/tarin2/reads_counts-Tarin_2-1_mismatch-120928_H114_0230.tab", header=FALSE)

column_names <- c('flowcell', 'pool-lane','sample','counts')
names(run1) <- column_names
names(run2) <- column_names
names(run3) <- column_names
names(run4) <- column_names
names(run5) <- column_names
  
# sum by elements of a column (just in case some samples are repeated). It does not need to be a factor
library('plyr')
sum_run1 <- ddply(run1, 'sample', summarize, sum_run1=sum(counts)) 
sum_run2 <- ddply(run2, 'sample', summarize, sum_run2=sum(counts))
sum_run3 <- ddply(run3, 'sample', summarize, sum_run3=sum(counts))
sum_run4 <- ddply(run4, 'sample', summarize, sum_run4=sum(counts))
sum_run5 <- ddply(run5, 'sample', summarize, sum_run5=sum(counts))

# sorting the factors
# http://stackoverflow.com/questions/3744178/ggplot2-sorting-a-plot
# http://wiki.stdout.org/rcookbook/Manipulating%20data/Changing%20the%20order%20of%20levels%20of%20a%20factor/

## selecting the tarin 2 pool4

t2p4m <- run3[run3$'pool-lane'=="TARIN2-Pool4-mod",]
t2p4m$sample_sorted <- factor(t2p4m$sample, levels=unique(as.character(t2p4m[rev(order(t2p4m$counts)),'sample'])) )

## plotting the read counts
ggplot(t2p4m, aes(x= t2p4m$sample_sorted, y=t2p4m$counts)) + geom_point() + opts(axis.text.x=theme_text(angle=90, hjust=1)) + xlab('Sample')+ ylab('read_counts') +opts(title='Tarin 2, pool-4_mod 120814_H114_0227_Ac127uacxx')

ggplot(t2p4m, aes(x= t2p4m$sample_sorted, y=t2p4m$counts)) + geom_point() + xlab('Sample')+ ylab('read_counts') +opts(title='Tarin 2, pool-4_mod 120814_H114_0227_Ac127uacxx') + coord_flip() 


## Selecting the tarin2 repesca
t2_remix <-  run3[run3$'pool-lane'=="TARIN2-Pool-mix",]
t2_remix$sample_sorted <- factor(t2_remix$sample, levels=unique(as.character(t2_remix[rev(order(t2_remix$counts)),'sample'])) )
ggplot(t2_remix, aes(x= t2_remix$sample_sorted, y=t2_remix$counts)) + geom_point() + opts(axis.text.x=theme_text(angle=90, hjust=1)) + xlab('Sample')+ ylab('read_counts') +opts(title='Tarin 2, pool-1-2-3 remix 120814_H114_0227_Ac127uacxx')


## Selecting the tarin2 repesca TARIN2-Poolmix-2
t2_pm2 <-  run4[run4$'pool-lane'=="TARIN2-Poolmix-2",]
t2_pm2$sample_sorted <- factor(t2_pm2$sample, levels=unique(as.character(t2_pm2[rev(order(t2_pm2$counts)),'sample'])) )
ggplot(t2_pm2, aes(x= t2_pm2$sample_sorted, y=t2_pm2$counts)) + geom_point() + opts(axis.text.x=theme_text(angle=90, hjust=1)) + xlab('Sample')+ ylab('read_counts') +opts(title='Tarin 2, pool remix2 120828_H114_0228_AD1B47ACXX')

## Selecting the tarin2 pool4 TARIN2-Pool4-mod2
t2_p4m2 <-  run4[run4$'pool-lane'=="TARIN2-Pool4-mod2",]
t2_p4m2$sample_sorted <- factor(t2_p4m2$sample, levels=unique(as.character(t2_p4m2[rev(order(t2_p4m2$counts)),'sample'])) )
ggplot(t2_p4m2, aes(x= t2_p4m2$sample_sorted, y=t2_p4m2$counts)) + geom_point() + opts(axis.text.x=theme_text(angle=90, hjust=1)) + xlab('Sample')+ ylab('read_counts') +opts(title='Tarin 2, pool-4m2 120828_H114_0228_AD1B47ACXX')



## all values
##  
##  run3$sample_sorted <- factor(run3$sample, levels=unique(as.character(run3[rev(order(run3$counts)),'sample'])) )
##  ggplot(run3, aes(x= run3$sample_sorted, y=run3$counts)) + geom_point() + xlab('Sample')+ ylab('read_counts') +opts(title='Tarin 2, 120814_H114_0227_Ac127uacxx') + coord_flip() 
##  
##  
##  all_counts <- merge(sum_run1, sum_run2, by='sample')
##  all_counts <- merge(all_counts, sum_run3, by='sample')
##  all_counts$total <- rowSums(all_counts[2:4])
##  all_counts <- all_counts[rev(order(all_counts$total)),]
##  
##  # add the pool name
##  all_counts_with_pool <- merge(all_counts, run2[,c('sample', 'pool-lane')])
##  all_counts_with_pool <- all_counts_with_pool[rev(order(all_counts_with_pool$total)),]
##  
##  
##  
##  write.table(all_counts_with_pool, file='tarin2_summary_counts-runs-1-2-3.tab',  row.names=FALSE, sep="\t", dec=',', quote=FALSE)
##  #write.csv2(all_counts_with_pool, file='tarin2_summary_counts.csv',  row.names=FALSE)
##  


##  all counts
###########################

all_counts <- merge(sum_run1, sum_run2, by='sample', all=TRUE)
all_counts <- merge(all_counts, sum_run3, by='sample', all=TRUE)
all_counts <- merge(all_counts, sum_run4, by='sample', all=TRUE)
all_counts <- merge(all_counts, sum_run5, by='sample', all=TRUE)


## remember to expand the all_counts one column when adding a new run
all_counts$total <- rowSums(all_counts[2:6], na.rm=TRUE)
all_counts <- all_counts[rev(order(all_counts$total)),]

# add the pool name
all_counts_with_pool <- merge(all_counts, run5[,c('sample', 'pool-lane')], all=TRUE)
all_counts_with_pool <- all_counts_with_pool[rev(order(all_counts_with_pool$total)),]

all_counts_with_pool$'pool-lane' <- ifelse(test=is.na(all_counts_with_pool$'pool-lane'), yes='previous_pool', no=as.character(all_counts_with_pool$'pool-lane'))

all_counts_with_pool$sample_sorted <- factor(all_counts_with_pool$sample, levels=unique(as.character(all_counts_with_pool[rev(order(all_counts_with_pool$total)),'sample'])) )
#ggplot(all_counts_with_pool, aes(x= all_counts_with_pool$sample_sorted, y=all_counts_with_pool$total)) + geom_point() + xlab('Sample')+ ylab('read_counts') +opts(title='Tarin 2 all three runs') + coord_flip()


ggplot(all_counts_with_pool, aes(x=sample_sorted, y=total, colour=all_counts_with_pool$'pool-lane')) + geom_point() + xlab('Sample')+ ylab('read_counts') +opts(title=paste('Tarin-2 read counts after run ',tarin2.last_run,' ',date, sep='')) + opts(axis.text.x=theme_text(angle=90, hjust=1)) + labs(colour='Pools') + geom_hline(yintercept=c(10000000, 7000000))

pdf(paste('tarin2-reads_count_afer_run-',tarin2.last_run,'_',date,'.pdf', sep=''),  width = 11.7, height = 8.3 ) #paper='a4r')
ggplot(all_counts_with_pool, aes(x=sample_sorted, y=total, colour=all_counts_with_pool$'pool-lane')) + geom_point() + xlab('Sample')+ ylab('read_counts') +opts(title='Tarin_2 runs 1-5') + opts(axis.text.x=theme_text(angle=90, hjust=1, size=8)) + labs(colour='Pools') + geom_hline(yintercept=c(10000000, 7000000))
dev.off()

write.table(all_counts_with_pool, file=paste('tarin2_summary_counts_afer_run-',tarin2.last_run,'_',date,'.tab', sep=''),  row.names=FALSE, sep="\t", dec=',', quote=FALSE)


## how many samples with reads less than 10 milion

table(subset(all_counts_with_pool, total<10000000)$'pool-lane')

## TARIN2-Pool1 TARIN2-Pool2 TARIN2-Pool3 TARIN2-Pool4 
##            7           10            9           24 

## list of samples/pool and how many reds they need for 100000

lt_10m <- subset(all_counts_with_pool, total<10000000)
lt_10m.sorted <- lt_10m[order(lt_10m$'pool-lane', lt_10m$total), c("pool-lane", "sample","total")]

write.table(lt_10m.sorted, file=paste('tarin2-samples_less_10_million_afer_run-',tarin2.last_run,'_',date,'.tab', sep=''),  row.names=FALSE, sep="\t", dec=',', quote=FALSE)



needed_reads.sorted <- data.frame(lt_10m.sorted$'pool-lane', lt_10m.sorted$sample, -1*(lt_10m.sorted$total - 10000000)/1000000)
colnames(needed_reads.sorted) = c("pool","sample","missed_reads")


write.table(needed_reads.sorted, file=paste('tarin2-number_of_reads_still_needed_by_sample_afer_run-',tarin2.last_run,'_',date,'.tab', sep=''),  row.names=FALSE, sep="\t", dec=',', quote=FALSE)


############################################
##   TARIN 1
############################################

# no mismatch
#flowcells

tarin1.last_run <- "120928_H114_0230"
#tarin1.last_run <- "120828_H114_0228"

setwd("/home/anab/javier_labo/tarin/tarin1")
run1 <- read.delim("/home/anab/javier_labo/tarin/tarin1/tarin_1_previous_runs.tab", header=FALSE)
run2 <- read.delim("/home/anab/javier_labo/tarin/tarin1/reads_counts-Tarin_1-1_mismatch-120814_H114_0227.tab", header=FALSE)
run3 <- read.delim("/home/anab/javier_labo/tarin/tarin1/reads_counts-Tarin_1-1_mismatch-120828_H114_0228.tab", header=FALSE)
run4 <- read.delim("/home/anab/javier_labo/tarin/tarin1/reads_counts-Tarin_1-1_mismatch-120928_H114_0230.tab", header=FALSE)

column_names <- c('flowcell', 'pool-lane','sample','counts')
names(run1) <- column_names
names(run2) <- column_names
names(run3) <- column_names
names(run4) <- column_names


# sum by factor
library('plyr')
sum_run1 <- ddply(run1, 'sample', summarize, sum_run1=sum(counts)) 
sum_run2 <- ddply(run2, 'sample', summarize, sum_run2=sum(counts))
sum_run3 <- ddply(run3, 'sample', summarize, sum_run3=sum(counts))
sum_run4 <- ddply(run4, 'sample', summarize, sum_run4=sum(counts))


# sorting the factors
# http://stackoverflow.com/questions/3744178/ggplot2-sorting-a-plot
# http://wiki.stdout.org/rcookbook/Manipulating%20data/Changing%20the%20order%20of%20levels%20of%20a%20factor/


## ## Selecting the tarin2 repesca
## t2_remix <-  run3[run3$'pool-lane'=="TARIN2-Pool-mix",]
## t2_remix$sample_sorted <- factor(t2_remix$sample, levels=unique(as.character(t2_remix[rev(order(t2_remix$counts)),'sample'])) )
## ggplot(t2_remix, aes(x= t2_remix$sample_sorted, y=t2_remix$counts)) + geom_point() + opts(axis.text.x=theme_text(angle=90, hjust=1)) + xlab('Sample')+ ylab('read_counts') +opts(title='Tarin 2, pool-1-2-3 remix 120814_H114_0227_Ac127uacxx')
## 

## get indexes for samples
##########################

tarin1_sample_sheet<-read.csv("/home/anab/javier_labo/tarin/tarin1-sample-index.csv", header=FALSE)
## head(tarin1_sample_sheet)
tarin1_sample_index <-  unique(tarin1_sample_sheet[,c("V3","V6")])
names(tarin1_sample_index) <- c("sample","index")
cat(paste("# number of samples unique for Tarin1:",nrow(tarin1_sample_index),"\n"))

##  all counts
###########################

all_counts <- merge(sum_run1, sum_run2, by='sample', all=TRUE)
all_counts <- merge(all_counts, sum_run3, by='sample', all=TRUE)
all_counts <- merge(all_counts, sum_run4, by='sample', all=TRUE)



all_counts$total <- rowSums(all_counts[2:5], na.rm=TRUE)
all_counts <- all_counts[rev(order(all_counts$total)),]

## add the pool name
all_counts_with_pool <- merge(all_counts, run4[,c('sample', 'pool-lane')], all=TRUE)
all_counts_with_pool <- all_counts_with_pool[rev(order(all_counts_with_pool$total)),]


## add the index
all_counts_with_pool <- merge(all_counts_with_pool, tarin1_sample_index, all=TRUE)

##
### sort by total,  plot and write table
##
## sort the factor sample
all_counts_with_pool$sample <- factor(all_counts_with_pool$sample, levels=unique(as.character(all_counts_with_pool[rev(order(all_counts_with_pool$total)),'sample'])) )

ggplot(all_counts_with_pool, aes(x=sample, y=total, colour=all_counts_with_pool$'pool-lane')) + geom_point() + xlab('Sample')+ ylab('read_counts') +opts(title=paste('Tarin_1 after-run ', tarin1.last_run, sep='')) + opts(axis.text.x=theme_text(angle=90, hjust=1)) + labs(colour='Pools') + geom_hline(yintercept=c(10000000, 7000000))

write.table(all_counts_with_pool, file=paste('tarin1_summary_counts-after_run-',tarin1.last_run,'_',date,'.tab', sep=''),  row.names=FALSE, sep="\t", dec=',', quote=FALSE)


##
### run2 and sort by counts
##

## get run2 dataframe
run2_counts_with_pool <- all_counts_with_pool[which(!is.na(all_counts_with_pool$sum_run2)),c(1,3,7)]
run2_counts_with_pool <-merge(run2_counts_with_pool,run2[,c(2,3)])

## sort sample by count
run2_counts_with_pool$sample <- factor(run2_counts_with_pool$sample, levels=unique(as.character(run2_counts_with_pool[rev(order(run2_counts_with_pool$sum_run2)),'sample'])) )

ggplot(run2_counts_with_pool, aes(x=sample, y=sum_run2, colour=run2_counts_with_pool$'pool-lane')) + geom_point() + xlab('Sample')+ ylab('read_counts') +opts(title='Tarin_1 last-runs') + opts(axis.text.x=theme_text(angle=90, hjust=1)) + labs(colour='Pools') + geom_hline(yintercept=c(10000000, 7000000))

run2_counts_with_pool$sample <- factor(run2_counts_with_pool$sample, levels=unique(as.character(run2_counts_with_pool[rev(order(run2_counts_with_pool$index)),'sample'])) )

ggplot(run2_counts_with_pool, aes(x=sample, y=sum_run2, colour=run2_counts_with_pool$index)) + geom_point() + xlab('Sample')+ ylab('read_counts') +opts(title='Tarin_1 last-runs') + opts(axis.text.x=theme_text(angle=90, hjust=1)) + labs(colour='index') + geom_hline(yintercept=c(10000000, 7000000))

##
### run3 and sort by counts
##

## get run3 dataframe
run3_counts_with_pool <- all_counts_with_pool[which(!is.na(all_counts_with_pool$sum_run3)),c(1,4,7)]
run3_counts_with_pool <-merge(run3_counts_with_pool,run3[,c(2,3)])

## sort sample by count
run3_counts_with_pool$sample <- factor(run3_counts_with_pool$sample, levels=unique(as.character(run3_counts_with_pool[rev(order(run3_counts_with_pool$sum_run3)),'sample'])) )

ggplot(run3_counts_with_pool, aes(x=sample, y=sum_run3, colour=run3_counts_with_pool$'pool-lane')) + geom_point() + xlab('Sample')+ ylab('read_counts') +opts(title='Tarin_1 run 2012-08-28') + opts(axis.text.x=theme_text(angle=90, hjust=1)) + labs(colour='Pools') + geom_hline(yintercept=c(10000000, 7000000))

## make boxplots
ggplot()

# sort smaple by index
run3_counts_with_pool$sample <- factor(run3_counts_with_pool$sample, levels=unique(as.character(run3_counts_with_pool[rev(order(run3_counts_with_pool$index)),'sample'])) )

ggplot(run3_counts_with_pool, aes(x=sample, y=sum_run3, colour=run3_counts_with_pool$index)) + geom_point() + xlab('Sample')+ ylab('read_counts') +opts(title='Tarin_1 run 2012-08-28') + opts(axis.text.x=theme_text(angle=90, hjust=1)) + labs(colour='index') + geom_hline(yintercept=c(10000000, 7000000))



## how many samples with reads less than 10 milion
all_counts_with_pool$'pool-lane' <- ifelse(test=is.na(all_counts_with_pool$'pool-lane'), yes='previous_pool', no=as.character(all_counts_with_pool$'pool-lane'))
table(subset(all_counts_with_pool, total<10000000)$'pool-lane')

## TARIN2-Pool1 TARIN2-Pool2 TARIN2-Pool3 TARIN2-Pool4 
##            7           10            9           24 

## list of samples/pool and how many reds they need for 100000

lt_10m <- subset(all_counts_with_pool, total<10000000)
lt_10m.sorted <- lt_10m[order(lt_10m$'pool-lane', lt_10m$total), c("pool-lane", "sample","total")]

write.table(lt_10m.sorted, file=paste('tarin1-samples_less_10_million-after_run',tarin1.last_run,'_',date,'.tab', sep=''),  row.names=FALSE, sep="\t", dec=',', quote=FALSE)



needed_reads.sorted <- data.frame(lt_10m.sorted$'pool-lane', lt_10m.sorted$sample, -1*(lt_10m.sorted$total - 10000000)/1000000)
colnames(needed_reads.sorted) = c("pool","sample","missed_reads")


write.table(needed_reads.sorted, file=paste('tarin1-number_of_reads_still_needed_by_sample-after_run',tarin1.last_run,'_',date,'.tab', sep=''),  row.names=FALSE, sep="\t", dec=',', quote=FALSE)

