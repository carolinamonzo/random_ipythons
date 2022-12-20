##
###
##



# FUNCTIONS ---------------------------------------------------------------
# lane is the first arg in order to use sapply
get_total_reads_per_lane <- function (lane,dat){sum(dat[dat$lane==lane,"total_sum"])}
get_percent <- function(lane,val){as.vector(100*val/totals[lane])}


############################################################################

library ('ggplot2')
require('plyr')
date <-Sys.Date()

## lane is the first arg in order to use sapply
# get_total_reads_per_lane <- function (lane,dat){sum(dat[dat$lane==lane,"sum"])}
# get_percent <- function(lane,val){as.vector(100*val/totals[lane])}

flowcell <- "130417_H114_0273_AD21MKACXX"  ##  <<====== var
flowcell.dir <-"130417"  ###<<<=== var
mismatch=1
## read_counts-130417_H114_0273_AD21MKACXX-0_mismatch.tab
counts_file <- paste0("read_counts-", flowcell,'-',mismatch,'_mismatch.tab')

#dir <- "/home/anab/javier_labo/workspace/exome"
dir <- paste0("/home/anab/javier_labo/illumina/illumina_runs/",flowcell.dir,"/fastq_counts")  ## <<==== var
data<- read.delim(file.path(dir,counts_file))
str(data)

unique(data$lane)

## summarize the samples (or multiply by 2) because each sample has two read_groups in each lane
#ggplot(data=data)+geom_point(aes(x=sample, y=count*2))+facet_wrap(~ lane, scales='free_x')+theme(axis.text.x  = element_text(angle=90, vjust=0.5))
data_sum <- ddply(data, .(project,sample,lane),summarize, read_count=sum(count))
str(data_sum)

wanted_lanes <- paste("L00", c(1,2,3,4,5,7,8),sep='')
ggplot(data=data_sum)+geom_point(aes(x=sample, y=read_count))+facet_wrap(~ lane, scales='free_x')+theme(axis.text.x  = element_text(angle=90, vjust=0.5))
fc<- unique(as.vector(data$flowcell))
fc
plot.path<-file.path(dir,paste0("counts_reads_",fc,'-lanes_',paste(wanted_lanes,collapse='-'),'-',date,".pdf"))
plot.path
ggsave(plot.path)

totals<-sapply(wanted_lanes, get_total_reads_per_lane, data_sum)

# add percent
data_sum <- ddply(data_sum,.(lane),transform,percent=round(x=100*read_count/sum(read_count),2))


######
## only some lanes
######

wanted_lanes <- c("L004", "L005", "L006", "L007", "L008")  ## <<<=== var 
corunya.dat <- data_sum[data_sum$lane %in% wanted_lanes,]
str(corunya.dat)
#############################
## subset for only some lanes
#############################
totals<-sapply(wanted_lanes, get_total_reads_per_lane, corunya.dat)
corunya.dat <- ddply(corunya.dat,.(lane),transform,percent=round(x=100*read_count/sum(read_count),2))

## total counts per lane
sum_counts <- ddply(corunya.dat, 
                    .(lane),
                    summarize, 
                    sum_total=sum(read_count), 
                    mean_sum_total=round(x=sum(read_count/length(sample)),2)
                    )
sum_counts
# lane sum_total mean_sum_total
# 1 L001 137067122       11422260
# 2 L002 149015970       12417998
# 3 L003 153062412       12755201


## get counts
ggplot(data=corunya.dat)+geom_point(aes(x=sample, y=read_count))+
  facet_wrap(~ lane, scales='free_x')+
  geom_hline(aes(yintercept=mean_sum_total),sum_counts) +
  theme(axis.text.x  = element_text(angle=90, vjust=0.5))
fc<- unique(as.vector(data$flowcell))
fc
plot.path.pdf<-file.path(dir,paste0("count_reads_",fc,'-lanes_',paste(wanted_lanes,collapse='-'),'-',date,".pdf"))
plot.path.pdf
plot.path.png<-file.path(dir,paste0("count_reads_",fc,'-lanes_',paste(wanted_lanes,collapse='-'),'-',date,".png"))
plot.path.png

ggsave(plot.path.pdf)
ggsave(plot.path.png)

## percentage
ggplot(data=corunya.dat)+geom_point(aes(x=sample, y=percent)) +
  facet_wrap(~ lane, scales='free_x') +
  geom_hline(aes(yintercept=100/6))+
  theme(axis.text.x  = element_text(angle=90, vjust=0.5))
fc<- unique(as.vector(data$flowcell))
fc
plot.path.pdf<-file.path(dir,paste0("percent_reads_",fc,'-lanes_',paste(wanted_lanes,collapse='-'),'-',date,".pdf"))
plot.path.pdf
plot.path.png<-file.path(dir,paste0("percent_reads_",fc,'-lanes_',paste(wanted_lanes,collapse='-'),'-',date,".png"))
plot.path.png

ggsave(plot.path.pdf)
ggsave(plot.path.png)


#<TO_DO> do with gsub
table.path <- file.path(dir,paste0("percent_reads_",fc,'-lanes_',paste(wanted_lanes,collapse='-'),'-',date,".tab"))
table.path
write.table(x=corunya.dat,file=table.path,quote=FALSE,row.names=FALSE)


## plot the data with the l1 and 8 combined


# process by projects (conbine lanes) -------------------------------------



#####
## add l1 nad l8 into HMC1

l1_and_l8 <- ddply(subset(corunya.dat, lane %in% c("L001", "L008")), .(sample), summarize, read_count=sum(read_count))

l1_and_l8$project <- rep("HMC1", nrow(l1_and_l8))
l1_and_l8
l1_and_l8 <- ddply (l1_and_l8, .(project), transform, percent=round(x=100*read_count/sum(read_count),2))[,c("project", "sample", "read_count", "percent")]
l3_l4 <- subset(corunya.dat, lane%in% c("L003", "L004"))[,c(1,2,4,5)]
l3_l4

hmc.dat <- rbind(l1_and_l8, l3_l4)
hmc.dat

# <to_do> add title
ggplot(data=hmc.dat)+geom_point(aes(x=sample, y=percent))+facet_wrap(~ project, scales='free_x')+theme(axis.text.x  = element_text(angle=90, vjust=0.5))
plot.path<-file.path(dir,paste0("percent_reads",'-projects_hmc1,2,3','-',date,".pdf"))
plot.path
ggsave(plot.path)

ggplot(data=hmc.dat)+geom_point(aes(x=sample, y=read_count))+facet_wrap(~ project, scales='free_x')+theme(axis.text.x  = element_text(angle=90, vjust=0.5))
plot.path<-file.path(dir,paste0("read_counts",'-projects_hmc1,2,3','-',date,".pdf"))
plot.path
ggsave(plot.path)



## seed(1234)
## a <- rep(letters[seq( from = 1, to = 10 )],2)
## b <- as.factor(rep(c(rep(1,10),rep(2,10)),2))
## b
## c <- sort(rnorm(mean=10,sd=0.7,n=40))
## d <- c(rep("set1",20), rep("set2", 20))
## test <- data.frame(a,b,c,d)
## test
## str(test)
## ggplot(data=test)+geom_point(aes(x=a, y=c))+facet_wrap(~ d, scales='free_x')+theme(axis.text.x  = element_text(angle=90))

## ******************************************************
##   plot the 0 mismatchs counts per lane and the undetermined
## ******************************************************

## run 130417


