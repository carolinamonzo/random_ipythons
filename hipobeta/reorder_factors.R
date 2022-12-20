##
### reorder factors
##


library('ggplot2')
meanprice <- tapply(diamonds$price, diamonds$cut, mean)
# 
# > meanprice
#      Fair      Good Very Good   Premium     Ideal 
#  4358.758  3928.864  3981.760  4584.258  3457.542 

cut <- factor(levels(diamonds$cut), levels = levels(diamonds$cut))

## unsorted
qplot(cut, meanprice, geom="bar", stat="identity")

## sorted by meanprice
qplot(reorder(cut, meanprice), meanprice, geom="bar", stat="identity")
