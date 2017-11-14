#!/usr/local/bin/Rscript
#Setting the location to read in test statistics 
#setwd("C:/Users/drhomsy/Desktop/Diana")
#Read in test statistics that are used to get pvalues 
statistic<-read.table("test_stats_mod.csv",header=TRUE, sep=",", dec=".") 
#Naming each column to analyze separately
f= statistic [,1]
m= statistic [,2]
non= statistic [,3]
#Getting pvlaues for each base with degrees of freedom of 188 
pvalue_f=2*pt(-abs(f),238)#188)
pvalue_m=2*pt(-abs(m),238)#188)
pvalue_non=2*pt(-abs(non),238)#188)
#Combining all pvalues into same data frame in order to export together 
finalpvalue<-data.frame(pvalue_f,pvalue_m,pvalue_non)
#Export data frame of pvalues into csv file
write.table(finalpvalue, file="pvalue_mod.csv",sep=",",row.names=F)
