#!/usr/local/bin/Rscript
#Setting the location to read in test statistics 
#setwd("C:/Users/drhomsy/Desktop/Diana")
#Read in test statistics that are used to get pvalues 
statistic<-read.table("test_stats.csv",header=TRUE, sep=",", dec=".") 
#Naming each column to analyze separately
A= statistic [,1]
U= statistic [,2]
C= statistic [,3]
G= statistic [,4]
#G= statistic [,3]
#C= statistic [,4]
#Getting pvlaues for each base with degrees of freedom of 188 
pvalueA=2*pt(-abs(A),161)#238)#182)#154)#136)#188)#318)#188)
pvalueU=2*pt(-abs(U),161)#238)#182)#154)#136)#188)#318)#188)
pvalueC=2*pt(-abs(C),161)#238)#182)#154)#136)#188)#318)#188)
pvalueG=2*pt(-abs(G),161)#238)#182)#154)#136)#188)#318)#188)                               
#Combining all pvalues into same data frame in order to export together 
finalpvalue<-data.frame(pvalueA,pvalueU,pvalueC,pvalueG)
#Export data frame of pvalues into csv file
write.table(finalpvalue, file="pvalue.csv",sep=",",row.names=F)
