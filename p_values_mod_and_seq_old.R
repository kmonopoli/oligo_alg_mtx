#!/usr/local/bin/Rscript
#Read in test statistics that are used to get pvalues 
statistic<-read.table("test_stats_mod_and_seq.csv",header=TRUE, sep=",", dec=".") 

#Naming each column to analyze separately
a= statistic [,1]
u= statistic [,2]
c= statistic [,3]
g= statistic [,4]

pyr_f= statistic [,5]
pyr_m= statistic [,6]
pyr_non= statistic [,7]

pur_f= statistic [,8]
pur_m= statistic [,9]
pur_non= statistic [,10]

#Getting pvlaues for each base with degrees of freedom of 188 
pvalue_a=2*pt(-abs(a),188)
pvalue_u=2*pt(-abs(u),188)
pvalue_c=2*pt(-abs(c),188)
pvalue_g=2*pt(-abs(g),188)


pvalue_pyr_f=2*pt(-abs(pyr_f),188)
pvalue_pyr_m=2*pt(-abs(pyr_m),188)
pvalue_pyr_non=2*pt(-abs(pyr_non),188)


pvalue_pur_f=2*pt(-abs(pur_f),188)
pvalue_pur_m=2*pt(-abs(pur_m),188)
pvalue_pur_non=2*pt(-abs(pur_non),188)

#Combining all pvalues into same data frame in order to export together 
finalpvalue<-data.frame(pvalue_a,pvalue_u,pvalue_c,pvalue_g,pvalue_pyr_f,pvalue_pyr_m,pvalue_pyr_non,pvalue_pur_f,pvalue_pur_m,pvalue_pur_non)

#Export data frame of pvalues into csv file
write.table(finalpvalue, file="pvalue_mod_and_seq.csv",sep=",",row.names=F)
