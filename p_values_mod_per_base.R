#!/usr/local/bin/Rscript
#Read in test statistics that are used to get pvalues 
statistic<-read.table("test_stats_mod_per_base.csv",header=TRUE, sep=",", dec=".") 

#Naming each column to analyze separately
A_f= statistic [,1]
A_m= statistic [,2]
A_non= statistic [,3]

U_f= statistic [,1]
U_m= statistic [,2]
U_non= statistic [,3]

C_f= statistic [,1]
C_m= statistic [,2]
C_non= statistic [,3]

G_f= statistic [,1]
G_m= statistic [,2]
G_non= statistic [,3]

#Getting pvlaues for each base with degrees of freedom of 238 
pvalue_A_f=2*pt(-abs(A_f),238)
pvalue_A_m=2*pt(-abs(A_m),238)
pvalue_A_non=2*pt(-abs(A_non),238)

pvalue_U_f=2*pt(-abs(U_f),238)
pvalue_U_m=2*pt(-abs(U_m),238)
pvalue_U_non=2*pt(-abs(U_non),238)

pvalue_C_f=2*pt(-abs(C_f),238)
pvalue_C_m=2*pt(-abs(C_m),238)
pvalue_C_non=2*pt(-abs(C_non),238)

pvalue_G_f=2*pt(-abs(G_f),238)
pvalue_G_m=2*pt(-abs(G_m),238)
pvalue_G_non=2*pt(-abs(G_non),238)

#Combining all pvalues into same data frame in order to export together 
finalpvalue<-data.frame(pvalue_A_f,pvalue_A_m,pvalue_A_non,pvalue_U_f,pvalue_U_m,pvalue_U_non,pvalue_C_f,pvalue_C_m,pvalue_C_non,pvalue_G_f,pvalue_G_m,pvalue_G_non)

#Export data frame of pvalues into csv file
write.table(finalpvalue, file="pvalue_mod_per_base.csv",sep=",",row.names=F)
