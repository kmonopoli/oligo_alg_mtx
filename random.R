#!/usr/local/bin/Rscript
match<-NULL
x=c("A","U","G","C")  ##  note this C and G are flipped vs weight matrix!
p=c(0.2861271676300578, 0.28034682080924855, 0.2398843930635838, 0.1936416184971098)
#p=c(0.24277456647398843, 0.3063583815028902, 0.2138728323699422, 0.23699421965317918)
#p=c(0.325144509,0.30867052,0.185982659,0.180202312)
for(i in 1:150)
  {
#    random=sample(x,size=346*1, replace=TRUE, prob=p) 
    random=sample(x,size=150*1, replace=TRUE, prob=p) 
    total<-sum(random=="A")
    match<-c( match, total)
  }
  
  mean(match) 
  median(match) 
  sd(match)


