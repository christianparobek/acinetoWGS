## R script to plot num variants vs PFGE status
## Started 1 April 2015
## Christian P

data <- read.table("variants_vs_PFGE_results.txt", header=TRUE)

data$PFGE_status <- factor(data$PFGE_status, c("I","C","P","D"))

par(mar=c(5,6,3,3), mgp=c(3,1.75,1))
stripchart(jitter(data$Nb_variants, factor=200) ~ data$PFGE_status, 
           method="jitter", 
           vertical=TRUE, 
           pch=19, 
           cex=0.2,
           axes=FALSE,
           xlab="PFGE Status",
           ylab="Nuc. Differences (thousands)",
           font=2)
axis(1, at=c(1,2,3,4), labels=c("I","C","P","D"), lwd=2, cex.lab=2)
axis(2, las=2, lwd=2, at=c(0,10000,20000), labels=c(0,10,20))