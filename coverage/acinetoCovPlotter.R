## For plotting Acinetobacter % Genome Coverage
## First run `bedtools genomecov` on Kure
## Then add names in Excel and import into R

### READ IN THE DATA AND FORMAT
genome10 <- read.table("coverage10.txt", header=FALSE)
genome10 <- genome10[with(genome10, order(V3)),]
genome10$V4 <- 1:48
genome11 <- genome10[-(1:2),] # remove the S. maltophilia


### PLOT IT

svg("coverage.svg", width = 7, height = 4.2)
tiff("coverage.tiff", width = 7, height = 4.2, units = "in", res = 300, compression = "lzw")

plot(genome11$V2[1:46] ~ genome11$V4[1:46], 
     axes=FALSE, 
     xlab="Sample ID", 
     ylab="Genome Fraction", 
     main="Short-read Coverage over A03 Long-Read Assembly",
     ylim=c(0,1), 
     col="black",
     cex.main = 1.1,
     type = "n")
axis(1, at=3:48, labels=genome11$V1, cex.axis=0.75, las=2)
axis(2, at=c(0.0,0.2,0.4,0.6,0.8,1.0), las=1)
#rect(1, 0, 2, 1, border = "white", col = "grey")

points(genome11$V2[5:46] ~ genome11$V4[5:46], pch=1) # add 10x for good Ab
points(genome11$V3[5:46] ~ genome11$V4[5:46], pch=19) # add 5x for good Ab
points(genome11$V2[1:4] ~ genome11$V4[1:4], pch=5) # add 10x for low cov Ab
points(genome11$V3[1:4] ~ genome11$V4[1:4], pch=18, cex = 1.5) # add 5x for low cov Ab

legend(30, 0.3, 
       legend = c(expression(Coverage >= 05*x), expression(Coverage >= 10*x)), 
                  pch = c(15, 0), bty = "n", pt.cex = 1.3)
dev.off()