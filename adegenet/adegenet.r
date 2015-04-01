# For genetic analysis of my WGS Pv population(s)
# Started 10 Dec 2014
# Basics Tutorial: http://adegenet.r-forge.r-project.org/files/tutorial-basics.pdf
# Genomics Tutorial: http://adegenet.r-forge.r-project.org/files/tutorial-genomics.pdf
# Extra Commands: http://www.inside-r.org/packages/cran/adegenet/docs/.rmspaces


### LOAD LIBRARIES ###
library(adegenet)

### READ IN THE DATA ###
data <- read.table("good42_UG_pass.str", skip=1) # read in data
sorted <- data[order(data[,1]),] # sort
inds <- sorted$V1 # grab the indiv names
pops <- sorted$V2 # grab the pop names
sorted <- sorted[-c(1,2)] # remove ind and pop columns from data frame
genlight <- new("genlight", sorted) # convert data frame into genlight object
indNames(genlight) <- inds # add back individual information
ploidy(genlight) <- 1 # add back population information

### DO THE PCA CALCULATIONS ###
pca1 <- glPca(genlight) # for genlight

### PLOT PCA EIGENVALUES ###
barplot(pca1$eig, xlab = "", ylab = "Variance", main = "P. vivax Eigenvalues") # barplot

### ADD JITTER ###
pca1$scores[,1] <- jitter(pca1$scores[,1], factor=400) # add jitter if overplotting
pca1$scores[,2] <- jitter(pca1$scores[,2], factor=400) # add jitter if overplotting
pca1$scores[,3] <- jitter(pca1$scores[,3], factor=400) # add jitter if overplotting
pca1$scores[,4] <- jitter(pca1$scores[,4], factor=400) # add jitter if overplotting

### PLOT PCA PICTURE ###
plot(pca1$scores[,1], pca1$scores[,2], 
     pch=19, 
     axes=FALSE, 
     xlab=paste("PC1 - ", round(pca1$eig[1]/sum(pca1$eig)*100), "% of the Variance", sep = ""),
     ylab=paste("PC2 - ", round(pca1$eig[2]/sum(pca1$eig)*100), "% of the Variance", sep = ""),
     main="A. baumannii PCA"
)
axis(1)
axis(2)

