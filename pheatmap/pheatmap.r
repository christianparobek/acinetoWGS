## An ape script to analyze Acineto WGS Data
## Started March 25, 2015
## Example: http://stackoverflow.com/questions/15153202/how-to-create-a-heatmap-with-a-fixed-external-hierarchical-cluster

################################################
########## LOAD THE LIBRARIES WE NEED ##########
################################################

library(ape)
library(phangorn)
library(adegenet)

################################################
######### DEFINE SOME USEFUL FUNCTIONS #########
################################################

name.cleaner <- function(x){
  library(stringr)
  str_extract(x, "A[0-9]*")
} # a function to clean up the sample names

a03.comparer <- function(x){
  sum(as.numeric(!mapply(identical, x, ord_mat[3,])))
} # a function to count snp diffs bt A03 and other indivs

diff.paster <- function(x){
  paste(rownames(ord_mat)[x], " (", snps[x], ")", sep="")
} # add the num of snp diffs to the names

################################################
################# READ IN DATA #################
################################################

good42_fasta <- read.FASTA("good42_UG_pass.fa") # read in fasta for tree
names(good42_fasta) <- name.cleaner(names(good42_fasta)) # change names

good42_df <- 
  read.table(
    "/run/user/1001/gvfs/sftp:host=kure.unc.edu,user=prchrist/proj/julianog/users/ChristianP/acinetoWGS/adegenet/good42_UG_pass.str", 
    skip=1) # read in data for the genlight object (for the heatmap)
ind_names <- good42_df$V1 # get the individual names
good42_df <- good42_df[-c(1,2)] # remove ind and pop columns from data frame
good42_gl <- new("genlight", good42_df) # convert df into genlight object
indNames(good42_gl) <- name.cleaner(ind_names) # get, clean, and add the ind names
ploidy(good42_gl) <- 1 # add back population information

################################################
################ CREATE THE TREE ###############
################################################

dm <- dist.dna(good42_fasta) # calcualte distance matrix
treeNJ <- NJ(dm) # makes a phylo object
plot(treeNJ, "phy", main="NJ") # Plot the tree
rooted <- root(treeNJ, outgroup="A42", resolve.root = TRUE) # root it
ultra <- compute.brlen(rooted, method="Grafen") # ultrametricize
plot(ultra, "phy", main="NJ") # Plot the tree again

hc <- as.hclust(ultra) # convert phylo to dendrogram
dend <- as.dendrogram(hc) # convert phylo to dendrogram

################################################
###### COUNT DIFFS BETWEEN A03 AND OTHERS ######
################################################

good42_matrix <- as.matrix(good42_gl) # convert gl to matrix
ord_mat <- good42_matrix[ultra$tip,] # ord matrix by phylo tree

snps <- apply(ord_mat, 1, a03.comparer) # calculate nuc diffs bt A03 and others
new_names <- unlist(lapply(1:nrow(ord_mat), diff.paster)) # make list of updated names
rownames(ord_mat) <- new_names # assign new names to ord_mat

################################################
############### MAKE THE HEATMAP ###############
################################################

heatmap(ord_mat, Rowv=dend, Colv=NA, labCol="", col=topo.colors(4), margins=c(1,20)) # draw pheatmap
text(0.29, -8, "35551 SNVs")

# ML example: http://cran.r-project.org/web/packages/phangorn/vignettes/Trees.pdf
# Need to run the treeNJ phylo object through the ML Process
# Use the above link for an example of how to do this
fit <- pml(treeNJ, data=as.phyDat(data))
