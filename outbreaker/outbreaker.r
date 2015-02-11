## An Outbreaker Script to analyze Acineto WGS Data
## Started December 22, 2014
## To be run on Kure
## In the Acinetobacter/outbreaker file
## Launch R using `R`
## Make sure it's running on an interactive node - NOT the login node
## Intro to outbreaker: https://a4948e14-a-62cb3a1a-s-sites.googlegroups.com/site/therepiproject/outbreaker-intro.pdf?attachauth=ANoY7cq4FB8ojIwxPWn-PY-T2JYjpo8U-HmU7gQoxB-NagMJxmWUk85QhSgnyOEjVqkWwwMntEJQnt3P9wqJvLCwlCrsuuYlIHscbtOQit8WZnI-YFijVPnEDEy7n_02oEZZa3u6-A8a7syUTgNkLJCWTkna2_cEmDReJpBEK3Pu-ZgEIDkN8H8qlcHFSxso0-Hg7E6pLYq29r-78zM-EbHvDAYdzRRYIUp49vjhg20tC4SvW_j5Q3s%3D&attredirects=0
## PDF of functions: http://cran.r-project.org/web/packages/outbreaker/outbreaker.pdf

################################################
########## LOAD THE LIBRARIES WE NEED ##########
################################################

library(ape)
library(adegenet)
library(outbreaker)

################################################
########### READ IN THE DATA WE NEED ###########
################################################

DNAbin <- fasta2DNAbin(file="good42.fa") # a multifasta

dates <- c("2008/3/20", "2008/5/6", "2007/10/7", "2008/3/17", "2008/2/24", "2008/3/8", "2008/3/31", "2008/3/10", "2007/12/20", "2008/3/4", "2008/6/17", "2008/2/9", "2008/7/25", "2008/4/14", "2008/12/8", "2008/11/27", "2009/2/8", "2009/4/13", "2009/3/23", "2009/5/25", "2009/2/17", "2009/4/13", "2009/4/7", "2009/8/18", "2009/9/22", "2009/9/8", "2009/9/8", "2009/8/25", "2009/9/1", "2010/7/20", "2009/11/30", "2009/6/17", "2009/9/15", "2010/5/18", "2010/5/18", "2010/5/25", "2010/2/8", "2010/6/2", "2009/12/15", "2010/1/5", "2010/2/3", "2010/3/23")
dates2 <- as.Date(dates) # convert to outbreaker's date format

################################################
############ GET DISTRIBUTION OF W #############
################################################

#### EXPONENTIAL DECAY ####
  # Used by Stoesser et al (2014) for their NDM Kleb paper

  x <- seq(0, 120, length = 121) # 101 days

  expDecay <- 200*exp(-0.3*x)
  plot(expDecay, type="l", lwd=2)

#### NORMAL DISTRIBUTION ####
  # Recommended by Jombart

  normDist <- dnorm(x, mean=20, sd=10) # mean 50, very wide SD
  plot(x, normDist, type="l", lwd=2)

#### LOG-NORMAL DISTRIBUTION ####
  # A variable might be modeled as log-normal 
  # if it can be thought of as the multiplicative 
  # product of many independent random variables 
  # each of which is positive.

  # Jombart also recommended modeling intervals
  # between successive cases.

  # The log-norm dist makes sense clinically

  sigma = 0.6
  mu = 3
  logNorm = dlnorm(x, mu, sigma)
  plot(x, logNorm, type = "l")

################################################
########### DO THE OUTBREAKER THING ############
################################################

res <- outbreaker(dna=DNAbin, dates=dates2, w.dens=expDecay, w.trunc=120, max.kappa=100, import.method='gen', init.tree='seqTrack') 
  # w.trunc is the nubmer of generations at which transmission potention falls to zero
  # max.kappa is an integer indicating the maximum number of generations between a case and its most recent sampled ancestor; defaults to 10.
  # only when i messed around with these values did we get something that worked...
  # do we need to specify an import.method?
  # needs more optimizing!!!

g <- transGraph(res, thres=0.08, annot="", edge.curved=FALSE)









## Transmission Tree Work
get.tTree(res, burnin=2e4, best=c("ancestries","tree"))
plot(res, y=NULL, edge.col="black", col.edge.by="prob", col.pal=NULL, annot=c("dist","n.gen","prob"), sep="/")
as.igraph(res, edge.col="black", col.edge.by="prob", col.pal=NULL, annot=c("dist","n.gen","prob"), sep="/")
plot(res, annot="dist", main="Data - transmission tree")

## GET CONSENSUS ANCESTRIES
tre <- get.tTree(res)
plot(tre, annot="", main="Consensus ancestries")


## Determine histogram of dates
## Read in the acinetoTimeDiffs.txt file
## Then run the following command
hist(acinetoTimeDiffs$V1, breaks=40)


## Not run:
## COMMAND LINES TO GENERATE SIMILAR DATA ##
#w <- c(0, 0.5, 1, 0.75)
w <- c(0, 0.25, 0.5, 0.25)
## note: this works only if outbreak has at least 30 case
dat <- simOutbreak(R0 = 2, infec.curve = w, n.hosts = 100)[1:30]
collecDates <- dat$onset + sample(0:3, size=30, replace=TRUE, prob=w)


## LOAD TOY DATASET
load("fakeOutbreak.RData")
## OR
#data(fakeOutbreak)

## MAKE OBJECTS WITHIN DATASET SEARCHABLE
attach(fakeOutbreak)

## VISUALIZE DYNAMICS
matplot(fakeOutbreak$dat$dynam, type="o", pch=20, lty=1,
        main="Outbreak dynamics", xlim=c(0,28))
legend("topright", legend=c("S","I","R"), lty=1, col=1:3)
graphics.off()

## VISUALIZE TRANSMISSION TREE
plot(dat, annot="dist", main="Data - transmission tree")
mtext(side=3, "arrow annotations are numbers of mutations")
graphics.off()

## examine MCMC
plotChains(res)
plotChains(res,type="dens")
plotChains(res,type="dens", what="mu1", burnin=2e4)
graphics.off()


## DISTRIBUTION OF INFECTION DATES
barplot(w, main="Generation time distribution", ylab="probability", xlab="days", names=0:3)
graphics.off()


## VISUALIZE RECONSTRUCTED TRANSMISSION TREES
library(igraph)
library(adegenet)
transGraph(res, thres=0.5, annot="")
graphics.off()




