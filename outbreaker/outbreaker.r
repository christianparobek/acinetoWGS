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
library(stringr)
library(igraph)

################################################
######### DEFINE SOME USEFUL FUNCTIONS #########
################################################

name.cleaner <- function(x){
  str_extract(x, "A[0-9]*")
} # a function to clean up the sample names

detectR <- function(x){
  neighbors(g, x, mode="out")
} # a function to calc R (infective) value for each samp

################################################
########### READ IN THE DATA WE NEED ###########
################################################

good42 <- fasta2DNAbin(file="good42_UG_pass.fa")
good42_dates <- as.Date(c("2008/3/20", "2008/5/6", "2007/10/7", "2008/3/17", "2008/2/24", "2008/3/8", "2008/3/31", "2008/3/10", "2007/12/20", "2008/3/4", "2008/6/17", "2008/2/9", "2008/7/25", "2008/4/14", "2008/12/8", "2008/11/27", "2009/2/8", "2009/4/13", "2009/3/23", "2009/5/25", "2009/2/17", "2009/4/13", "2009/4/7", "2009/8/18", "2009/9/22", "2009/9/8", "2009/9/8", "2009/8/25", "2009/9/1", "2010/7/20", "2009/11/30", "2009/6/17", "2009/9/15", "2010/5/18", "2010/5/18", "2010/5/25", "2010/2/8", "2010/6/2", "2009/12/15", "2010/1/5", "2010/2/3", "2010/3/23"))

names <- name.cleaner(rownames(good42))

################################################
############ GET DISTRIBUTION OF W #############
################################################

x <- seq(0, 360, length = 361) # 361 days

#### EXPONENTIAL DECAY ####
  # Used by Stoesser et al (2014) for their NDM Kleb paper
  expDecay <- 200*exp(-0.3*x)
  plot(expDecay, type="l", lwd=2)

#### NORMAL DISTRIBUTION ####
  # Recommended by Jombart
  normDist <- dnorm(x, mean=20, sd=10) # mean 50, very wide SD
  plot(x, normDist, type="l", lwd=2)

#### LOG-NORMAL DISTRIBUTION ####
  # A variable can be modeled as log-norm if it can be thought of as the multiplicative 
  # product of many independent random variables each of which is positive.
  # The log-norm dist makes sense clinically
  sigma = 0.8
  mu = 3
  logNorm = dlnorm(x, mu, sigma)
  plot(x, logNorm, type = "l", xlim=c(0,50))

################################################
########### RUN OUTBREAKER FUNCTION ############
################################################

## Run outbreaker
resp <- outbreaker.parallel(n.runs=4, parallel=TRUE, n.iter=400000, dna=good42, dates=good42_dates, w.dens=logNorm, max.kappa=50)
  # w.trunc is the nubmer of generations at which transmission potential falls to zero
  # max.kappa is an integer indicating the maximum number of generations between a case and its most recent sampled ancestor; defaults to 10.
  # only when i messed around with these values did we get something that worked...

## Do chains converge?
plotChains(resp)
plotChains(resp, burnin=2e4)

## Look at infective time
plotOutbreak(resp, burnin=399999, lwd.arrow=0, annot="", axes=FALSE, xlab="Days Since Index Case Presentation", col.pal=colorRampPalette(c("gray50","gray50")))
axis(1, at=c(13800,14000,14200,14400,14600,14800,15000,15200), labels=c(0,200,400,600,800,1000,1200,1400))
axis(2, at=1:42, labels=names, las=2, cex.axis=0.9)

######################################
##### DRAW THE TRANSMISSION TREE #####
######################################

g <- transGraph(resp, thresh=0.3, annot="") # draw tx tree

R <- lapply(V(g), detectR) # get nb "out" transmissions per sample
case_size <- 10+unlist(lapply(R, length))*3 # make the sizing vector
Tinf <- resp$chains[resp$chains$step>=2e4, grep("Tinf", names(resp$chains))] # color on date
case.color <- any2col(apply(Tinf,2,mean), col.pal=spectral) # color on date

V(g)$color <- case.color$col # assign color
V(g)$size <- case_size # assign size
V(g)$label <- names # assign name

tkplot(g, canvas.width=600, canvas.height=600) # adjust to my liking # layout=coords, 
coords <- tkplot.getcoords(2) # get the coordinates

par(mar=c(0,0,0,0))
plot(g, layout=coords, vertex.label.color="black", edge.arrow.size=0.70) # make final plot

# Add the legend
## Plotting the colors for infection dates
par(mar=c(1,4,6,10), mfrow=c(1,2))
my.colors = colorRampPalette(c("#3288BD", "#479FB3", "#4BA4B1", "#64C0A5", 
                               "#7AC9A4", "#8AD0A4", "#96D4A4", "#A1D9A4", 
                               "#CDEB9C", "#D2ED9C", "#E9F69D", "#EFF8A7", 
                               "#F7FCB3", "#FBFDBA", "#FEF8B4", "#FEE99B", 
                               "#FEE797", "#FDB768", "#FCAA5F", "#F99555", 
                               "#F7854E", "#F06744", "#F06744", "#E1504A", "#D53E4F"))
z=matrix(1:100,nrow=1)
x=1
y=seq(0,96,len=100) # supposing 3 and 2345 are the range of your data
image(x,y,z,col=my.colors(100),axes=FALSE,xlab="",ylab="")
axis(4, at=(0:24)*4, labels=c("2010 - Jul", "", "2010 - May", "", "2010 - Feb", 
                              "", "2009 - Dec", "", "2009 - Sep", "", 
                              "2009 - Jun", "", "2009 - Apr", "", "2009 - Feb", 
                              "", "2008 - Nov", "", "2008 - Jun", "", 
                              "2008 - Apr", "", "2008 - Feb", "", "2007 - Oct"), 
     line=0.5, lwd=2, las=1, cex.axis=1.25)
mtext("Date of\nInfection", side=3, at=2.4, line=1, cex=1.5, font=2)


my.colors = colorRampPalette(c("white", "black"))
z=matrix(1:100,nrow=1)
x=1
y=seq(0,1,len=100) # supposing 3 and 2345 are the range of your data
image(x,y,z,col=my.colors(100),axes=FALSE,xlab="",ylab="")
axis(4, at=c(0,0.5,1), 
     line=0.5, lwd=2, las=1, cex.axis=1.25)
mtext("Posterior\nSupport", side=3, at=1.5, line=1, cex=1.5, font=2)

# then save the graph as an SVG etc, with big dimensions so that the font fits inside the bubbles


###############################
######## MUTATION RATE ########
###############################

## Determine and plot mutation rates
## I believe the get.mu() function expresses mutations in terms of unit time, in our case this is days
## Also, for mutation rate, we need to analyze the two outbreak strains separately
## Analyzing them together gave us a mutation rate that was an order of magnitude too high :)

## Read in the Yellow and Green datasets
yellow <- fasta2DNAbin(file="yellow.fa")
yellow_dates <- as.Date(c("2008/3/20", "2008/5/6", "2007/10/7", "2008/3/17", "2008/2/24", 
                          "2008/3/8", "2008/3/31", "2008/3/10", "2007/12/20", "2008/3/4", 
                          "2008/2/9", "2008/7/25", "2008/4/14")) 

green <- fasta2DNAbin(file="green.fa")
green_dates <- as.Date(c("2008/12/8", "2008/11/27", "2009/2/8", "2009/4/13", "2009/3/23", 
                         "2009/5/25", "2009/4/13", "2009/4/7", "2009/8/18", "2009/9/22", 
                         "2009/9/8", "2009/9/8", "2009/9/1", "2010/7/20", "2009/11/30", 
                         "2009/9/15", "2010/5/18", "2010/5/18", "2010/5/25", "2010/2/8", 
                         "2009/12/15", "2010/1/5", "2010/2/3", "2010/3/23"))                     

## Run outbreaker on each dataset
yellow_res <- outbreaker.parallel(n.runs=4, parallel=TRUE, n.iter=100000, dna=yellow, dates=yellow_dates, w.dens=logNorm, max.kappa=50)
green_res <- outbreaker.parallel(n.runs=4, parallel=TRUE, n.iter=100000, dna=green, dates=green_dates, w.dens=logNorm, max.kappa=50)

# Change so that mutation rate is expressed in mutations/year
genome_size <- 4461520
yellow_mut_per_year <- density(yellow_mu)$x*genome_size*365
green_mut_per_year <- density(green_mu)$x*genome_size*365

# Plot graph of the mutation rate
par(xpd=NA, mgp = c(3.75, 2.75, 2))
plot(green_mut_per_year, density(green_mu)$y, xlim=c(0,150), type="none", axes=FALSE, xlab="Mutations per Genome per Year", ylab="")
lines(green_mut_per_year, density(green_mu)$y, lwd=3, col="black")
lines(yellow_mut_per_year, density(yellow_mu)$y, lwd=3, col="gray50")
axis(1, at=c(0,25,50,75,100,125,150), labels=c("0","","50","","100","","150"))
axis(2, labels=FALSE, at=c(0, 260000000))
mtext("Probability Density", side=2, line=2.5)

# Add confidence intervals
segments(quantile(green_mu, 0.025)*genome_size*365, 
         -15000000, 
         quantile(green_mu, 0.975)*genome_size*365, 
         -15000000, 
         col="black", 
         lwd=3,
         lty=3)
points(quantile(green_mu, 0.5)*genome_size*365, -15000000, pch=19, col="black", cex=1.5)
segments(quantile(yellow_mu, 0.025)*genome_size*365, 
         -25000000, 
         quantile(yellow_mu, 0.975)*genome_size*365, 
         -25000000, 
         col="gray50", 
         lwd=3,
         lty=3)
points(quantile(yellow_mu, 0.5)*genome_size*365, -25000000, pch=19, col="gray50", cex=1.5)

# Add a legend
legend(55, 90000000, 
       legend=c("Probability Density", "Mean Mutation Rate", "95% Conf. Interval"), 
       lty=c(1, NA, 3), 
       pch=c(NA, 19, NA), 
       lwd=2, col=c("gray25"),
       pt.cex=c(NA, 1.5, NA),
       bty="n")




## Transmission Tree Work
get.tTree(resp, burnin=2e4, best=c("ancestries","tree"))
plot(resp, y=NULL, edge.col="black", col.edge.by="prob", col.pal=NULL, annot=c("dist","n.gen","prob"), sep="/")
as.igraph(res, edge.col="black", col.edge.by="prob", col.pal=NULL, annot=c("dist","n.gen","prob"), sep="/")
plot(resp, annot="dist", main="Data - transmission tree")

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


## VISUALIZE DYNAMICS
matplot(fakeOutbreak$dat$dynam, type="o", pch=20, lty=1,
        main="Outbreak dynamics", xlim=c(0,28))
legend("topright", legend=c("S","I","R"), lty=1, col=1:3)
graphics.off()

