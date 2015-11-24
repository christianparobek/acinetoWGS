## Plotting a circos plot for pairwise FST between Acinetobacter outbreaks
## Christian Parobek
## Started 22 November 2015

##########
### INITIALIZE LIBRARIES
##########
library(PopGenome)
library(circlize)
library(stringr)


##########
### DEFINE FUNCTIONS
##########

bedMaker <- function(GENOME_obj, index){
  df <- as.data.frame(cbind(GENOME_obj@region.names, GENOME_obj@nuc.F_ST.pairwise[index,])) ## make df w/chr info
  c1 <- str_extract(df$V1, "chr\\d+.vcf") ## need to clean it up to get it into bed format for circlize
  c2 <- str_replace(c1, ".vcf", "|quiver")
  chromosome <- str_replace(c2, "chr", "unitig_")
  start <- as.integer(str_extract(str_extract(df$V1, " \\d*"), "[0-9]+")) # match start value and trim whitespace
  end <- as.integer(str_extract(str_extract(df$V1, "- \\d*"), "[0-9]+")) # match end value and trim dash/whitespace
  fst_bed <- cbind.data.frame(chromosome, start, end) # make the data take bed file format
  fst_bed$fst <- as.numeric(as.character(df$V2)) # bet there's a better way to do this, but couldn't find it
  fst_bed <- fst_bed[!end < start,] ## clean entries where end < start, probably due to a bug in PopGenome
  return(fst_bed)
} ## make a bed file from the FST data to draw lines

fstDrawer <- function(region, value, ...){
  i = getI(...)
  circos.genomicPoints(region, value, col = i, ..., pch = 1, cex = 0.25, lwd = 0.01)
} ## a function to draw lines in a circos.genomicTrackPlotRegion

chrNamer <- function(region, value, ...){
  circos.genomicText(region, value, labels.column =1, facing = "clockwise")
} ## a function to draw lines in a circos.genomicTrackPlotRegion

trackPlotter <- function(bed, panel_fun){
  circos.genomicTrackPlotRegion(bed, panel.fun = panel_fun, ylim = c(0, 1), bg.col = c(rep("gray90", 14), "white"), bg.border = NA, track.height = 0.12)
} ## integrates the fstDrawer function to draw fst tracks

axisAdder <- function(){
  circos.genomicTrackPlotRegion(ylim = c(0, 1), bg.border = NA, 
                                track.height=0.05,
                                panel.fun = function(region, value, ...) {
                                  circos.axis(h = 0, major.at = NULL, labels = NULL, 
                                              labels.cex = 0.3 * par("cex"), labels.facing = "clockwise", 
                                              major.tick.percentage = 0.2)})} ## add axis without text

chrNamer <- function(){
  circos.genomicTrackPlotRegion(ylim=c(0,1),
                                bg.col=NULL, bg.border=NA, track.height=0.05,
                                panel.fun = function(region, value, ...){
                                  circos.text(mean(get.cell.meta.data("xlim")), 
                                              0.5, labels = get.cell.meta.data("sector.numeric.index"))})
} ## add chromosome names wrt sector index

leadingZeroScrubber <- function(bed){
  first_zero_index <- match(TRUE, bed$fst > 0)
  bed[-(1:first_zero_index),]
} ## For some reason there are leading zeros that get plotted outside the plot area

naScrubber <- function(bed){
  vectorToRemove <- is.na(bed$fst)
  return(bed[!vectorToRemove,])
}

##########
### READ IN DATA
##########

## read in genome data
setwd("/run/user/1000/gvfs/sftp:host=kure.unc.edu,user=prchrist/proj/julianog/users/ChristianP/acinetoWGS/circlize/")
acineto <- readData("vcf/", format="VCF")

## read in groupings
setwd("/run/user/1000/gvfs/sftp:host=kure.unc.edu,user=prchrist/proj/julianog/users/ChristianP/acinetoWGS/circlize/acineto_groups/")
outbreak1 <- read.table("outbreak1.txt")
outbreak2 <- read.table("outbreak2.txt")
outbreak3 <- read.table("outbreak3.txt")
pops <- set.populations(acineto, list(as.character(outbreak1$V1), as.character(outbreak2$V1), as.character(outbreak3$V1)))

## read chromosome info for circlize
setwd("/run/user/1000/gvfs/sftp:host=kure.unc.edu,user=prchrist/proj/julianog/users/ChristianP/acinetoWGS/circlize/")
chr_info <- read.table("acineto_chr_data.txt", header=TRUE)
chr_info$name <- factor(chr_info$name, levels = c("Ab_names",
                                                  "unitig_35|quiver",
                                                  "unitig_3|quiver",
                                                  "unitig_1|quiver",
                                                  "unitig_9|quiver",
                                                  "unitig_33|quiver",
                                                  "unitig_7|quiver",
                                                  "unitig_5|quiver",
                                                  "unitig_34|quiver"))

#fsttest <- read.table("pf_fst_data.txt", header=TRUE)

##########
### Calculate pairwise Fst
##########

#pops_slide <- sliding.window.transform(pops, width = 100, jump = 10, type = 1, whole.data = FALSE)
pops_slide <- sliding.window.transform(pops, width = 10000, jump = 1000, type = 2, whole.data = FALSE)
pops_slide <- F_ST.stats(pops_slide)

## Make beds from the data ... i used to have leadingZeroScrubber( around these
outbreak12_bed <- naScrubber(leadingZeroScrubber(bedMaker(pops_slide, 1)))
outbreak13_bed <- naScrubber(leadingZeroScrubber(bedMaker(pops_slide, 2)))
outbreak23_bed <- naScrubber(leadingZeroScrubber(bedMaker(pops_slide, 3)))

##########
### SET PLOT PARAMS AND PLOT DATA
##########

#svg("circlize_fst.svg", width = 7, height = 7)
tiff("circlize_fst.tiff", width = 7, height = 7, units = "in", res = 300, compression = "lzw")

## set params
circos.par(start.degree = 137, gap.degree = 0, unit.circle.segments = 100)

## plot data
circos.genomicInitialize(chr_info, plotType = NULL) 
draw.sector(45.5, 90.8,  col = "#FF000080", border = NA, rou1 = 0.94, rou2 = 0.50, clock.wise = FALSE)
axisAdder()
trackPlotter(outbreak12_bed, fstDrawer)
trackPlotter(outbreak13_bed, fstDrawer)
trackPlotter(outbreak23_bed, fstDrawer)

## add the track names
circos.updatePlotRegion(sector.index = "Ab_names", track.index = 1, bg.border = NA)
circos.text(mean(get.cell.meta.data("xlim")), 0.5, labels = "")

circos.updatePlotRegion(sector.index = "Ab_names", track.index = 2, bg.border = NA)
circos.text(mean(get.cell.meta.data("xlim")), 0.5, labels = "Outbreak 1 - Outbreak 2", cex = 1.17, facing = "bending.inside")

circos.updatePlotRegion(sector.index = "Ab_names", track.index = 3, bg.border = NA)
circos.text(mean(get.cell.meta.data("xlim")), 0.5, labels = "Outbreak 1 - Outbreak 3", cex = 0.98, facing = "bending.inside")

circos.updatePlotRegion(sector.index = "Ab_names", track.index = 4, bg.border = NA)
circos.text(mean(get.cell.meta.data("xlim")), 0.5, labels = "Outbreak 2 - Outbreak 3", cex = 0.78, facing = "bending.inside")

circos.clear()
text(0, 0, expression(paste("Between-Outbreak ", italic(F)[ST])), cex = 1.5)

dev.off()

