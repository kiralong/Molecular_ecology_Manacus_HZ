#Original code for gghybrid written by Richard Ian Bailey 17 April 2018#
#Adapted code viewed here written by Kira M Long and Angel G Rivera-Colon

#Install the package from GitHub#
#install.packages("devtools")
#library(devtools)
#devtools::install_github("ribailey/gghybrid")

#Attach it#
library(gghybrid)
library(data.table)

# User defined variables
work.dir <- "/projects/aces/kira2/telomeres_collab/populations_runs/populations_telomeres_2910_p3_r0.8_maf0.03_hwe_structure_36k_wl"
input.csv <- "populations.p3.r80.maf03.36k.structure.csv"
par1.id <- "020SS"
par2.id <- "100CG"
pop.ids <- c("020SS", "090PR", "100CG")
legend.ids <- c("Can", "Hyb", "Vit")
# Optional
maf <- 0.1
min.diff <- 0.1
n.iterations <- 10000
n.burnin <- 1000


#Take a look at the function help files#
#?read.data	#Read in a data file in structure format or similar#
#?data.prep	#Prepare the data for analysis#
#?esth	#Run hybrid index estimation#
#?plot_h	#Plot estimated hybrid indices and credible intervals#
#?ggcline	#Run genomic cline estimation#
#?plot_clinecurve	#Plot one or more fitted cline curves, with or without individual data points#
#?compare.models	#compare two models run on the same data set using the widely applicable information criterion#

#(Note: gghybrid relies on the data.table package for data manipulation)#

#Example: Using a data file in structure format but with a complete header row.
#The file contains 1 column per marker (equivalent to ONEROW=0 in structure) for a set
#of diploid markers, plus haploid (mitochondrial) ND2. For ND2, the second (non-existent)
#allele copy is coded as missing data. Other file formats would require more information
#when running 'read.data' (see ?read.data).

#The data file contains two mandatory columns to the left of the marker columns, named
#INDLABEL (the individual references) and POPID (the sample population). These columns
#don't need to be named in the header row (or could have different names), but it makes it easier if they are.

#Set working directory
#UIUC Campus Cluster
setwd(work.dir)

#Read in the data file (This is the simplest file format for reading in)#

dat=read.data(input.csv, nprecol=2, MISSINGVAL=NA)

#Data preparation and filtering. Here I'm filtering out loci that have a minor allele
#frequency greater than 0.1 in both parental reference sets. There are also options for
#filtering by difference in parental allele frequencies, and for number of allele copies
#in each parental reference set (variable among loci due to missing data).

#The function uses objects produced by 'read.data'#

prepdata=data.prep(data=dat$data,
                   loci=dat$loci,
                   alleles=dat$alleles,
                   S0=c(par1.id), #POPID names for the first parental reference set#
                   S1=par2.id, #POPID names for the second parental reference set#
                   precols=dat$precols,
                   max.S.MAF = maf,	#Filtering by parental minor allele frequency# @KML: Original cutoff in gghybrid example script was at 10%
                   min.diff = min.diff,
                   return.genotype.table=T,
                   return.locus.table=T)

write.table(prepdata$locus.data,"locus_clines.txt",quote = FALSE, row.names = FALSE)

#'return.genotype.table=T' makes an optional table of genotypes, where for each locus
#an individual's genotype (assuming diploidy) will be 0 (two copies of the allele with
#relatively higher frequency in the 'S0' parental set), 1 (heterozygote), 2 (two copies
#of the designated 'S1' allele). This table isn't needed in downstream functions, but
#could be useful e.g. for estimating parental linkage disequilibria (associations of
#alleles from the same parent species).

#'return.locus.table=T' is also optional and not needed downstream. It's just a table
#with one row per marker, giving some information on parental allele frequencies, sample
#sizes etc.

#Next, run hybrid index estimation#

#This function uses objects produced by both the previous functions#

hindlabel= esth(data.prep.object = prepdata$data.prep,
                read.data.precols = dat$precols,
                include.Source = TRUE,	#Set to TRUE if you want hybrid indices for the parental reference individuals#
                nitt=n.iterations,
                burnin=n.burnin)

#'esth' has more functionality - this above just shows the basics#

#data.tables sometimes have a strange habit of not showing up the first
#time - if that happens just run the above line again.

###

#Plot a subset of the estimated hybrid indices (the resulting object 'abc' is useful for making a legend)#
pdf('./hybrid_index.pdf',8,8) #To make a pdf of the graph
#par(mar=c(5,5,4,2)) #To resize the margins of the graph on the pdf

setkey(hindlabel$hi,POPID)	#function from data.table, for rapid sorting and subsetting#

par(mar=c(5,5,4,5)) #To resize the margins of the graph for hybrid index
#
abc = plot_h(data=hindlabel$hi[pop.ids],#Subset of POPIDs#
             test.subject=hindlabel$test.subject,
             mean.h.by="POPID",	#Calculate the mean hybrid index for each value of the "POPID" column#
             sort.by=c("mean_h","POPID","h_posterior_mode"),	#Order test subjects along the x axis by the mean hybrid
             #index calculated above and also by individual hybrid index#
             col.group="POPID",
             group.sep="POPID",
             fill.source=TRUE,
             basic.lines=FALSE,
             source.col=c("dodgerblue","firebrick"),
             source.limits=c("dodgerblue","firebrick"),
             cex=1,pch=16,
             cex.lab=1.5,cex.main=1.5,ylim=c(0,1),
             las=1)
#

#Reshape the plot window as you want#

#Add a legend using the 'plot_h' object 'abc'#

setkey(abc,rn)		#Order data by row number#
legend("bottomright",	#Place the legend in the top left of the figure# @KML: Overwrote legend position to outside the plot
       legend.ids,   # @KML: Override the legend names to have the right IDs for CG and SS. Careful as it overrides the data itself.
       # abc[,POPID], 		#Name of the field by which data point colours are grouped# (@KML: These are the original legend labels)
       bg="white",			#Background colour#
       text.col=c("black"), #Text colour#
       pch=22, 				#Text size#
       col=abc[,col.Dark2], #Name of the field containing colour information#
       pt.bg=abc[,col.Dark2],	#Name of the field containing colour information#
       ncol=1,				#Number of columns for the group names#
       cex=1, pt.cex=1,
       xpd = TRUE,
       title = "Population")
dev.off() #to turn off the pdf, if we were separating the graphs to their own pdfs
###

write.table(hindlabel,file = "HI_table.txt",quote = FALSE, row.names = FALSE)
