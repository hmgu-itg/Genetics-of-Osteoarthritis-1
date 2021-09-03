#!/software/bin/Rscript --vanilla

#########################################################
#### Genetics of Osteoarthritis Post.Metal.QC      #####
#########################################################
### GO consortium
### By Cindy G. Boer Erasmus MC
### c.boer@erasmusmc.nl
### Version 2.1
### Last-updated: 06/02/2019
#########################################################
###                 Dependancies                      ###
#########################################################
### library "data.table
# on farm fwrite does not work, fread does
library("data.table")
setDTthreads(threads = 1)

#library("qqman") starts only when plots are made

#########################################################
###                   Arguments                       ###
#########################################################
args <- commandArgs(trailingOnly=TRUE)

input.file  <- args[1] # Meta.tbl file
out         <- args[2] # ID/name file out/phenotype
path.out    <- args[3] # Path-output

#########################################################
###             Defining output.files                 ###
#########################################################

tophits     <- paste0(path.out,"GO.final.meta.results.",out, ".TopSignals.file", ".txt")
tophitsQC   <- paste0(path.out,"GO.final.meta.results.",out, ".TopSignalsQC.file", ".txt")

qqplot      <- paste0(path.out,"GO.final.meta.results.",out, ".qq.plot.freq", ".png")
manplot     <- paste0(path.out,"GO.final.meta.results.",out, ".man.plot.freq", ".png")

metaR       <- paste0(path.out,"GO.RAW.final.meta.results.",out,".txt")
metaF       <- paste0(path.out,"GO.FUMA.final.meta.results.",out,".txt")
metaFF	    <- paste0(path.out,"GO.FILTER.final.meta.results.",out,".txt")	

#########################################################
###                      Functions                    ###
#########################################################
### Meta Filter
MetaFilter <- function(q){
  q[,':='(Plus  = sapply(regmatches(Direction, gregexpr("\\+", Direction)), length))]
  q[,':='(Minus = sapply(regmatches(Direction, gregexpr("\\-", Direction)), length))]
  q[,':='(None  = sapply(regmatches(Direction, gregexpr("\\?", Direction)), length))]
  q[,':='(Total = Plus + Minus + None)]
  y <- q[which((q$Plus>= (q$Total/2))|(q$Minus >= (q$Total/2))),]
  return(y)
}

### MAF Filter
MAFFilter <- function(a){
  a[,':='(MAF = ifelse(EAF <0.5,EAF,1-EAF))]
  q <- a[which(a$MAF >= 0.005),]
  return(q)
}

#### Beta to OR calculator
OrCalculator <- function(a){
                a[,':='(OR = exp(BETA))]
                a[,':='(OR_U95 = ((exp(BETA+(1.959964*SE)))))]
                a[,':='(OR_L95 = ((exp(BETA-(1.959964*SE)))))]
                return(a)
}

#### TopHit Extract function:
BestHits <- function(a){
  z<-a[which(a$P <= 1e-6),]
  return(z)
}

#### CHR POS Function
pos <- function(x){
  y<- unlist(strsplit(x, split=':',fixed=TRUE))[2]
}

chr <- function(x){
  y<- unlist(strsplit(x, split=':',fixed=TRUE))[1]
}

#########################################################
###                   Post-Metal                      ###
#########################################################
# Reading Input files:
print((paste0("Starting Post-MetalQC, Reading:",input.file)))
df      <- fread(input.file, header = TRUE, data.table = TRUE, stringsAsFactors = FALSE)
x       <- fread("zcat Mapping_File_CPTID_SNP_AF_HRCr1.1_GRCH37.nochrX.txt.gz",
               header = TRUE, data.table = TRUE, stringsAsFactors = TRUE)

# Adding rsID's to Meta-analysis file:
print("Starting Post-MetalQC - adding rsIDs")
x2 <-as.data.frame(x)
x3       <- x2[, c("CPTID", "SNP", "CHR","POS")]
setnames(df, c("CPTID", "EA", "NEA", "EAF", "FreqSE", "MinFreq", "MaxFreq", "BETA", "SE", "P", "Direction", "HetISq", "HetChiSq", "HetDf", "HetPVal", "NCASES","NCONTROLS","N"))
m.data  <- merge(df,x3, by="CPTID", all.x = TRUE, all.y = FALSE, allow.cartesian=TRUE)

# clean-up to preserve memory:
rm(x)
rm(x2)
rm(x3)
rm(df)

# Make sure all variants have CHR and POS data
m.data$CHR<-apply(m.data, 1, chr)
m.data$POS<-apply(m.data, 1, pos)

# Adding OR to meta-analysis data
df <- OrCalculator(m.data)
head(df)

print((paste0("Saving meta-analysis results RAW file:",metaR)))
write.table(df, metaR, row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")
out <- paste0("gzip ",metaR)
system(out)

# Extracting Top signals
print("Extracting TOP Signals")
top   <- BestHits(df)
write.table(top,tophits, row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")


print("Performing Filers in TopHitsFileQC: MAF>=0.005 & direction effect same in half of cohorts...")
df_1 <- MAFFilter(df)
df_2 <- MetaFilter(df_1)
rm(df_1)
rm(df)

print((paste0("Saving meta-analysis results FIlter file:",metaFF)))
write.table(df_2, metaFF, row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")
out <- paste0("gzip ",metaFF)
system(out)

# Extracting Top signals
print("Extracting TOPQC Signals")
top   <- BestHits(df_2)
write.table(top,tophitsQC, row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")

print((paste0("Saving meta-analysis results FUMA file:",metaR)))
x <- as.data.frame(df_2)
fuma     <-x[,c("SNP","CHR","POS","EA", "NEA","EAF", "BETA", "SE", "P", "N")]
fuma_f1  <-fuma[which(fuma$P < 5e-2),]
write.table(fuma_f1, metaF, row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")
outF <- paste0("gzip ",metaF)
system(outF)

# Clean
rm(df_2)
rm(x)

########################################################
###       Post-Metal Plots: QQplot [MAF]             ###
#########################################################
### Plots QQ plot, stratified by MAF, and Lambda by MAF
print((paste0("Plotting MAF QQ-plot METAQC data:", qqplot)))
## stratify by MAF
x <- na.omit(fuma)
pvals_lo1=subset(x, ( x[,"EAF"] > 0.05 & x[,"EAF"] < 0.95 ))
pvals_lo4=subset(x, ( x[,"EAF"] < 0.05 & x[,"EAF"] > 0.01 ) | ( x[,"EAF"] > 0.95 & x[,"EAF"] < 0.99 ))
pvals_lo5=subset(x, ( x[,"EAF"] < 0.01 & x[,"EAF"] > 0.00 ) | ( x[,"EAF"] > 0.99 & x[,"EAF"] < 1.00 ))

pv<-sort(x[,"P"])
y.range=range(-log10(pv))
y.range[2]=y.range[2]+1
z=qnorm(x[,"P"]/2)
z_lo1=qnorm(pvals_lo1[,"P"]/2)
z_lo4=qnorm(pvals_lo4[,"P"]/2)
z_lo5=qnorm(pvals_lo5[,"P"]/2)

## calculates lambda
lambda = round(median(z^2)/qchisq(0.5,df=1),3)
lmb = as.expression(substitute(paste('Observed, ', lambda, '= ', lam, ' N=', (nrow(z))),list(lam = lambda)))

lambda1 = round(median(z_lo1^2)/qchisq(0.5,df=1),3)
lmb1 = as.expression(substitute(paste('MAF > 0.05, ', lambda, '= ', lam, ' N=',(nrow(pvals_lo1))),list(lam = lambda1)))

lambda4 = round(median(z_lo4^2)/qchisq(0.5,df=1),3)
lmb4 = as.expression(substitute(paste('0.01 < MAF < 0.05, ', lambda, '= ', lam, ' N=',(nrow(pvals_lo4))),list(lam = lambda4)))

lambda5 = round(median(z_lo5^2)/qchisq(0.5,df=1),3)
lmb5 = as.expression(substitute(paste('0 < MAF < 0.01, ', lambda, '= ', lam, ' N=',(nrow(pvals_lo5))),list(lam = lambda5)))

## Plots axes and null distribution
qqname <- paste0("GO.QQ-Plot.",out)
bitmap(qqplot, width = 10, height = 10)
plot(y.range, y.range, col = "red", lwd = 3, type = "l", xlab = "Expected Distribution (-log10 of P value)",
     ylab = "Observed Distribution (-log10 of P value)", xlim = c(0,8), ylim = y.range,
     las = 1, xaxs = "i", yaxs="i", bty = "l", main = qqname)

##function
plotQQ <- function(z,color,cex){
  p <- 2*pnorm(-abs(z))
  p <- sort(p)
  expected <- c(1:length(p))
  lobs <- -(log10(p))
  y.range<-range(lobs)
  lexp <- -(log10(expected / (length(expected)+1)))
  
  # plots all points with p < 1e-3
  p_sig = subset(p,p<0.001)
  points(lexp[1:length(p_sig)], lobs[1:length(p_sig)], pch=23, cex=.3, col=color, bg=color)
  
  # samples 2500 points from p > 1e-3
  n=2501
  i<- c(length(p)- c(0,round(log(2:(n-1))/log(n)*length(p))),1)
  lobs_bottom=subset(lobs[i],lobs[i] <= 3)
  lexp_bottom=lexp[i[1:length(lobs_bottom)]]
  
  if (length(lobs_bottom) != length(lexp_bottom)) {
    lobs_bottom <- lobs_bottom[1:length(lexp_bottom)]
  }
  
  points(lexp_bottom, lobs_bottom, pch=23, cex=cex, col=color, bg=color)
}


## plots data
plotQQ(z,"black",0.4);
plotQQ(z_lo1,"#018571",0.3);
plotQQ(z_lo4,"#f4a582",0.3);
plotQQ(z_lo5,"#ca0020",0.3);
legy<-y.range[2]-1

## provides legend
legend(.25,legy,bty="n", legend=c("Expected (null)",lmb,lmb1,lmb4,lmb5),
  pch=c((vector("numeric",5)+1)*23), cex=c((vector("numeric",5)+0.7)),
  pt.bg=c("red", "black", "#018571","#f4a582","#ca0020"))
dev.off()

#########################################################
###       Post-Metal Plots: Manhattan plot            ###
#########################################################
## Manhattan plot
x$CHR <- as.numeric(x$CHR)
x$POS <- as.numeric(x$POS)
library("qqman")
print((paste0("Plotting manhattan plot:",manplot)))
plot1 <- manplot
manname <- paste0("GO.Man-Plot.",out)
bitmap(plot1, width = 16, height = 10)
manhattan (x, chr = "CHR", bp = "POS", p = "P", snp = "SNP", col =
             c("gray82","gray51", "royalblue1"), ymax = NULL, suggestiveline = -log10(1e-06), genomewideline = -log10(5e-08),
           highlight = NULL, main = manname )
dev.off()

############################################################
###                  SCRIPT DONE                         ###
############################################################
print("Post-meta-analysis script done")
q()
