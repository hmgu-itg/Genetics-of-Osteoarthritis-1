###########################################
###      FileCheck-PreQC .README        ###
###########################################
### Creator: Cindy G. Boer
### Last Updated: 27/07/2018
### For questions: c.boer@erasmumc.nl/cb39@sanger.ac.uk
###########################################

###### SCRIPT INFORMATION
### minimum: R version (3.0.0)
### dependancies: library("data.table")
### Needs ~20-30 g memmory for good running
### Runtime 1 file, 30g mem ~ 20-30 min

### Run as:
# Run Local
# ./FileCheck-PreQC.R [PATH/filename.Uploaded.GWAS.data]

## RUN for directory with files[FARM]
# for i in $(ls PATH/FOLDER/*); do echo /FileCheck-PreQC.R  $i; done | ~ag15/array 30g NAME

### Input:
### Takes uploaded GWAS.summary Stats [as defined in analysis.plan, includes direct SNPtest output]. File can be in .txt or in .gz format.

### SCRIPT:
## FileCheck-PreQC.R
## R-script takes argument from command line [filename], needs library (data.table), runs pre-QC checks:
# - filename correct
# - NCASES/Ncotrols present in file name
# - Harmonizing headernaming [SNPtest -output to desired file format]
# - Calculation BETA from OR, if BETA is not present
# - Adding NCASES/NCONTROLS and total N to file
# - Reordering columns [for easy harmonizing across files]
# - Saving PreQC_filename.txt which is compressed into .txt.gz
# - output file suitable for Go.EasyQC scripts [Step_02]
#
# Output contains following data:
# - CTPID 	: unique identifier variant: CHR:POS
# - CHR		: chromosome position, HG reference can be found in output file name
# - POS		: base pair position of variant, HG reference can be found in output file name 
# - EA		: effect allele [A/T/C/G]
# - NEA		: non effect allele [A/T/C/G]
# - EAF		: effect allele frequency 
# - BETA	: beta of GWAS results for variant
# - SE		: Standard Error of the BETA
# - P		: Pvalue of BETA
# - NCASES	: nr of CASE individuals in analysis
# - NCONTROLS	: nr of CONTROL individuals in analysis
# - INFO	: metric of imputation quality [0-1]
# - N		: Total amount of individuals in test [NCASES + NCONTROLS]
# NO Filtering Steps or Cleaning steps are performed in this script, NO variant are excluded or tossed.
# File output has "PreQC_" prefix

#!/software/bin/Rscript --vanilla

#########################################################
#### Genetics of Osteoarthritis FileChecK-PreQC     #####
#########################################################
### GO consortium
### By Cindy G. Boer Erasmus MC
### c.boer@erasmusmc.nl / cb39@sanger.ac.uk
### Version 4.1
### Updates; With fudgefactor, EAF SNPtest, print warnings
### Last-updated: 18/10/2018
### Run: ./FileCheck-PreQC.r [arg1:input.file] 
#########################################################
###                 Dependancies                      ###
#########################################################
### library "data.table
# on farm fwrite does not work, fread does
library("data.table")
setDTthreads(threads = 1)

#########################################################
###                   Arguments                       ###
#########################################################
args <- commandArgs(trailingOnly=TRUE)

### Input data
input.file <- args[1] #GWAS summary stats file
#name output data is based on input data name

#########################################################
###                      Functions                    ###
#########################################################

### File Name Checker: Case/Control information at the right position
FileNameCheck <- function(a){
  if ((grepl("CASE", a) & grepl("CONTROL", a)) == TRUE) {
    print("N of Cases and Controls present in file name, moving to next step..")
  } else { if (((is.numeric(as.numeric(sapply(strsplit(basename(a), "_"), "[[", 4)))) &
                  (is.numeric(as.numeric(sapply(strsplit(basename(a), "_"), "[[", 5))))) == TRUE ) {
                  print("N of Cases and Controls present in file name, moving to next step..")
    } else {
      print("No N of cases and controls present in file name. Stopping FileChecK-PreQC")
      q()
    }
  }
}

### File Reader: loads in files .txt, or .gz
FileReader <- function(a){
  if ((grepl("*.gz",a))== TRUE) {
    b <- as.character(paste0("zcat ",a))
    y <-fread(b, header= TRUE, stringsAsFactors = FALSE, data.table = TRUE)
       return(y)
  } else {
	y <-fread(a, header = TRUE, stringsAsFactors = FALSE, data.table = TRUE)
    	return(y)
  }
}

### header nameing harmonization tool
HeaderHarmonization <- function(a){
  fl <- paste((names(a)), collapse = ' ')
  if (((grepl("SNPTEST", (basename(as.character(input.file)))) | grepl("snptest", (basename(as.character(input.file))))) &(!grepl("EAF|eaf",fl))) == TRUE){
  colnames(a) <- gsub("pos.*|POS.*|POSB37", "POS", colnames(a))
  colnames(a) <- gsub("chr.*|CHR.*", "CHR", colnames(a))
  colnames(a) <- gsub("all_OR|ALL_OR|all_or", "OR", colnames(a))
  colnames(a) <- gsub("frequentist_add_pvalue", "P", colnames(a))
  colnames(a) <- gsub("frequentist_add_beta_1", "BETA", colnames(a))
  colnames(a) <- gsub("frequentist_add_se_1", "SE", colnames(a))
  colnames(a) <- gsub("info", "INFO", colnames(a))
  colnames(a) <- gsub("alleleB", "EA", colnames(a))
  colnames(a) <- gsub("alleleA", "NEA", colnames(a))
  colnames(a) <- gsub("pval.*|Pval.*|PVAL.*", "P", colnames(a))
  colnames(a) <- gsub("se.*", "SE", colnames(a))
  colnames(a) <- gsub("beta.*", "BETA", colnames(a))
  colnames(a) <- gsub("info.*", "INFO", colnames(a))
  colnames(a) <- gsub("effect_allele|Effect_Allele", "EA", colnames(a))
  colnames(a) <- gsub("other_allele|Other_Allele", "NEA", colnames(a))
  colnames(a) <- gsub("eaf", "EAF", colnames(a))
  colnames(a) <- gsub("ea", "EA", colnames(a))
  colnames(a) <- gsub("nea", "NEA", colnames(a))
  a$CHR <- gsub(":.*", "", a$CHR)
  return(a)
  } else {
  fl <- paste((names(a)), collapse = ' ')  
  if ((grepl("\\bCHR\\b", fl) & grepl("\\bPOS\\b", fl) & grepl("\\bEA\\b", fl) &
       grepl("\\bNEA\\b", fl) & grepl("\\bP\\b", fl) & grepl("INFO", fl) &
       ((grepl("\\bBETA\\b", fl) & grepl("\\bSE\\b", fl))|(grepl("\\bOR\\b", fl) & grepl("\\bOR_L95\\b", fl) & grepl("\\bOR_U95\\b", fl)))) == TRUE){
      return(a)
      } else{
        colnames(a) <- gsub("frequentist_add_pvalue", "P", colnames(a))
        colnames(a) <- gsub("frequentist_add_beta_1", "BETA", colnames(a))
  	colnames(a) <- gsub("frequentist_add_se_1", "SE", colnames(a))
        colnames(a) <- gsub("pos.*|POS.*|POS.*|POSB37", "POS", colnames(a))
        colnames(a) <- gsub("chr.*|CHR.*", "CHR", colnames(a))
        colnames(a) <- gsub("all_OR|ALL_OR|all_OR", "OR", colnames(a))
        colnames(a) <- gsub("all_OR_lower|ALL_OR_lower|ALL_OR_Lower|L95", "OR_L95", colnames(a))
        colnames(a) <- gsub("all_OR_upper|ALL_OR_upper|ALL_OR_Upper|U95", "OR_U95", colnames(a))
        colnames(a) <- gsub("pval.*|Pval.*|PVAL.*", "P", colnames(a))
	colnames(a) <- gsub("se.*", "SE", colnames(a))
	colnames(a) <- gsub("beta.*", "BETA", colnames(a))
	colnames(a) <- gsub("info.*", "INFO", colnames(a))
	colnames(a) <- gsub("effect_allele|Effect_Allele", "EA", colnames(a))
  	colnames(a) <- gsub("other_allele|Other_Allele", "NEA", colnames(a))
	colnames(a) <- gsub("eaf", "EAF", colnames(a))
	colnames(a) <- gsub("ea", "EA", colnames(a))
	colnames(a) <- gsub("nea", "NEA", colnames(a))
        return(a)
      }
  }
}

### header checker: are all columns present?
HeaderCheck <- function(a){
  fl <- paste((names(a)), collapse = ' ')
  if ((grepl("\\bCHR\\b", fl) & grepl("\\bPOS\\b", fl) & grepl("\\bEA\\b", fl) &
       grepl("\\bNEA\\b", fl) & grepl("\\bP\\b", fl) & grepl("\\bINFO\\b", fl) &
       ((grepl("\\bBETA\\b", fl) & grepl("\\bSE\\b", fl))|(grepl("\\bOR\\b", fl) & grepl("\\bOR_L95\\b", fl) & grepl("\\bOR_U95\\b", fl)))
  ) == TRUE) {
    print("All headers present with correct naming, moving to next step...")
  } else{
    print("One or more header(s) missing or incorrectly named. Stopping FileChecK-PreQ .")
    q()
  }
}

### CaseNrExtraction : get n of Cases from header
CaseNExtract <- function(a){
  if (grepl("CONTROL", a) == TRUE) {
    x <- sapply(strsplit((basename(as.character(a))),"_"), "[[",4)  
    matches <- regmatches(x, gregexpr("[[:digit:]]+", x))
    as.numeric(unlist(matches))	 
  } else {
    as.numeric(sapply(strsplit((basename(as.character(a))), "_"), "[[", 4))
  }
} 

### ControlNrExtraction : get n of Cases from header
ControlNExtract <- function(a){
  if (grepl("CONTROL", a) == TRUE) {
    x <- sapply(strsplit((basename(as.character(a))),"_"), "[[",5)  
    matches <- regmatches(x, gregexpr("[[:digit:]]+", x))
    as.numeric(unlist(matches))	 
  } else {
    as.numeric(sapply(strsplit((basename(as.character(a))), "_"), "[[", 5))
  }
} 
  
### BETA Checker and OR to Beta calculator

fudge <- function(b){
        b[,':='(fudgefactor = 1/((NCASES/N)*(1-(NCASES/N))))]
	b[,':='(BETA = BETA*fudgefactor)] 
        b[,':='(SE = SE*fudgefactor)]
        return(b)
}

OrToBeta <- function(a){
	if ((exists("BETA",a) == TRUE) & (exists("SE",a) == TRUE)){
  	return(a)
	  } else {
    		a[,':='(BETA =log(OR))]
    		a[,':='(SE = ((log(OR_U95)-(log(OR)))/1.95996))]
    		return(a)
	}
}

### EAF check SNPtest ######################
EAFcheck <- function(a){
  if (exists("EAF",a) == TRUE){
    return(a)
  } else {
    a[,':='(EAF = (2*(all_BB)+all_AB)/(2*(all_AB + all_BB + all_AA)))]
    return(a)
  }
}

#########################################################
###                   FileChecK-PreQ                  ###
#########################################################
### Step 1. : File naming check
# Check is file name contains case/control numbers
FileNameCheck(input.file)

### Step 2. : Reading file
print(paste0("Reading file: ",basename(input.file)))
df <-FileReader(input.file)
print("File read in as:")
print(head(df))

### Step 3. : Header check
data.file <- HeaderHarmonization(df)
print("Header Harmonized in file, moving to next step...")
rm(df)
HeaderCheck(data.file)

### Step 4. : add CASES/CONTROLS from header (file.name)

data.file[,':='(NCASES = CaseNExtract(input.file))]
data.file[,':='(NCONTROLS = ControlNExtract(input.file))]
data.file[,':='(N = (NCASES + NCONTROLS))]
data.file[,':='(CPTID = paste0(CHR,":",POS))]

print("Added Cases/Control and CPTID to data, moving to next step...")

### Step 5. : EAF calculation/check for 
y <-EAFcheck(data.file)

### Step 6. : Beta/SE calculation 
if((grepl("BOLTLMM", input.file) == TRUE) | (grepl('BOLT',input.file) ==TRUE) | (grepl("GEMMA",input.file) == TRUE)){
	x <- fudge(y)
	rm(data.file)
} else{
	x <- OrToBeta(y)
	rm(y)
}

print("BETA and SE in dataset, moving to next step...")


### Step 7. : subsetting and ordering data
x2 <- as.data.frame(x)
rm(x)
data.out  <- x2[, c("CPTID", "CHR", "POS", "EA", "NEA", "EAF", "P", "BETA", "SE","INFO", "NCASES", "NCONTROLS", "N")]
print("DataFrame made and subsetted, writing output...")
rm(x2)

### Step 8. : Writing output
name 	 <- basename(as.character(input.file)) 
out.name <-paste0("PreQC_",sapply(strsplit(name, "_"),"[[", 1),
                   "_",sapply(strsplit(name, "_"), "[[", 2), 
                   "_", sapply(strsplit(name, "_"), "[[", 3),
                   "_",sapply(strsplit(name, "_"), "[[", 6),
                   "_",sapply(strsplit(name, "_"), "[[", 7), 
                   "_", sapply(strsplit(name, "_"), "[[", 8),
                   ".txt")

# We create the .txt file
print(paste0("Output is being writen to:",out.name))
print("Output header:")
print(head(data.out))

write.table(data.out, out.name, row.names = FALSE, quote = FALSE, sep = "\t")
#fwrite(data.out, file = out.name, row.names = FALSE, append = FALSE, quote = FALSE, col.names = TRUE, sep = "\t")
rm(data.out)
print(".txt file witten, compressing to .gzip...")
# Save as .txt.gzip
out <- paste0("gzip ",out.name)
system(out)

print("PreQC-check done!")
warnings()

#########################################################
###           Exiting FileCheck-PreQC                 ###
#########################################################
#q()
#########################################################
