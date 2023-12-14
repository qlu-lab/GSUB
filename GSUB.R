options(stringsAsFactors = FALSE)

repos <- "https://cloud.r-project.org"
# Ensure user has all necessary packages
if (!require(tidyverse)) {
  install.packages("tidyverse", repos = repos)
  library(tidyverse)
}
if (!require(data.table)) {
  install.packages("data.table", repos = repos)
  library(data.table)
}
if (!require(optparse)) {
  install.packages("optparse", repos = repos)
  library(optparse)
}
if (!require(GenomicSEM)) {
  if (!require(devtools)) {
    install.packages("devtools", repos = repos)
    library(devtools)
  }
  install_github("GenomicSEM/GenomicSEM")
  library(GenomicSEM)
}

# Read the arguments into R
option_list <- list(
  # required
  make_option("--sumstat_files", action = "store", default = NULL, type = "character"),
  make_option("--output_path", action = "store", default = NULL, type = "character"),
  make_option("--correction", action = "store", default = "standard", type = "character"),
  # optional
  make_option("--N", action = "store", default = NULL, type = "character"),
  make_option("--info.filter", action = "store", default = 0.9, type = "numeric"),
  make_option("--maf.filter", action = "store", default = 0.01, type = "numeric"),
  make_option("--sample.prev", action = "store", default = NULL, type = "character"),
  make_option("--population.prev", action = "store", default = NULL, type = "character"),
  make_option("--se.logit", action = "store", default = NULL, type = "character"),
  make_option("--OLS", action = "store", default = NULL, type = "character"),
  make_option("--linprob", action = "store", default = NULL, type = "character"),
  make_option("--keep.indel", action = "store", default = FALSE, type = "logical"),
  make_option("--hm3", action = "store", default = "./w_hm3.noMHC.snplist", type = "character"),
  make_option("--ld", action = "store", default = "./eur_w_ld_chr/", type = "character"),
  make_option("--wld", action = "store", default = "./eur_w_ld_chr/", type = "character"),
  make_option("--ref", action = "store", default = "./reference.1000G.maf.0.005.txt", type = "character")
)

### step 0. specify the inputs: ============
opt <- parse_args(OptionParser(option_list=option_list))
# Required inputs
sumstats_file <- as.character(unlist(strsplit(opt$sumstat_files, split = ",")))
num_sumstats_files <- length(sumstats_file)
# throw error if the number of sumstats files provided isn't 2
if (num_sumstats_files != 2) {
  stop(paste0("You provided ", num_sumstats_files, " sumstats files. GSUB only works for 2 sumstats files (GWAS a minus b)."))
}
outdir <- opt$output_path
correction <- opt$correction

# Helper function for parsing arguments which are vectors of characters
# sed for N, sample.prev, population.prev, se.logit, OLS, and linprob
# vector size is based on how many sumstats files there are
get_vector_option <- function(option, type, default) {
  if (!is.null(option)) {
    if (type == "numeric") {
      v <- suppressWarnings(as.numeric(unlist(lapply(strsplit(option, ","), trimws))))
    } else if (type == "logical") {
      v <- as.logical(unlist(lapply(strsplit(option, ","), trimws)))
    } else {
      v <- as.character(unlist(lapply(strsplit(option, ","), trimws)))
    }
    if (length(v) != num_sumstats_files) {
      stop(paste0("The length of ", option, " you provided does not equal the required number 2 - the number of sumstat files."))
    }
  } else {
    v <- rep(default, num_sumstats_files)
  }
  return(v)
}

# Optional inputs
# sample sizes, if sumstats file do not contain N column:
# if want to use N in sumstats file, use NA here like this: N = c(NA, 100000)
N <- get_vector_option(opt$N, "numeric", NA) # default is NA
# Set up munged file directory with correct amount of file names
munged_file <- c()
for (i in 1:num_sumstats_files) {
  munged_file <- append(munged_file, paste(opt$output_path, "/file_", i, ".munge", sep = ""))
}
# quality controls
info.filter <- opt$info.filter
maf.filter <- opt$maf.filter
# sample and population prevalence to compute LDCS outputs in liability scale
sample.prev <- get_vector_option(opt$sample.prev, "numeric", NA) # default is NA
population.prev <- get_vector_option(opt$population.prev, "numeric", NA) # default is NA
# settings for cleaning and reformatting sumstats
se.logit <- get_vector_option(opt$se.logit, "logical", TRUE) # default is TRUE
OLS <- get_vector_option(opt$OLS, "logical", FALSE) # default is FALSE
linprob <- get_vector_option(opt$linprob, "logical", FALSE) # default is FALSE
# references and resources - in most cases users should use the defaults
hm3 <- opt$hm3
ld <- opt$ld
wld <- opt$wld
ref <- opt$ref
# parallel processing to speed up GSUB
keep.indel <- opt$keep.indel # default is FALSE
# parallel is true by default with 2 cores because this speeds up GSUB
# and we can use a maximum of 2 cores because we have a maximum of 2
# sumstats files
parallel <- TRUE
cores <- 2

### step 1. munge sumstats: ===================
#
# munge() computes z-score from p values and its sign is the sign of effect size

cat("\n\nMunge input GWAS sumstats ... \n")
munge(files = sumstats_file,
      hm3 = hm3, trait.names = munged_file,
      N = N,
      info.filter = info.filter,
      maf.filter = maf.filter,
      parallel = parallel, cores = cores,
      log.name = paste0(outdir, "/gsem"))

# log file will be outdir/gsem_munge.log
# will write a "outdir/trait.names.sumstats.gz" file for each input sumstats

# be careful not overwrite the original GWAS sumstats!
system(paste0("cat ", outdir, "/gsem_munge.log"))


### step 2. bivariate LDSC: ========================

cat("\n\nBivariate LDSC ... \n")

LDSCoutput <- ldsc(traits = paste0(munged_file, ".sumstats.gz"),
                  sample.prev = sample.prev,
                  population.prev = population.prev,
                  ld = ld, wld = wld,
                  stand = TRUE, ldsc.log = paste0(outdir, "/gsem"))

# check results:
print(LDSCoutput)

# save LDSC results for later use:
save(LDSCoutput, file = paste0(outdir, "/LDSCoutput.RData"))

# log file will be at outdir/gsem_ldsc.log
system(paste0("cat ", outdir, "/gsem_ldsc.log"))


### step 3. Run the SEM model without SNPs ======================
#
# the results from this step won't be used in later steps
# not sure exactly what's the purpose of this step
# maybe this is some kind of sanity check: this step should return
# reasonable results, otherwise model with SNPs will not make sense?
# 
# note, this step is not needed for applying analytical solutions
# it is included here only for the sake of completeness

cat("\n\nRunning model without SNP ... \n")

# GSEM model withou using SNP:
model_noSNP <- 'F1=~NA*V2 + start(0.4)*V1
           F2=~NA*V2
         
           F2~~1*F2
           F1~~1*F1
           F1~~0*F2

           V1~~0*V2
           V1~~0*V1
           V2~~0*V2'

output <- usermodel(LDSCoutput, estimation = "DWLS", model = model_noSNP)

cat("\n\nModel without SNP results:\n")
print(output)

# save the model without SNP results:
save(output, file = paste0(outdir, "/Model.NoSNP.output.Rdata"))


### step 4. Clean GWAS sumstats: =====================================
#
# this returns a data.frame with aligned sumstats
# this step may give a warning about its log file due to bug in its code
# which is fine
#
# note that the sumstats() function requires that the variables be listed in the same order 
# as they are listed for the ldsc() function above.
#
# coefficients and their SEs are then further transformed such that
# they are scaled relative to unit-variance scaled phenotypes
#
# when in-sample MAF is not available, the function uses reference MAF to back out 
# a partially standardized beta, and mismatch between in-sample and reference MAF 
# may cause some minimal amount of bias in the relative scaling of the betas across summary statistics
#
# Note that it is possible to back out a beta and SE when only Z-statistics are provided 
# as long as you have the total sample size for a continuous trait
# or the sum of effective sample sizes across the contributing cohorts for a binary trait.

cat("\n\nCleaning input GWAS sumstats ... \n")
p_sumstats <- sumstats(files = sumstats_file,
                      ref = ref,
                      se.logit = se.logit,
                      OLS = OLS,
                      linprob = linprob,
                      N = N,
                      info.filter = info.filter,
                      maf.filter = maf.filter,
                      keep.indel = keep.indel,
                      parallel = parallel,
                      cores = cores)

# this results a data.fame with aligned sumstats
#
# files: The name of the summary statistics files. This should be the same
#        as the name of the files used for the munge function in Step 1 above
#        and the files should be in the same listed order used for the ldsc function in step 2.
#
# ref: The reference file used to calculate SNP variance across traits.
#
# log file: "_sumstats.log" in the current directory since this function does
#           not allow use to specify the log filename and directory
#     therefore, let's re-name it and move it to the "outdir" directory
#     and check its contents:
cat("\n\nThere are ", format(nrow(p_sumstats), big.mark = ","), " SNPs left after cleaning.\n", sep = "")

system(paste0("mv _sumstats.log ", outdir, "/gsem_clean_sumstats.log"))
system(paste0("cat ", outdir, "/gsem_clean_sumstats.log"))

# export the entire clean'ed sumstats as one file for future reference:
cat("\n\nExporting the cleaned (QC and aligned) sumstats ... \n")
fwrite(p_sumstats, paste0(outdir, "/gsem_sumstats.cleaned.txt.gz"),
       col.names = TRUE, row.names = FALSE, sep = "\t", quote = FALSE, na = NA)

cat("\nThe cleaned (QC'ed and aligned) sumstats were written to ", outdir, "/gsem_sumstats.cleaned.txt.gz\n", sep = "")


### step 5. Run the GWAS-by-subtraction: ====================
#
# note, in this step, we'll use analytical solutions
# to estimate the loading of the two edges
#
# apply analytical solutions:
begin.time <- Sys.time()
cat("\n\nApplying GWAS-by-subtraction (GSUB) using analytical solutions ... begins at ", format(begin.time, format = "%F %R %Z"), ":", sep = "", "\n\n")

h2_1 <- LDSCoutput[["S"]][1]
h2_2 <- LDSCoutput[["S"]][4]

gcov <- LDSCoutput[["S"]][2]
gcov_int <- LDSCoutput[["I"]][2]

h2_1_int <- LDSCoutput[["I"]][1, 1]
h2_2_int <- LDSCoutput[["I"]][2, 2]

lambda11 <- sqrt(h2_1)
lambda12 <- gcov / lambda11

hh <- h2_2 - lambda12^2
if (hh < 0) {
  cat("\nlambda22 is set to NA due to negative value under square root.\n")
  lambda22 <- NA
} else {
  lambda22 <- sqrt(hh)
  cat("\nlambda22 = ", lambda22, "\n", sep = "")
}

# apply correction to standard error
apply_correction <- function(se, h2_int) {
  if (correction == "standard") {
    return(se * sqrt(h2_int))
  } else if (correction == "conserv") {
    return(se * h2_int)
  }
  return(se)
}

p_sumstats <- p_sumstats %>%
  mutate(lambda11 = lambda11,
         lambda12 = lambda12,
         lambda22 = lambda22) %>%
  mutate(se.1 = apply_correction(se.1, h2_1_int),
         se.2 = apply_correction(se.2, h2_2_int)) %>%
  mutate(gamma1 = beta.1 / lambda11,
         se.gamma1 = se.1 / lambda11) %>%
  mutate(z.gamma1 = gamma1 / se.gamma1,
         p.gamma1 = pnorm(abs(gamma1 / se.gamma1), lower.tail = FALSE) * 2,
         gamma2 = (beta.2 - gamma1 * lambda12)/lambda22,
         se.gamma2 = sqrt(se.2^2 + (se.1 * lambda12 / lambda11)^2 - 2 * gcov_int * se.1 * se.2 * lambda12 / lambda11) / lambda22) %>%
  mutate(z.gamma2 = gamma2 / se.gamma2,
         p.gamma2 = pnorm(abs(gamma2 / se.gamma2), lower.tail = FALSE) * 2)

# add effective sample size:
# https://github.com/GenomicSEM/GenomicSEM/wiki/5.-User-Specified-Models-with-SNP-Effects
cat("\nCompute effective sample size using 1/(2p(1-p)SE^2) for common SNPs ...\n")

tt <- p_sumstats %>% filter(MAF <= 0.4 & MAF >= 0.1)
neff1 <- as.integer(mean(1 / ((2 * tt$MAF * (1 - tt$MAF)) * tt$se.gamma1^2), na.rm = TRUE))
neff2 <- as.integer(mean(1/((2 * tt$MAF * (1 - tt$MAF)) * tt$se.gamma2^2), na.rm = TRUE))

p_sumstats$n.gamma1 <- ifelse(p_sumstats$MAF > 0.4 | p_sumstats$MAF < 0.1, NA, neff1)
p_sumstats$n.gamma2 <- ifelse(p_sumstats$MAF > 0.4 | p_sumstats$MAF < 0.1, NA, neff2)

cat("Effective sample size: ", format(neff1, big.mark = ","), " for loading 1 and ", format(neff2, big.mark = ","), " for loading 2.\n", sep = "")

cat("\nSummary of the results:\n")
print(summary(p_sumstats))

# export:
cat("\n\nWrite results to ", outdir, "\nWill also return the results as a data.frame\n", sep = "")

fwrite(p_sumstats %>% select(CHR, SNP, BP, A1, A2, 
                             BETA = gamma1, SE = se.gamma1, P = p.gamma1, N = n.gamma1),
       file = paste0(outdir, "/gsub_analytical_loading1.txt.gz"),
       col.names = TRUE, row.names = FALSE, quote = FALSE, sep = "\t", na = NA)

fwrite(p_sumstats %>% select(CHR, SNP, BP, A1, A2, 
                             BETA = gamma2, SE = se.gamma2, P = p.gamma2, N = n.gamma2),
       file = paste0(outdir, "/gsub_analytical_loading2.txt.gz"),
       col.names = TRUE, row.names = FALSE, quote = FALSE, sep = "\t", na = NA)

# report computing times:
end.time <- Sys.time()
total.time <- difftime(time1 = end.time, time2 = begin.time, units = "sec")
hrs <- floor(floor(total.time) / 3600)
mins <- floor(floor(floor(total.time) - hrs * 3600) / 60)
secs <- total.time - hrs * 3600 - mins * 60
cat("\nGSUB using analytical solutions ended at ", format(end.time, format = "%F %R %Z"), ". \nTotal time: ", hrs, " hours, ", mins, " minutes and ", secs, " seconds.\n", sep = "")

# export final gsub analytical results:
cat("\nWriting the GSUB using analytical solutions results to ", outdir, "/gsub_analytical_all.txt.gz\n", sep = "")
fwrite(p_sumstats,
       file = paste0(outdir, "/gsub_analytical_all.txt.gz"),
       col.names = TRUE, row.names = FALSE, quote = FALSE, sep = "\t", na = NA)