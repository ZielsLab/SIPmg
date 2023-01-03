pkgname <- "SIPmg"
source(file.path(R.home("share"), "R", "examples-header.R"))
options(warn = 1)
library('SIPmg')

base::assign(".oldSearch", base::search(), pos = 'CheckExEnv')
base::assign(".old_wd", base::getwd(), pos = 'CheckExEnv')
cleanEx()
nameEx("DESeq2_l2fc")
### * DESeq2_l2fc

flush(stderr()); flush(stdout())

### Name: DESeq2_l2fc
### Title: Calculating log2 fold change for HTS-SIP data.
### Aliases: DESeq2_l2fc

### ** Examples

data(physeq_S2D2)
## Not run: 
##D df_l2fc = DESeq2_l2fc(physeq_S2D2, density_min=1.71, density_max=1.75, design=~Substrate)
##D head(df_l2fc)
## End(Not run)




cleanEx()
nameEx("HRSIP")
### * HRSIP

flush(stderr()); flush(stdout())

### Name: HRSIP
### Title: (MW-)HR-SIP analysis
### Aliases: HRSIP

### ** Examples

data(physeq_S2D2_l)

## Not run: 
##D # HR-SIP on just 1 treatment-control comparison
##D ## 1st item in list of phyloseq objects
##D physeq = physeq_S2D2_l[[1]]
##D ## HR-SIP
##D ### Note: treatment-control samples differentiated with 'design=~Substrate'
##D df_l2fc = HRSIP(physeq, design=~Substrate)
##D head(df_l2fc)
##D 
##D ## Same, but multiple BD windows (MW-HR-SIP) & run in parallel
##D ### Windows = 1.7-1.73 & 1.72-1.75
##D doParallel::registerDoParallel(2)
##D dw = data.frame(density_min=c(1.7, 1.72), density_max=c(1.73, 1.75))
##D df_l2fc = HRSIP(physeq_S2D1_l[[1]],
##D                 design=~Substrate,
##D                 density_windows=dw,
##D                 parallel=TRUE)
##D head(df_l2fc)
## End(Not run)




### * <FOOTER>
###
cleanEx()
options(digits = 7L)
base::cat("Time elapsed: ", proc.time() - base::get("ptime", pos = 'CheckExEnv'),"\n")
grDevices::dev.off()
###
### Local variables: ***
### mode: outline-minor ***
### outline-regexp: "\\(> \\)?### [*]+" ***
### End: ***
quit('no')
