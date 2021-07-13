# PURPOSE ----
# Reformat needed values for use with CausalPath from the results of a mass 
# spectrometry experiment. This needs customization for each project. 

library(readxl)
library(stringr)

# READ DATA ----
# Supplementary File 2 from https://www.sciencedirect.com/science/article/pii/S0006497120799349
dat <- read_excel("./gpvi_mass_spec/1-s2.0-S0006497120799349-mmc2.xlsx", 
                  sheet="Cond #2 TiO2 Brief", skip=11)

# INITIALIZE RESULTS DATA.FRAME ----
reformatted_dat <- data.frame(ID=character(0), Symbols=character(0),
                              Sites=character(0), Effect=character(0),
                              SignedP=numeric(0), stringsAsFactors=FALSE)

# EXTRACT VALUES ----
pb <- txtProgressBar(min=1, max=nrow(dat), style=3)
for(i in 1:nrow(dat)) {
  setTxtProgressBar(pb, i)

  gene_symbol <- dat$`UniProt Gene Name`[i]
  
  # Format as 1-letter amino acid abbreviation and site number (e.g., Y7)
  position <- dat$`Site List`[i]
  positions <- str_split(position, "; ")[[1]]
  sites <- paste(positions, collapse="|")
  sites_id <- paste(positions, collapse="_")
  
  id <- paste0(gene_symbol, "_", sites_id)
  
  t1 <- data.frame(ID=id, Symbols=gene_symbol,
                   Sites=sites, Effect="",
                   SignedP=sign(dat$logFC[i])*dat$FDR[i], 
                   stringsAsFactors=FALSE)
  
  reformatted_dat <- rbind(reformatted_dat, t1)
}

# The original and reformatted data should have the same row count
stopifnot(nrow(dat) == nrow(reformatted_dat))

# SAVE RESULTS ----
write.table(reformatted_dat, "data_causalpath.txt", sep="\t", 
            row.names=FALSE, quote=FALSE)
