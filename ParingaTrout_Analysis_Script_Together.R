### Bioinformatics Script 2 for: Insights into the ecological impact of trout introduction in an 
### oligotrophic lake using sedimentary environmental DNA 


##################################################################################################################################################################
############     SHORT COI    ############################################################################################################################
#########################################################################################################################


path <- "/home/john/Documents/Georgia/Chapter_2_Plate_4_sequences/CO1"

list.files(path)

###new code

fnFs <- sort(list.files(path, pattern = "R1_001.fastq.gz", full.names = TRUE))
fnRs <- sort(list.files(path, pattern = "R2_001.fastq.gz", full.names = TRUE))

FWD <- "YTCHACWAAYCAYAARGAYATYGG"

REV <- "ARYCARTTHCCRAAHCCHCC"

allOrients <- function(primer) {
  # Create all orientations of the input sequence
  require(Biostrings)
  dna <- DNAString(primer)  # The Biostrings works w/ DNAString objects rather than character vectors
  orients <- c(Forward = dna, Complement = complement(dna), Reverse = reverse(dna), 
               RevComp = reverseComplement(dna))
  return(sapply(orients, toString))  # Convert back to character vector
}
FWD.orients <- allOrients(FWD)
REV.orients <- allOrients(REV)
FWD.orients

primerHits <- function(primer, fn) {
  # Counts number of reads in which the primer is found
  nhits <- vcountPattern(primer, sread(readFastq(fn)), fixed = FALSE)
  return(sum(nhits > 0))
}
rbind(FWD.ForwardReads = sapply(FWD.orients, primerHits, fn = fnFs[[10]]), 
      REV.ReverseReads = sapply(REV.orients, primerHits, fn = fnRs[[10]]))

###cutadapt

cutadapt <- "/home/john/miniconda3/bin/cutadapt"

system2(cutadapt, args = "--version") 


###create path

path.cut <- file.path(path, "cutadapt")
if(!dir.exists(path.cut)) dir.create(path.cut)
fnFs.cut <- file.path(path.cut, basename(fnFs))
fnRs.cut <- file.path(path.cut, basename(fnRs))

# Trim FWD off of R1 (forward reads) - 
R1.flags <- paste0("-g", " ^", FWD) 
# Trim REV off of R2 (reverse reads)
R2.flags <- paste0("-G", " ^", REV) 
# Run Cutadapt
for(i in seq_along(fnFs)) {
  system2(cutadapt, args = c("-e 0.08 --discard-untrimmed", R1.flags, R2.flags,"-m", 1,
                             "-o", fnFs.cut[i], "-p", fnRs.cut[i], # output files
                             fnFs[i], fnRs[i])) # input files
}

rbind(FWD.ForwardReads = sapply(FWD.orients, primerHits, fn = fnFs.cut[[1]]), 
      REV.ReverseReads = sapply(REV.orients, primerHits, fn = fnRs.cut[[1]]))

# Forward and reverse fastq filenames have the format:
cutFs <- sort(list.files(path.cut, pattern = "R1_001.fastq.gz", full.names = TRUE))
cutRs <- sort(list.files(path.cut, pattern = "R2_001.fastq.gz", full.names = TRUE))

if(length(cutFs) == length(cutRs)) print("Forward and reverse files match. Go forth and explore")
if (length(cutFs) != length(cutRs)) stop("Forward and reverse files do not match. Better go back and have a check")

# Extract sample names, assuming filenames have format:
get.sample.name <- function(fname) strsplit(basename(fname), "_")[[1]][1]
sample.names <- unname(sapply(cutFs, get.sample.name))
head(sample.names)

if(length(cutFs) <= 20) {
  fwd_qual_plots <- plotQualityProfile(cutFs) + 
    scale_x_continuous(breaks=seq(0,250,10)) + 
    scale_y_continuous(breaks=seq(0,40,2)) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1))
  rev_qual_plots <- plotQualityProfile(cutRs) + 
    scale_x_continuous(breaks=seq(0,250,10)) + 
    scale_y_continuous(breaks=seq(0,40,2)) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1))
} else {
  rand_samples <- sample(size = 20, 1:length(cutFs)) # grab 20 random samples to plot
  fwd_qual_plots <- plotQualityProfile(cutFs[rand_samples]) + 
    scale_x_continuous(breaks=seq(0,250,10)) + 
    scale_y_continuous(breaks=seq(0,40,2)) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1))
  rev_qual_plots <- plotQualityProfile(cutRs[rand_samples]) + 
    scale_x_continuous(breaks=seq(0,250,10)) + 
    scale_y_continuous(breaks=seq(0,40,2)) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1))
}

fwd_qual_plots
rev_qual_plots

jpeg(file="Quality.Plot.F.jpg",res=300, width=15, height=8, units="in")
fwd_qual_plots
dev.off()

jpeg(file="Quality.Plot.R.jpg",res=300, width=15, height=8, units="in")
rev_qual_plots
dev.off()

####up to here
filtpathF <- file.path(path.cut, "filtered", basename(cutFs))
filtpathR <- file.path(path.cut, "filtered", basename(cutRs))

###for RV have trialled lots of different trunclen and 110,110 is the best

out <- filterAndTrim(cutFs, filtpathF, cutRs, filtpathR,
                     truncLen=c(150,110), maxEE=c(2,2), truncQ=2, maxN=0, rm.phix=TRUE,
                     compress=TRUE, verbose=TRUE, multithread=TRUE)

out2 <- as.data.frame(out)
str(out2)
out2$perc <- (out2$reads.out/out2$reads.in)*100
out2

sample.names <- sapply(strsplit(basename(filtpathF), "_"), `[`, 1) # Assumes filename = samplename_XXX.fastq.gz
sample.namesR <- sapply(strsplit(basename(filtpathR), "_"), `[`, 1) # Assumes filename = samplename_XXX.fastq.gz
if(identical(sample.names, sample.namesR)) {print("Files are still matching.....congratulations")
} else {stop("Forward and reverse files do not match.")}
names(filtpathF) <- sample.names
names(filtpathR) <- sample.namesR

set.seed(100) # set seed to ensure that randomized steps are replicatable

# Learn forward error rates
errF <- learnErrors(filtpathF, nbases=1e8, multithread=TRUE, verbose = TRUE)
## 100285500 total bases in 668570 reads from 24 samples will be used for learning the error rates.

jpeg(file="Error.F.Plot.jpg",res=300, width=15, height=8, units="in")
plotErrors(errF, nominalQ=TRUE)
dev.off()

# Learn reverse error rates
errR <- learnErrors(filtpathR, nbases=1e8, multithread=TRUE, verbose = TRUE)

jpeg(file="Error.R.Plot.jpg",res=300, width=15, height=8, units="in")
plotErrors(errR, nominalQ=TRUE)
dev.off()


derepF <- derepFastq(filtpathF, verbose=TRUE)
derepR <- derepFastq(filtpathR, verbose=TRUE)

dadaF.pseudo <- dada(derepF, err=errF, multithread=TRUE, pool="pseudo")
dadaR.pseudo <- dada(derepR, err=errR, multithread=TRUE, pool="pseudo")

mergers <- mergePairs(dadaF.pseudo, derepF, dadaR.pseudo, derepR, maxMismatch = 1, minOverlap = 10, 
                      verbose=TRUE)

head(mergers)

seqtab <- makeSequenceTable(mergers)

split.dir <- sapply(strsplit(basename(path), "-"), `[`, 1:2)

split.dir.name  <- paste(split.dir[1], split.dir[2], sep=".")

saveRDS(seqtab, "seqtab2.rds")

getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(dadaF.pseudo, getN), sapply(dadaR.pseudo, getN), sapply(mergers, getN))
colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged")
rownames(track) <- sample.names

write.csv(track,"track.csv")

track

###Look at the length distribution of the sequences

trimtable <- as.data.frame(table(nchar(getSequences(seqtab))))
colnames(trimtable) <- c("Length.bp", "Frequency")
trimtable$Frequency <- as.numeric(trimtable$Frequency)
str(trimtable)

g <- ggplot(trimtable, aes(x = Length.bp, y = Frequency)) +
  geom_bar(stat="identity")

jpeg(file="Chim.Dist.Plot.R.jpg",res=300, width=15, height=8, units="in")
g

dev.off()


write.csv(trimtable, "trimtable.csv")



##checking for chimeras (sequences outside the expected size range) Look at table to figure out spread (min,max for amplicon size)

seqtab2 <- seqtab[,nchar(colnames(seqtab)) %in% seq(177,190)]

seqtab.nochim <- removeBimeraDenovo(seqtab2, multithread=TRUE, verbose=TRUE)

###Produce a table which shows the number or reads at each stage

getN <- function(x) sum(getUniques(x))
track2 <- cbind(out, sapply(dadaF.pseudo, getN), sapply(dadaR.pseudo, getN), sapply(mergers, getN), rowSums(seqtab.nochim))
colnames(track2) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
rownames(track2) <- sample.names

write.csv(track2, "track2.csv")  


Tutorial.dada.all.nochimera.forsav <- as.data.frame(seqtab.nochim)
write.csv(Tutorial.dada.all.nochimera.forsav, "seqtab.nochim.csv", row.names = T)


####assign taxonomy #####

tutorial.nochimera.tax <- assignTaxonomy(seqtab.nochim, "BOLD.NCBI.short.genus_dada2.fasta", multithread=TRUE)

##Save cleaned sequences

Tutorial.all.nochimera.tax.forsav <- as.data.frame(tutorial.nochimera.tax)
write.csv(Tutorial.all.nochimera.tax.forsav, "tutorial.nochimera.tax.csv", row.names = T)









################# INPUT FILES ################

## ALL short CO1
seqtab.nochim <- readRDS("Lakes380_seqtab.zooplankton2.nochim.rds")
seqtab.nochim.tax <- readRDS("Georgia_shortCOI_downcore.taxonomy.rds")
Tutorial.all.nochimera.forsav <- as.data.frame(seqtab.nochim)
#write.csv(Tutorial.all.nochimera.forsav, "tutorial.nochimera.csv", row.names = T)
Tutorial.all.nochimera.tax.forsav <- as.data.frame(seqtab.nochim.tax)
#write.csv(Tutorial.all.nochimera.tax.forsav, "tutorial.nochimera.tax.csv", row.names = T)

#####phyloseq#####

####Import sequence table into phyloseq

###otu-table
Tutorial.dada.all.nochimera <- fread("tutorial.nochimera.csv", header = TRUE)
#head(Tutorial.dada.all.nochimera)

otu_table = as.data.frame(Tutorial.dada.all.nochimera)
otu_table2 <- subset(otu_table, select = -c(`V1`) )
row.names(otu_table2) = otu_table$"V1"

ps_otu_table = otu_table(otu_table2, taxa_are_rows = F)
#head(ps_otu_table)

####tax-table
Tutorial.all.nochimera.tax <- fread("tutorial.nochimera.tax.csv", header = TRUE)

tax_table = as.data.frame(Tutorial.all.nochimera.tax)
tax_table2 <- subset(tax_table, select = -c(`V1`) )
row.names(tax_table2) = tax_table$"V1"

ps_tax = tax_table(tax_table2)
#head(ps_tax)
#adding ID's back
rownames(ps_tax)=rownames(tax_table2)
colnames(ps_tax)=colnames(tax_table2)


###get sample data

sample.data <-fread("zooplankton.metadata2.csv")
names(sample.data)

sample.data$Lake <- as.factor(sample.data$Lake)
levels(sample.data$Lake)

sample.data$Type<- as.factor(sample.data$Type)
levels(sample.data$Type)

sample.data$Core<- as.factor(sample.data$Core)
levels(sample.data$Core)

sample.data$Age<- as.factor(sample.data$Age)
levels(sample.data$Age)

ps_map = sample_data(sample.data)
#head(ps_map)
rownames(ps_map)=sample.data$'SampleID'
#head(ps_map)

##new phyloseq

phyloseq <- phyloseq(ps_otu_table, ps_map, ps_tax)

phyloseq

#library(microViz)
phyloseq <- phyloseq %>% ps_arrange(Depth)
#write.csv(tax_table(phyloseq), file="phyloseqtax.csv")

tax.clean <- data.frame(tax_table(phyloseq))

for (i in 1:6){ tax.clean[,i] <- as.character(tax.clean[,i])}

tax.clean[is.na(tax.clean)] <- ""

for (i in 1:nrow(tax.clean)){
  if (tax.clean[i,2] == ""){
    kingdom <- paste("Unclassified_Kingdom_", tax.clean[i,1], sep = "")
    tax.clean[i, 2:6] <- kingdom
  } else if (tax.clean[i,3] == ""){
    phylum <- paste("Unclassified_Phylum_", tax.clean[i,2], sep = "")
    tax.clean[i, 3:6] <- phylum
  } else if (tax.clean[i,4] == ""){
    class <- paste("Unclassified_Class_", tax.clean[i,3], sep = "")
    tax.clean[i, 4:6] <- class
  } else if (tax.clean[i,5] == ""){
    order <- paste("Unclassified_Order_", tax.clean[i,4], sep = "")
    tax.clean[i, 5:6] <- order
  } else if (tax.clean[i,6] == ""){
    tax.clean$Genus[i] <- paste("Unclassified_Family",tax.clean$Family[i], sep = "_")
  } 
}

tax_table(phyloseq) <- as.matrix(tax.clean)
#write.csv(tax_table(phyloseq), file="pholoseqshorttax.csv")


# Make genus IDs unique

tax = tax_table(phyloseq)
taxa_table = as.data.frame(tax)
taxa_table$Genus <- sub('[.]', '_', make.names(taxa_table$Genus, unique=TRUE))
tax.table2 = as.data.frame(taxa_table)
tax = tax_table(taxa_table)

#head(tax)
rownames(tax)=rownames(tax.table2)
colnames(tax)=colnames(tax.table2)
#head(tax)


phyloseq <- phyloseq(ps_otu_table, ps_map, tax)

phylo_paringa = subset_samples(phyloseq, Lake %in% c("Paringa", "Water"))
phylo_paringa1 = prune_taxa(taxa_sums(phylo_paringa) > 0, phylo_paringa)

#extract OTU's that are in blank

####subtraction of contamination

Controls = subset_samples(phylo_paringa1, Type %in% c("Control"))
Controls = filter_taxa(Controls, function(x) sum(x) > 0, TRUE)
sample_sums(Controls)

###removing ASV contamination - ASV total removal
Controls.ASVs <- taxa_names(Controls)

allTaxa = taxa_names(phylo_paringa1)
allTaxa.filtered <- allTaxa[!(allTaxa %in% Controls.ASVs)]
phylo_paringa1.nocont = prune_taxa(allTaxa.filtered, phylo_paringa1)

phylo_paringa1.7 = subset_samples(phylo_paringa1.nocont, Type != "Control")
phylo_paringa1.5 = subset_samples(phylo_paringa1.7, Core != "Paringa SC")

phylo_paringa1.6 = prune_taxa(taxa_sums(phylo_paringa1.5) > 5, phylo_paringa1.5)
ntaxa(phylo_paringa1.6)
# 10837

saveRDS(phylo_paringa1.6, file = "SHORT_phylo_paringa1.6.rds")

### RAREFY

rare = ggrare(phylo_paringa1.6, step = 50,
              plot = TRUE, parallel = FALSE, se = FALSE, color = "Depth", label = "Depth") 
rare + theme_bw() 

set.seed(100)
phylo_paringa1.62 = rarefy_even_depth(phylo_paringa1.6, sample.size = 25000, 
                                      replace = FALSE, trimOTUs = TRUE, verbose = TRUE)

rare = ggrare(phylo_paringa1.62, step = 50,
              plot = TRUE, parallel = FALSE, se = FALSE, color = "Depth", label = "Depth") 
rare + theme_bw() 

##paringa - remove terrestrial groups and bacteria
phylo_paringa1.63 = subset_taxa(phylo_paringa1.62, Kingdom == "Eukaryota")
ntaxa(phylo_paringa1.63)
# # 7792

phylo_paringa3 = subset_taxa(phylo_paringa1.63, (Phylum!= "Chordata")|is.na(Phylum))
phylo_paringa3.1 = subset_taxa(phylo_paringa3, (Phylum!= "Streptophyta")|is.na(Phylum))
phylo_paringa3.5 = subset_taxa(phylo_paringa3.1, (Family!= "Megascolecidae")|is.na(Family))
ntaxa(phylo_paringa3.5)
# 7521
# 1168 unclassified removed
nsamples(phylo_paringa3.5)
# 27

#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~
#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~
#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~

######## FINAL SHORT COI EUKARYOTES
saveRDS(phylo_paringa3.5, file = "UPDATED25000RARE_COIShortEuksFinal.rds")

######## FINAL SHORT COI PREY
phylo_paringaINVER = subset_taxa(phylo_paringa3.5, Phylum %in% c("Arthropoda", "Mollusca",
                                                                 "Rotifera"))
any(taxa_sums(phylo_paringaINVER) == 0)
any(sample_sums(phylo_paringaINVER) == 0)
ntaxa(phylo_paringaINVER)
# 869
nsamples(phylo_paringaINVER)
# 27
saveRDS(phylo_paringaINVER, file = "UPDATED_25000RARE_COIShortInvertsFinal.rds")

######## FINAL SHORT COI ALGAE
phylo_paringaALGAE = subset_taxa(phylo_paringa3.5, Phylum %in% c("Chlorophyta",
                                                                 "Ochrophyta",
                                                                 "Bacillariophyta",
                                                                 "Rhodophyta"))

any(taxa_sums(phylo_paringaALGAE) == 0)
any(sample_sums(phylo_paringaALGAE) == 0)
ntaxa(phylo_paringaALGAE)
# 142
nsamples(phylo_paringaALGAE)
# 27
saveRDS(phylo_paringaALGAE, file = "UPDATED_25000RARE_COIShortAlgaeFinal.rds")

#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~
#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~
#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~

##################################################################################################################################################################
############################# LONG COI #################################################################################################################
##################################################################################################################################################################



library(dada2)
packageVersion("dada2")
library(ShortRead)
packageVersion("ShortRead")
library(Biostrings)
packageVersion("Biostrings")
library(ggplot2)

path <- "~/Documents/Sequencing_results/N2205433_30-708477447_Meta/MS220621-1806/Plate3/COI/"  ## CHANGE ME to the directory containing the fastq files.

list.files(path)


fnFs <- sort(list.files(path, pattern = "1.fq.gz", full.names = TRUE))
fnRs <- sort(list.files(path, pattern = "2.fq.gz", full.names = TRUE))

FWD <- "ACDGGDTGRACHGTNTAYCC"

REV <- "TCDGGRTGNCCRAARAAYCA"

allOrients <- function(primer) {
  # Create all orientations of the input sequence
  require(Biostrings)
  dna <- DNAString(primer)  # The Biostrings works w/ DNAString objects rather than character vectors
  orients <- c(Forward = dna, Complement = complement(dna), Reverse = reverse(dna), 
               RevComp = reverseComplement(dna))
  return(sapply(orients, toString))  # Convert back to character vector
}
FWD.orients <- allOrients(FWD)
REV.orients <- allOrients(REV)
FWD.orients

primerHits <- function(primer, fn) {
  # Counts number of reads in which the primer is found
  nhits <- vcountPattern(primer, sread(readFastq(fn)), fixed = FALSE)
  return(sum(nhits > 0))
}

rbind(FWD.ForwardReads = sapply(FWD.orients, primerHits, fn = fnFs[[1]]), 
      REV.ReverseReads = sapply(REV.orients, primerHits, fn = fnRs[[1]]))


cutadapt <- "/Users/johnpearman/opt/miniconda3/bin/cutadapt" # CHANGE ME to the cutadapt path on your machine
system2(cutadapt, args = "--version")

path.cut <- file.path(path, "cutadapt")
if(!dir.exists(path.cut)) dir.create(path.cut)
fnFs.cut <- file.path(path.cut, basename(fnFs))
fnRs.cut <- file.path(path.cut, basename(fnRs))


# Trim FWD off of R1 (forward reads) - 
R1.flags <- paste0("-g", " ^", FWD) 
# Trim REV off of R2 (reverse reads)
R2.flags <- paste0("-G", " ^", REV) 
# Run Cutadapt
for(i in seq_along(fnFs)) {
  system2(cutadapt, args = c("-e 0.05 --discard-untrimmed", R1.flags, R2.flags,
                             "-o", fnFs.cut[i], "-p", fnRs.cut[i], # output files
                             fnFs[i], fnRs[i])) # input files
}
  
rbind(FWD.ForwardReads = sapply(FWD.orients, primerHits, fn = fnFs.cut[[1]]), 
        REV.ReverseReads = sapply(REV.orients, primerHits, fn = fnRs.cut[[1]]))
  
  
cutFs <- sort(list.files(path.cut, pattern = "1.fq.gz", full.names = TRUE))
cutRs <- sort(list.files(path.cut, pattern = "2.fq.gz", full.names = TRUE))
  
  
if(length(cutRs) == length(cutRs)) print("Forward and reverse files match. Go forth and explore")
if (length(cutRs) != length(cutRs)) stop("Forward and reverse files do not match. Better go back and have a check")
  
  
  # Extract sample names, assuming filenames have format:
get.sample.name <- function(fname) strsplit(basename(fname), "_")[[1]][1]
sample.names <- unname(sapply(cutFs, get.sample.name))
head(sample.names)
  
  
  if(length(cutFs) <= 20) {
    fwd_qual_plots <- plotQualityProfile(cutFs) + 
      scale_x_continuous(breaks=seq(0,250,10)) + 
      scale_y_continuous(breaks=seq(0,40,2)) +
      theme(axis.text.x = element_text(angle = 90, hjust = 1))
    rev_qual_plots <- plotQualityProfile(cutRs) + 
      scale_x_continuous(breaks=seq(0,250,10)) + 
      scale_y_continuous(breaks=seq(0,40,2)) +
      theme(axis.text.x = element_text(angle = 90, hjust = 1))
  } else {
    rand_samples <- sample(size = 20, 1:length(cutFs)) # grab 20 random samples to plot
    fwd_qual_plots <- plotQualityProfile(cutFs[rand_samples]) + 
      scale_x_continuous(breaks=seq(0,250,10)) + 
      scale_y_continuous(breaks=seq(0,40,2)) +
      theme(axis.text.x = element_text(angle = 90, hjust = 1))
    rev_qual_plots <- plotQualityProfile(cutRs[rand_samples]) + 
      scale_x_continuous(breaks=seq(0,250,10)) + 
      scale_y_continuous(breaks=seq(0,40,2)) +
      theme(axis.text.x = element_text(angle = 90, hjust = 1))
}
  
fwd_qual_plots
rev_qual_plots
  
jpeg(file="Quality.Plot.F.jpg",res=300, width=15, height=8, units="in")
fwd_qual_plots
dev.off()
  
  
jpeg(file="Quality.Plot.R.jpg",res=300, width=15, height=8, units="in")
rev_qual_plots
dev.off()
  
  
filtpathF <- file.path(path.cut, "filtered", basename(cutFs))
filtpathR <- file.path(path.cut, "filtered", basename(cutRs))
  
out <- filterAndTrim(cutFs, filtpathF, cutRs, filtpathR,
                     truncLen=c(228,230), maxEE=c(4,6), truncQ=2, maxN=0, rm.phix=TRUE,
                     compress=TRUE, verbose=TRUE, multithread=TRUE)

out

sample.names <- sapply(strsplit(basename(filtpathF), "_"), `[`, 1) # Assumes filename = samplename_XXX.fastq.gz
sample.namesR <- sapply(strsplit(basename(filtpathR), "_"), `[`, 1) # Assumes filename = samplename_XXX.fastq.gz
if(identical(sample.names, sample.namesR)) {print("Files are still matching.....congratulations")
} else {stop("Forward and reverse files do not match.")}
names(filtpathF) <- sample.names
names(filtpathR) <- sample.namesR

library(magrittr)

loessErrfun_mod <- function (trans) {
  qq <- as.numeric(colnames(trans))
  est <- matrix(0, nrow = 0, ncol = length(qq))
  for (nti in c("A", "C", "G", "T")) {
    for (ntj in c("A", "C", "G", "T")) {
      if (nti != ntj) {
        errs <- trans[paste0(nti, "2", ntj), ]
        tot <- colSums(trans[paste0(nti, "2", c("A",
                                                "C", "G", "T")), ])
        rlogp <- log10((errs + 1)/tot)
        rlogp[is.infinite(rlogp)] <- NA
        df <- data.frame(q = qq, errs = errs, tot = tot,
                         rlogp = rlogp)
        mod.lo <- loess(rlogp ~ q, df, weights = log10(tot),span = 2)
        pred <- predict(mod.lo, qq)
        maxrli <- max(which(!is.na(pred)))
        minrli <- min(which(!is.na(pred)))
        pred[seq_along(pred) > maxrli] <- pred[[maxrli]]
        pred[seq_along(pred) < minrli] <- pred[[minrli]]
        est <- rbind(est, 10^pred)
      }
    }
  }
  MAX_ERROR_RATE <- 0.25
  MIN_ERROR_RATE <- 1e-07
  est[est > MAX_ERROR_RATE] <- MAX_ERROR_RATE
  est[est < MIN_ERROR_RATE] <- MIN_ERROR_RATE
  err <- rbind(1 - colSums(est[1:3, ]), est[1:3, ], est[4,
  ], 1 - colSums(est[4:6, ]), est[5:6, ], est[7:8, ], 1 -
    colSums(est[7:9, ]), est[9, ], est[10:12, ], 1 - colSums(est[10:12,
    ]))
  rownames(err) <- paste0(rep(c("A", "C", "G", "T"), each = 4),
                          "2", c("A", "C", "G", "T"))
  colnames(err) <- colnames(trans)
  return(err)
}


set.seed(100) # set seed to ensure that randomized steps are replicatable

# Learn forward error rates
errF1 <- learnErrors(filtpathF, nbases=1e8, errorEstimationFunction = loessErrfun_mod, multithread=TRUE, verbose = TRUE)

errR1 <- learnErrors(filtpathR, nbases=1e8, errorEstimationFunction = loessErrfun_mod, multithread=TRUE, verbose = TRUE)

errF_plot <- plotErrors(errF1, nominalQ=TRUE)
errF_plot

errR_plot <- plotErrors(errR1, nominalQ=TRUE)
errR_plot

derepF <- derepFastq(filtpathF, verbose=TRUE)
derepR <- derepFastq(filtpathR, verbose=TRUE)

dadaF.pseudo <- dada(derepF, err=errF1, multithread=TRUE, pool="pseudo")
dadaR.pseudo <- dada(derepR, err=errR1, multithread=TRUE, pool="pseudo")

mergers <- mergePairs(dadaF.pseudo, derepF, dadaR.pseudo, derepR, maxMismatch = 1, minOverlap = 10, verbose=TRUE)
seqtab <- makeSequenceTable(mergers)

saveRDS(seqtab, "Paringa_DC_coi_seqtab.rds")

getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(dadaF.pseudo, getN), sapply(dadaR.pseudo, getN), sapply(mergers, getN))
colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged")
rownames(track) <- sample.names

write.csv(track, "Paringa_DC_coi_track.csv")

rownames(seqtab) <- c("PAR13-5BF", "PAR15-5BF", "PAR22BF", "PAR_SC7_BF", "PAR_SC9_BF", "PAR5-5BF", "PAR7-5BF", "PAR11-5BF", "PAR24BF", "PAR41BF", "PAR44-5BF", "PAR45-5BF", "PAR27BF", "PAR30-5BF", "PAR32-5BF", "PAR34-5BF", "PAR36BF", "PAR37BF", "PAR39BF", "PAR40BF", "PAR46-5BF", "PAR85-5BF", "Par-B1BF", "Par-B2BF", "PAR47-5BF", "PAR48-5BF", "PAR54BF", "PAR55BF", "PAR60BF", "PAR64BF", "PAR70-5BF", "PAR78BF")

table(nchar(getSequences(seqtab)))
hist(nchar(getSequences(seqtab)), main="Distribution of sequence lengths")

seqtab.2 <- seqtab[,nchar(colnames(seqtab)) %in% seq(316,316)]

seqtab.nochim <- removeBimeraDenovo(seqtab.2, multithread=TRUE, verbose=TRUE)

saveRDS(seqtab.nochim, "Paringa_DC_seqtab.nochim.rds")

taxonomy <- assignTaxonomy(seqtab.nochim, "~/Documents/ReferenceDB/COI.database.fasta", minBoot = 50, multithread=T)

saveRDS(taxonomy, "Paringa_DC_seqtab.nochim.tax.rds")



################# INPUT FILES ################

## ALL long CO1

seqtab.nochim <- readRDS("Paringa_DC_seqtab.nochim.rds")
seqtab.nochim.tax <- readRDS("Paringa_DC_seqtab.nochim.tax.rds")

Tutorial.all.nochimera.forsav <- as.data.frame(seqtab.nochim)
#write.csv(Tutorial.all.nochimera.forsav, "tutorial.nochimeraBF.csv", row.names = T)
Tutorial.all.nochimera.tax.forsav <- as.data.frame(seqtab.nochim.tax)
#write.csv(Tutorial.all.nochimera.tax.forsav, "tutorial.nochimeraBF.tax.csv", row.names = T)

#####phyloseq#####

####Import sequence table into phyloseq

###otu-table
Tutorial.dada.all.nochimera <- fread("tutorial.nochimeraBF.csv", header = TRUE)
#head(Tutorial.dada.all.nochimera)

otu_table = as.data.frame(Tutorial.dada.all.nochimera)
otu_table2 <- subset(otu_table, select = -c(`V1`) )
row.names(otu_table2) = otu_table$"V1"

ps_otu_table = otu_table(otu_table2, taxa_are_rows = F)
#head(ps_otu_table)

####tax-table
Tutorial.all.nochimera.tax <- fread("tutorial.nochimeraBF.tax.csv", header = TRUE)

tax_table = as.data.frame(Tutorial.all.nochimera.tax)
tax_table2 <- subset(tax_table, select = -c(`V1`) )
row.names(tax_table2) = tax_table$"V1"

ps_tax = tax_table(tax_table2)
#head(ps_tax)
#adding ID's back
rownames(ps_tax)=rownames(tax_table2)
colnames(ps_tax)=colnames(tax_table2)


###get sample data

sample.data <-fread("Sample data Bact.BF.csv")
names(sample.data)

sample.data$Lake <- as.factor(sample.data$Lake)
levels(sample.data$Lake)

sample.data$Type<- as.factor(sample.data$Type)
levels(sample.data$Type)

sample.data$Core<- as.factor(sample.data$Core)
levels(sample.data$Core)

sample.data$Age<- as.factor(sample.data$Age)
levels(sample.data$Age)

ps_map = sample_data(sample.data)
#head(ps_map)
rownames(ps_map)=sample.data$'SampleID'
#head(ps_map)

##new phyloseq

phyloseq <- phyloseq(ps_otu_table, ps_map, ps_tax)

phyloseq



tax.clean <- data.frame(tax_table(phyloseq))

for (i in 1:6){ tax.clean[,i] <- as.character(tax.clean[,i])}

tax.clean[is.na(tax.clean)] <- ""

for (i in 1:nrow(tax.clean)){
  if (tax.clean[i,2] == ""){
    kingdom <- paste("Unclassified_Kingdom_", tax.clean[i,1], sep = "")
    tax.clean[i, 2:6] <- kingdom
  } else if (tax.clean[i,3] == ""){
    phylum <- paste("Unclassified_Phylum_", tax.clean[i,2], sep = "")
    tax.clean[i, 3:6] <- phylum
  } else if (tax.clean[i,4] == ""){
    class <- paste("Unclassified_Class_", tax.clean[i,3], sep = "")
    tax.clean[i, 4:6] <- class
  } else if (tax.clean[i,5] == ""){
    order <- paste("Unclassified_Order_", tax.clean[i,4], sep = "")
    tax.clean[i, 5:6] <- order
  } else if (tax.clean[i,6] == ""){
    tax.clean$Genus[i] <- paste("Unclassified_Family",tax.clean$Family[i], sep = "_")
  } 
}

tax_table(phyloseq) <- as.matrix(tax.clean)
#write.csv(tax_table(phyloseq), file="phyloseqtax.csv")


# Make genus IDs unique

tax = tax_table(phyloseq)
taxa_table = as.data.frame(tax)
taxa_table$Genus <- sub('[.]', '_', make.names(taxa_table$Genus, unique=TRUE))
tax.table2 = as.data.frame(taxa_table)
tax = tax_table(taxa_table)

#head(tax)
rownames(tax)=rownames(tax.table2)
colnames(tax)=colnames(tax.table2)
#head(tax)


phyloseq <- phyloseq(ps_otu_table, ps_map, tax)

phylo_paringa = subset_samples(phyloseq, Lake %in% c("Paringa", "Water"))
phylo_paringa1 = prune_taxa(taxa_sums(phylo_paringa) > 0, phylo_paringa)

#extract OTU's that are in blank

####subtraction of contamination

Controls = subset_samples(phylo_paringa1, Type %in% c("Control"))
Controls = filter_taxa(Controls, function(x) sum(x) > 0, TRUE)
sample_sums(Controls)

###removing ASV contamination - ASV total removal
Controls.ASVs <- taxa_names(Controls)

allTaxa = taxa_names(phylo_paringa1)
allTaxa.filtered <- allTaxa[!(allTaxa %in% Controls.ASVs)]
phylo_paringa1.nocont = prune_taxa(allTaxa.filtered, phylo_paringa1)

phylo_paringa1.7 = subset_samples(phylo_paringa1.nocont, Type != "Control")
phylo_paringa1.5 = subset_samples(phylo_paringa1.7, Core != "Paringa SC")

phylo_paringa1.6 = prune_taxa(taxa_sums(phylo_paringa1.5) > 5, phylo_paringa1.5)
ntaxa(phylo_paringa1.6)
# 4054

saveRDS(phylo_paringa1.6, file = "LONG_phylo_paringa1.6.rds")


rare = ggrare(phylo_paringa1.6, step = 50,
              plot = TRUE, parallel = FALSE, se = FALSE, color = "Depth", label = "Depth") 
rare + theme_bw() 

set.seed(150)
phylo_paringa1.62 = rarefy_even_depth(phylo_paringa1.6, sample.size = 4000, 
                                      replace = FALSE, trimOTUs = TRUE, verbose = TRUE)

rare = ggrare(phylo_paringa1.62, step = 50,
              plot = TRUE, parallel = FALSE, se = FALSE, color = "Depth", label = "Depth") 
rare + theme_bw() 

##paringa - remove terrestrial groups and bacteria
phylo_paringa1.63 = subset_taxa(phylo_paringa1.62, Kingdom == "Eukaryota")
ntaxa(phylo_paringa1.63)
# 4002

# 1704
phylo_paringa3 = subset_taxa(phylo_paringa1.63, (Phylum!= "Chordata")|is.na(Phylum))
phylo_paringa3.1 = subset_taxa(phylo_paringa3, (Phylum!= "Streptophyta")|is.na(Phylum))
phylo_paringa3.5 = subset_taxa(phylo_paringa3.1, (Family!= "Megascolecidae")|is.na(Family))
ntaxa(phylo_paringa3.5)
# 3948
# remove unclassified 1652
nsamples(phylo_paringa3.5)
#  28

#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~
#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~
#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~

######## FINAL LONG COI EUKARYOTES
saveRDS(phylo_paringa3.5, file = "COILongEuksFinal.rds")
COILongEuksFinal = readRDS("COILongEuksFinal.rds")

######## FINAL LONG COI PREY
long_paringaINVER = subset_taxa(COILongEuksFinal, Phylum %in% c("Arthropoda", "Mollusca",
                                                                "Rotifera"))
any(taxa_sums(long_paringaINVER) == 0)
any(sample_sums(long_paringaINVER) == 0)
ntaxa(long_paringaINVER)
# 536
nsamples(long_paringaINVER)
# 28
saveRDS(long_paringaINVER, file = "COILongInvertsFinal.rds")
write.csv(tax_table(long_paringaINVER), file="long_paringaINVERtax.csv")

######## FINAL LONG COI ALGAE
long_paringaALGAE = subset_taxa(phylo_paringa3.5, Phylum %in% c("Chlorophyta",
                                                                "Ochrophyta",
                                                                "Bacillariophyta",
                                                                "Rhodophyta"))

any(taxa_sums(long_paringaALGAE) == 0)
any(sample_sums(long_paringaALGAE) == 0)
ntaxa(long_paringaALGAE)
# 83
nsamples(long_paringaALGAE)
# 28
saveRDS(long_paringaALGAE, file = "COILongAlgaeFinal.rds")

#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~
#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~
#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~


####################################################################################################################################################################
############### 18S ######################################################################################################################
#####################################################################################################################################################################

library(dada2)
packageVersion("dada2")
library(ShortRead)
packageVersion("ShortRead")
library(Biostrings)
packageVersion("Biostrings")
library(ggplot2)

path <- "~/Documents/Sequencing_results/AG0246-89-313806324/18S/"
list.files(path)

fnFs <- sort(list.files(path, pattern = "R1_001.fastq.gz", full.names = TRUE))
fnRs <- sort(list.files(path, pattern = "R2_001.fastq.gz", full.names = TRUE))

FWD <- "AGGGCAAKYCTGGTGCCAGC"

REV <- "RCGGTATCTRATCGYCTT"


allOrients <- function(primer) {
  # Create all orientations of the input sequence
  require(Biostrings)
  dna <- DNAString(primer)  # The Biostrings works w/ DNAString objects rather than character vectors
  orients <- c(Forward = dna, Complement = complement(dna), Reverse = reverse(dna), 
               RevComp = reverseComplement(dna))
  return(sapply(orients, toString))  # Convert back to character vector
}
FWD.orients <- allOrients(FWD)
REV.orients <- allOrients(REV)
FWD.orients

primerHits <- function(primer, fn) {
  # Counts number of reads in which the primer is found
  nhits <- vcountPattern(primer, sread(readFastq(fn)), fixed = FALSE)
  return(sum(nhits > 0))
}
rbind(FWD.ForwardReads = sapply(FWD.orients, primerHits, fn = fnFs[[5]]), 
      REV.ReverseReads = sapply(REV.orients, primerHits, fn = fnRs[[5]]))


cutadapt <- "/Users/johnpearman/Documents/miniconda2/bin/cutadapt" # CHANGE ME to the cutadapt path on your machine
system2(cutadapt, args = "--version")


path.cut <- file.path(path, "cutadapt18S")
if(!dir.exists(path.cut)) dir.create(path.cut)
fnFs.cut <- file.path(path.cut, basename(fnFs))
fnRs.cut <- file.path(path.cut, basename(fnRs))


# Trim FWD off of R1 (forward reads) - 
R1.flags <- paste0("-g", " ^", FWD) 
# Trim REV off of R2 (reverse reads)
R2.flags <- paste0("-G", " ^", REV) 
# Run Cutadapt
for(i in seq_along(fnFs)) {
  system2(cutadapt, args = c("-e 0.08 --discard-untrimmed", R1.flags, R2.flags,
                             "-o", fnFs.cut[i], "-p", fnRs.cut[i], # output files
                             fnFs[i], fnRs[i])) # input files
}



rbind(FWD.ForwardReads = sapply(FWD.orients, primerHits, fn = fnFs.cut[[1]]), 
      REV.ReverseReads = sapply(REV.orients, primerHits, fn = fnRs.cut[[1]]))

cutFs <- sort(list.files(path.cut, pattern = "R1_001.fastq.gz", full.names = TRUE))
cutRs <- sort(list.files(path.cut, pattern = "R2_001.fastq.gz", full.names = TRUE))


if(length(cutRs) == length(cutRs)) print("Forward and reverse files match. Go forth and explore")
if (length(cutRs) != length(cutRs)) stop("Forward and reverse files do not match. Better go back and have a check")


# Extract sample names, assuming filenames have format:
get.sample.name <- function(fname) strsplit(basename(fname), "_")[[1]][1]
sample.names <- unname(sapply(cutFs, get.sample.name))
head(sample.names)

if(length(cutFs) <= 30) {
  fwd_qual_plots <- plotQualityProfile(cutFs) + 
    scale_x_continuous(breaks=seq(0,250,10)) + 
    scale_y_continuous(breaks=seq(0,40,2)) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1))
  rev_qual_plots <- plotQualityProfile(cutRs) + 
    scale_x_continuous(breaks=seq(0,250,10)) + 
    scale_y_continuous(breaks=seq(0,40,2)) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1))
} else {
  rand_samples <- sample(size = 30, 1:length(cutFs)) # grab 20 random samples to plot
  fwd_qual_plots <- plotQualityProfile(cutFs[rand_samples]) + 
    scale_x_continuous(breaks=seq(0,250,10)) + 
    scale_y_continuous(breaks=seq(0,40,2)) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1))
  rev_qual_plots <- plotQualityProfile(cutRs[rand_samples]) + 
    scale_x_continuous(breaks=seq(0,250,10)) + 
    scale_y_continuous(breaks=seq(0,40,2)) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1))
}

fwd_qual_plots
rev_qual_plots

jpeg(file="Quality.Plot.18S.F.jpg",res=300, width=15, height=8, units="in")
fwd_qual_plots
dev.off()

jpeg(file="Quality.Plot.18S.R.jpg",res=300, width=15, height=8, units="in")
rev_qual_plots
dev.off()


filtpathF <- file.path(path.cut, "filtered", basename(cutFs))
filtpathR <- file.path(path.cut, "filtered", basename(cutRs))


out <- filterAndTrim(cutFs, filtpathF, cutRs, filtpathR,
                     truncLen=c(226,228), maxEE=c(2,4), truncQ=2, maxN=0, rm.phix=TRUE,
                     compress=TRUE, verbose=TRUE, multithread=TRUE)

out


sample.names <- sapply(strsplit(basename(filtpathF), "_"), `[`, 1) # Assumes filename = samplename_XXX.fastq.gz
sample.namesR <- sapply(strsplit(basename(filtpathR), "_"), `[`, 1) # Assumes filename = samplename_XXX.fastq.gz
if(identical(sample.names, sample.namesR)) {print("Files are still matching.....congratulations")
} else {stop("Forward and reverse files do not match.")}
names(filtpathF) <- sample.names
names(filtpathR) <- sample.namesR

set.seed(100) # set seed to ensure that randomized steps are replicatable

# Learn forward error rates
errF <- learnErrors(filtpathF, nbases=1e8, multithread=TRUE, verbose = TRUE)
## 100285500 total bases in 668570 reads from 24 samples will be used for learning the error rates.

# Learn reverse error rates
errR <- learnErrors(filtpathR, nbases=1e8, multithread=TRUE, verbose = TRUE)


errF_plot <- plotErrors(errF, nominalQ=TRUE)
errF_plot

errR_plot <- plotErrors(errR, nominalQ=TRUE)
errR_plot


derepF <- derepFastq(filtpathF, verbose=TRUE)
derepR <- derepFastq(filtpathR, verbose=TRUE)

dadaF.pseudo <- dada(derepF, err=errF, multithread=TRUE, pool="pseudo")
dadaR.pseudo <- dada(derepR, err=errR, multithread=TRUE, pool="pseudo")

mergers <- mergePairs(dadaF.pseudo, derepF, dadaR.pseudo, derepR, maxMismatch = 1, minOverlap = 10, verbose=TRUE)
seqtab <- makeSequenceTable(mergers)

saveRDS(seqtab, "Katie.2022.18S.rds")

getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(dadaF.pseudo, getN), sapply(dadaR.pseudo, getN), sapply(mergers, getN))
colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged")
rownames(track) <- sample.names

write.csv(track, "Zooplankton18S.track.csv")


Katie.2022.18S <- readRDS("Katie.2022.18S.rds")

table(nchar(getSequences(Katie.2022.18S)))

hist(nchar(getSequences(Katie.2022.18S)), main="Distribution of sequence lengths")

Katie.2022.18S <- Katie.2022.18S[,nchar(colnames(Katie.2022.18S)) %in% seq(400,443)]

Katie.2022.18S.nochim <- removeBimeraDenovo(Katie.2022.18S, multithread=TRUE, verbose=TRUE)

saveRDS(Katie.2022.18S.nochim, "Katie.2022.18S.nochim.rds")


taxonomy <- assignTaxonomy(seqtab.nochim, "~/Documents/ReferenceDB/pr2_v14.fasta", minBoot = 70, multithread=T)

saveRDS(taxonomy, "Katie.2022.18S.tax.rds")


Katie.18S.otutab <- readRDS("Katie.2022.18S.nochim.rds")
Katie.18S.tax <- readRDS("Katie.2022.18S.tax.rds")
map <- read.csv("metadata.csv", h=T, row.names = 1)


Katie.18S.ps <- phyloseq(otu_table(Katie.18S.otutab, taxa_are_rows=FALSE),
                         sample_data(map), 
                         tax_table(Katie.18S.tax))


Batch1_neg = subset_samples(Katie.18S.ps, Type == "Control")

Batch1_neg_max <- apply(data.frame(as.matrix(t(otu_table(Batch1_neg)))), 1, max)

Batch1_neg_max_vec <- as.vector(Batch1_neg_max)

B1 = as(otu_table(Katie.18S.ps), "matrix")
B1df = as.data.frame(B1)

B1df[,1:length(B1df)] <- sweep(B1df[,1:length(B1df)],2,Batch1_neg_max_vec)


B1df <- replace(B1df, B1df < 0, 0)

B1df_noneg <- B1df[rowSums(B1df)!=0, ]


Katie.2022.18S.cleaned.ps <- phyloseq(otu_table(B1df_noneg, taxa_are_rows=FALSE), 
                                      sample_data(map), 
                                      tax_table(Katie.18S.tax))



org.ss <- as.data.frame(sample_sums(Katie.18S.ps))
org.ss$Names <- rownames(org.ss)

new.ss <- as.data.frame(sample_sums(Katie.2022.18S.cleaned.ps))
new.ss$Names <- rownames(new.ss)

control.rem.sums <- dplyr::left_join(org.ss, new.ss, by="Names")

colnames(control.rem.sums) <- c("Original", "Names", "New")

control.rem.sums <- control.rem.sums %>%
  mutate(perc = New/Original*100)

control.rem.sums


tax.clean <- data.frame(tax_table(Katie.2022.18S.cleaned.ps))

for (i in 1:6){ tax.clean[,i] <- as.character(tax.clean[,i])}

tax.clean[is.na(tax.clean)] <- ""

for (i in 1:nrow(tax.clean)){
  if (tax.clean[i,2] == ""){
    kingdom <- paste("Unclassified_Kingdom_", tax.clean[i,1], sep = "")
    tax.clean[i, 2:6] <- kingdom
  } else if (tax.clean[i,3] == ""){
    phylum <- paste("Unclassified_Phylum_", tax.clean[i,2], sep = "")
    tax.clean[i, 3:6] <- phylum
  } else if (tax.clean[i,4] == ""){
    class <- paste("Unclassified_Class_", tax.clean[i,3], sep = "")
    tax.clean[i, 4:6] <- class
  } else if (tax.clean[i,5] == ""){
    order <- paste("Unclassified_Order_", tax.clean[i,4], sep = "")
    tax.clean[i, 5:6] <- order
  } else if (tax.clean[i,6] == ""){
    tax.clean$Genus[i] <- paste("Unclassified_Family",tax.clean$Family[i], sep = "_")
  } 
}

tax_table(Katie.2022.18S.cleaned.ps) <- as.matrix(tax.clean)





#### Paringa Downcore eDNA Eukaryotic 18S
# These are phyloseq input files from Katie.2022.18S.cleaned.ps.rds put into the depth sequence order, and 
# with before/after trout added to sample data

taxW = fread("/Users/Lena/Desktop/OneDrive/Documents/Hydrosphere/2021/Lakes380 Cawthron 2021/Analysis/18S Euks/TaxW.csv", header = TRUE)
taxW = as.data.frame(taxW)
row.names(taxW)=taxW$Sample
taxW$`Sample`=NULL
tax = as.matrix(taxW)
taxa = tax_table(tax)

#write.csv(otu_table(Katie.2022.18S.cleaned.ps), file="C:/Users/Lena/OneDrive - University of Otago/Picocyanobacteria/Hydrosphere/Lakes380 Cawthron 2021/Analysis/18S Euks/OTU_ar.csv")
otuW = fread("/Users/Lena/Desktop/OneDrive/Documents/Hydrosphere/2021/Lakes380 Cawthron 2021/Analysis/18S Euks/OTU_ar.csv", header = TRUE)
otuW = as.data.frame(otuW)
row.names(otuW)=otuW$"V1"
otuW$`V1`=NULL
otutab = as.matrix(otuW)
otutabl = otu_table(otutab, taxa_are_rows = F)

#write.csv(sample_data(Katie.2022.18S.cleaned.ps), file="C:/Users/Lena/OneDrive - University of Otago/Picocyanobacteria/Hydrosphere/Lakes380 Cawthron 2021/Analysis/Bacteria/fixedtax1_Sam.csv")
sam_data = fread("/Users/Lena/Desktop/OneDrive/Documents/Hydrosphere/2021/Lakes380 Cawthron 2021/Analysis/18S Euks/sam_ar.csv", header = TRUE)
sam_data = as.data.frame(sam_data)
row.names(sam_data) = sam_data$"SampleID" #only now possible to replace the rownnames with the first column (OTU ID)
#sam_data$`SampleID`=NULL #the first column "OTU-ID can then be deleted"

ps_map = sample_data(sam_data)

fixedtax <- merge_phyloseq(otutabl, taxa, ps_map) 

ntaxa(fixedtax)
# 18988


# renaming

tax.clean <- data.frame(tax_table(fixedtax))

for (i in 1:6){ tax.clean[,i] <- as.character(tax.clean[,i])}

tax.clean[is.na(tax.clean)] <- ""

for (i in 1:nrow(tax.clean)){
  if (tax.clean[i,2] == ""){
    kingdom <- paste("Unclassified_Kingdom_", tax.clean[i,1], sep = "")
    tax.clean[i, 2:6] <- kingdom
  } else if (tax.clean[i,6] == ""){
    tax.clean$Genus[i] <- paste("",tax.clean$Family[i], sep = "")
  } 
}

tax_table(fixedtax) <- as.matrix(tax.clean)

# Make strain IDs unique

tax = tax_table(fixedtax)
taxa_table = as.data.frame(tax)
taxa_table$Genus <- sub('[.]', '_', make.names(taxa_table$Genus, unique=TRUE))
tax.table2 = as.data.frame(taxa_table)
tax = tax_table(taxa_table)

#head(tax)
rownames(tax)=rownames(tax.table2)
colnames(tax)=colnames(tax.table2)
#head(tax)

downcore_18Sfixed = phyloseq(otu_table(fixedtax), tax, sample_data(fixedtax))
ntaxa(downcore_18Sfixed)
# 18988 at ASV with everything in
nsamples(downcore_18Sfixed)
# 58

#saveRDS(downcore_18Sfixed, file = "FINALdowncore_18Sfixed.rds")
#downcore_18Sfixed = readRDS("FINALdowncore_18Sfixed.rds")

desired_order <- list("PARIN-SC3M-1A-2-5","PARIN-SC3M-1A-3-5","PARIN-SC3M-1A-7",
                      "PARIN-SC3M-1A-9","PARIN-SC3M-1A-11-5","PARIN-1M-1A-5","PARIN-1M-1A-5-5","PARIN-1M-1A-7-5",
                      "PARIN-1M-1A-11-5","PARIN-1M-1A-12-5","PARIN-1M-1A-13-5","PARIN-1M-1A-14-5","PARIN-1M-1A-15-5",
                      "PARIN-1M-1A-20-5","PARIN-1M-1A-22","PARIN-1M-1A-24","PARIN-1M-1A-25","PARIN-1M-1A-26",
                      "PARIN-1M-1A-27","PARIN-1M-1A-29-5","PARIN-1M-1A-30-5","PARIN-1M-1A-31-5","PARIN-1M-1A-32-5",
                      "PARIN-1M-1A-34-5","PARIN-1M-1A-36","PARIN-1M-1A-37","PARIN-1M-1A-39","PARIN-1M-1A-40",
                      "PARIN-1M-1A-41","PARIN-1M-1A-44-5","PARIN-1M-1A-45-5","PARIN-1M-1A-46-5","PARIN-1M-1A-47-5",
                      "PARIN-1M-1A-48-5","PARIN-1M-1A-54","PARIN-1M-1A-55","PARIN-1M-1A-60","PARIN-1M-1A-64",
                      "PARIN-1M-1A-70-5","PARIN-1M-1A-78","PARIN-1M-1A-85-5","PARIN-1M-1A-88","PARIN-1M-1A-91",
                      "PARIN-1M-1A-94","PARIN-1M-2B-100","PARIN-1M-2B-104","PARIN-1M-2B-106","PARIN-1M-2B-108",
                      "PARIN-1M-2B-110","PARIN-1M-2B-114-5","PARIN-1M-2B-163","PARIN-1M-2B-168-5",
                      "PARIN-1M-2B-173-5","PARIN-1M-2B-183","PARIN-1M-3A-212","PARIN-1M-3A-219-5",
                      "PARIN-1M-3A-228","PARIN-1M-3A-237-5")

# remove samples after 100cm depth

Samples_toRemove <- c("PARIN-SC3M-1A-2-5","PARIN-SC3M-1A-3-5","PARIN-SC3M-1A-7",
                      "PARIN-SC3M-1A-9","PARIN-SC3M-1A-11-5", "PARIN-1M-2B-100","PARIN-1M-2B-104","PARIN-1M-2B-106","PARIN-1M-2B-108",
                      "PARIN-1M-2B-110","PARIN-1M-2B-114-5","PARIN-1M-2B-163","PARIN-1M-2B-168-5",
                      "PARIN-1M-2B-173-5","PARIN-1M-2B-183","PARIN-1M-3A-212","PARIN-1M-3A-219-5",
                      "PARIN-1M-3A-228","PARIN-1M-3A-237-5")
#To see what samples get removed, run the following; note, I have a column called "SampleID"
subset_samples(downcore_18Sfixed, SampleID %in% Samples_toRemove)
#This will return a ps object that contains the samples you want to remove

#To remove those from your phyloseq object
downcore_18S_NEW = subset_samples(downcore_18Sfixed, !(SampleID %in% Samples_toRemove))
#This will return a ps object with the samples removed

######## PRUNE LESS THAN 5

ntaxa(downcore_18S_NEW)
# 18988

downcore_18S_NEWpruned = prune_taxa(taxa_sums(downcore_18S_NEW) > 5, downcore_18S_NEW)
ntaxa(downcore_18S_NEWpruned)
# 3498
nsamples(downcore_18S_NEWpruned)
# 39

###### Rarefy

rare = ggrare(downcore_18S_NEWpruned, step = 50,
              plot = TRUE, parallel = FALSE, se = FALSE, color = "Depth", label = "Depth") 
rare + theme_bw() 

rare20000_18s=rarefy_even_depth(downcore_18S_NEWpruned, sample.size = 20000, rngseed=FALSE, replace=F, trimOTUs = TRUE, verbose = TRUE)
ntaxa(rare20000_18s)
#3436

rare = ggrare(rare20000_18s, step = 50,
              plot = TRUE, parallel = FALSE, se = FALSE, color = "Depth", label = "Depth") 
rare + theme_bw()


######################################

### Add correct depth to sample data here:

sam_data1 = fread("/Users/Lena/OneDrive - University of Otago/Picocyanobacteria/Hydrosphere/Lakes380 Cawthron 2021/Analysis/18S Euks/Samdata_Aquatic_ASVsFixed.csv", header = TRUE)
sam_data1 = as.data.frame(sam_data1)
row.names(sam_data1) = sam_data1$"SampleID" #only now possible to replace the rownnames with the first column (OTU ID)
ps_map1 = sample_data(sam_data1)

######### --------------- CREATE THE PHYLOSEQ FILE ------------------- #########

FINAL_rare20000_18sFIXED <- merge_phyloseq(otu_table(rare20000_18s), tax_table(rare20000_18s), ps_map1) 
saveRDS(FINAL_rare20000_18sFIXED, file = "FINAL_rare20000_18sFIXED.rds")

####### adding sample interval to age model

sam_data1 = fread("C:/Users/Lena/OneDrive - University of Otago/Picocyanobacteria/Hydrosphere/Lakes380 Cawthron 2021/Analysis/18S Euks/FINAL_rare20000_18sFIXED_SampleDataFIXED.csv", header = TRUE)
sam_data1 = as.data.frame(sam_data1)
row.names(sam_data1) = sam_data1$"SampleID" #only now possible to replace the rownnames with the first column (OTU ID)
sam_data1$"V1"=NULL #the first column "OTU-ID can then be deleted"
ps_map1 = sample_data(sam_data1)

######### --------------- CREATE THE PHYLOSEQ FILE ------------------- #########

FINAL_rare20000_18sFIXEDAges <- merge_phyloseq(otu_table(FINAL_rare20000_18sFIXED), tax_table(FINAL_rare20000_18sFIXED), ps_map1) 
write.csv(tax_table(FINAL_rare20000_18sFIXEDAges), file="FINAL_rare20000_18sFIXEDAges_tax.csv")

NoTroutNoPlants = subset_taxa(FINAL_rare20000_18sFIXEDAges, Class != "Embryophyceae")
NoTroutNoPlants = subset_taxa(NoTroutNoPlants, Phylum != "Chordata")
NoTroutNoPlants = subset_taxa(NoTroutNoPlants, Family != "Megascolecidae")
ntaxa(NoTroutNoPlants)
# 3123


#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~
#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~
#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~

# FINAL 18S EUKS
saveRDS(NoTroutNoPlants, file = "Euks18SFinal.rds")


### FINAL 18S PREY
Invert = subset_taxa(NoTroutNoPlants, Phylum == "Arthropoda" | Phylum == "Rotifera" | Phylum == "Mollusca")
any(taxa_sums(Invert) == 0)
any(sample_sums(Invert) == 0)
# Invert = prune_taxa(taxa_sums(Invert) > 0, Invert)
# Invert = prune_samples(sample_sums(Invert) > 0, Invert)
ntaxa(Invert)
# 239

saveRDS(Invert, file = "Invert18SFinal.rds")


### FINAL 18S Fish NO TROUT
FishNOTROUT = subset_taxa(FINAL_rare20000_18sFIXEDAges, Class == "Actinopteri")
FishNOTROUT = subset_taxa(FishNOTROUT, Genus != "Salmo")
any(taxa_sums(FishNOTROUT) == 0)
any(sample_sums(FishNOTROUT) == 0)
# FishNOTROUT = prune_taxa(taxa_sums(FishNOTROUT) > 0, FishNOTROUT)
# FishNOTROUT = prune_samples(sample_sums(FishNOTROUT) > 0, FishNOTROUT)
ntaxa(FishNOTROUT)
# 8

saveRDS(FishNOTROUT, file = "FishNOTROUT18SFinal.rds")


### FINAL 18S Fish TROUT

Fish = subset_taxa(FINAL_rare20000_18sFIXEDAges, Class == "Actinopteri")
any(taxa_sums(Fish) == 0)
any(sample_sums(Fish) == 0)
# Fish = prune_taxa(taxa_sums(Fish) > 0, Fish)
# Fish = prune_samples(sample_sums(Fish) > 0, Fish)
ntaxa(Fish)
# 9

saveRDS(Fish, file = "FishwithSalmo18SFinal.rds")


### FINAL 18S Algae

Phyto = subset_taxa(NoTroutNoPlants, Phylum == "Chlorophyta" | Phylum == "Rhodophyta" | Phylum == "Ochrophyta" | Phylum == "Dinoflagellata" | Phylum == "Cryptophyta")
any(taxa_sums(Phyto) == 0)
any(sample_sums(Phyto) == 0)
# Phyto = prune_taxa(taxa_sums(Phyto) > 0, Phyto)
# Phyto = prune_samples(sample_sums(Phyto) > 0, Phyto)
ntaxa(Phyto)
# 40
# Most samples less than 100 reads
saveRDS(Phyto, file = "Algae18SFinal.rds")



##############################################################################################################################
######## 12S rRNA #################################################################################################
############################################################################################################


###shift to HPC
setwd("..Chapter_2_Plate_4_sequences/12S")

path <- "..Chapter_2_Plate_4_sequences/12S"

list.files(path)

###new code

fnFs <- sort(list.files(path, pattern = "R1_001.fastq.gz", full.names = TRUE))
fnRs <- sort(list.files(path, pattern = "R2_001.fastq.gz", full.names = TRUE))

FWD <- "TTAGATACCCCACTATGC"

REV <- "TAGAACAGGCTCCTCTAG"

allOrients <- function(primer) {
  # Create all orientations of the input sequence
  require(Biostrings)
  dna <- DNAString(primer)  # The Biostrings works w/ DNAString objects rather than character vectors
  orients <- c(Forward = dna, Complement = complement(dna), Reverse = reverse(dna), 
               RevComp = reverseComplement(dna))
  return(sapply(orients, toString))  # Convert back to character vector
}
FWD.orients <- allOrients(FWD)
REV.orients <- allOrients(REV)
FWD.orients

primerHits <- function(primer, fn) {
  # Counts number of reads in which the primer is found
  nhits <- vcountPattern(primer, sread(readFastq(fn)), fixed = FALSE)
  return(sum(nhits > 0))
}
rbind(FWD.ForwardReads = sapply(FWD.orients, primerHits, fn = fnFs[[10]]), 
      REV.ReverseReads = sapply(REV.orients, primerHits, fn = fnRs[[10]]))

###cutadapt

cutadapt <- "../miniconda3/bin/cutadapt"

system2(cutadapt, args = "--version") 


###create path

path.cut <- file.path(path, "cutadapt")
if(!dir.exists(path.cut)) dir.create(path.cut)
fnFs.cut <- file.path(path.cut, basename(fnFs))
fnRs.cut <- file.path(path.cut, basename(fnRs))

# Trim FWD off of R1 (forward reads) - 
R1.flags <- paste0("-g", " ^", FWD) 
# Trim REV off of R2 (reverse reads)
R2.flags <- paste0("-G", " ^", REV) 
# Run Cutadapt
for(i in seq_along(fnFs)) {
  system2(cutadapt, args = c("-e 0.08 --discard-untrimmed", R1.flags, R2.flags,
                             "-o", fnFs.cut[i], "-p", fnRs.cut[i], # output files
                             fnFs[i], fnRs[i])) # input files
}

rbind(FWD.ForwardReads = sapply(FWD.orients, primerHits, fn = fnFs.cut[[1]]), 
      REV.ReverseReads = sapply(REV.orients, primerHits, fn = fnRs.cut[[1]]))

# Forward and reverse fastq filenames have the format:
cutFs <- sort(list.files(path.cut, pattern = "R1_001.fastq.gz", full.names = TRUE))
cutRs <- sort(list.files(path.cut, pattern = "R2_001.fastq.gz", full.names = TRUE))

if(length(cutFs) == length(cutRs)) print("Forward and reverse files match. Go forth and explore")
if (length(cutFs) != length(cutRs)) stop("Forward and reverse files do not match. Better go back and have a check")

# Extract sample names, assuming filenames have format:
get.sample.name <- function(fname) strsplit(basename(fname), "_")[[1]][1]
sample.names <- unname(sapply(cutFs, get.sample.name))
head(sample.names)

if(length(cutFs) <= 20) {
  fwd_qual_plots <- plotQualityProfile(cutFs) + 
    scale_x_continuous(breaks=seq(0,250,10)) + 
    scale_y_continuous(breaks=seq(0,40,2)) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1))
  rev_qual_plots <- plotQualityProfile(cutRs) + 
    scale_x_continuous(breaks=seq(0,250,10)) + 
    scale_y_continuous(breaks=seq(0,40,2)) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1))
} else {
  rand_samples <- sample(size = 20, 1:length(cutFs)) # grab 20 random samples to plot
  fwd_qual_plots <- plotQualityProfile(cutFs[rand_samples]) + 
    scale_x_continuous(breaks=seq(0,250,10)) + 
    scale_y_continuous(breaks=seq(0,40,2)) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1))
  rev_qual_plots <- plotQualityProfile(cutRs[rand_samples]) + 
    scale_x_continuous(breaks=seq(0,250,10)) + 
    scale_y_continuous(breaks=seq(0,40,2)) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1))
}

fwd_qual_plots
rev_qual_plots

jpeg(file="Quality.Plot.F.jpg",res=300, width=15, height=8, units="in")
fwd_qual_plots
dev.off()

jpeg(file="Quality.Plot.R.jpg",res=300, width=15, height=8, units="in")
rev_qual_plots
dev.off()

####up to here
filtpathF <- file.path(path.cut, "filtered", basename(cutFs))
filtpathR <- file.path(path.cut, "filtered", basename(cutRs))

###for RV have trialled lots of different trunclen and 110,110 is the best

out <- filterAndTrim(cutFs, filtpathF, cutRs, filtpathR,
                     truncLen=c(110,110), maxEE=c(2,2), truncQ=2, maxN=0, rm.phix=TRUE,
                     compress=TRUE, verbose=TRUE, multithread=TRUE)

out2 <- as.data.frame(out)
str(out2)
out2$perc <- (out2$reads.out/out2$reads.in)*100
out2

sample.names <- sapply(strsplit(basename(filtpathF), "_"), `[`, 1) # Assumes filename = samplename_XXX.fastq.gz
sample.namesR <- sapply(strsplit(basename(filtpathR), "_"), `[`, 1) # Assumes filename = samplename_XXX.fastq.gz
if(identical(sample.names, sample.namesR)) {print("Files are still matching.....congratulations")
} else {stop("Forward and reverse files do not match.")}
names(filtpathF) <- sample.names
names(filtpathR) <- sample.namesR

set.seed(100) # set seed to ensure that randomized steps are replicatable

# Learn forward error rates
errF <- learnErrors(filtpathF, nbases=1e8, multithread=TRUE, verbose = TRUE)
## 100285500 total bases in 668570 reads from 24 samples will be used for learning the error rates.

jpeg(file="Error.F.Plot.jpg",res=300, width=15, height=8, units="in")
plotErrors(errF, nominalQ=TRUE)
dev.off()

# Learn reverse error rates
errR <- learnErrors(filtpathR, nbases=1e8, multithread=TRUE, verbose = TRUE)

jpeg(file="Error.R.Plot.jpg",res=300, width=15, height=8, units="in")
plotErrors(errR, nominalQ=TRUE)
dev.off()


derepF <- derepFastq(filtpathF, verbose=TRUE)
derepR <- derepFastq(filtpathR, verbose=TRUE)

dadaF.pseudo <- dada(derepF, err=errF, multithread=TRUE, pool="pseudo")
dadaR.pseudo <- dada(derepR, err=errR, multithread=TRUE, pool="pseudo")

mergers <- mergePairs(dadaF.pseudo, derepF, dadaR.pseudo, derepR, maxMismatch = 1, minOverlap = 10, 
                      verbose=TRUE)

head(mergers)

seqtab <- makeSequenceTable(mergers)

split.dir <- sapply(strsplit(basename(path), "-"), `[`, 1:2)

split.dir.name  <- paste(split.dir[1], split.dir[2], sep=".")

saveRDS(seqtab, "~/../Chapter_2_Plate_4_sequences/12S/seqtab2.rds")

getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(dadaF.pseudo, getN), sapply(dadaR.pseudo, getN), sapply(mergers, getN))
colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged")
rownames(track) <- sample.names

write.csv(track,"~/../Chapter_2_Plate_4_sequences/12S/track.csv")

track

###Look at the length distribution of the sequences

trimtable <- as.data.frame(table(nchar(getSequences(seqtab))))
colnames(trimtable) <- c("Length.bp", "Frequency")
trimtable$Frequency <- as.numeric(trimtable$Frequency)
str(trimtable)

g <- ggplot(trimtable, aes(x = Length.bp, y = Frequency)) +
  geom_bar(stat="identity")

jpeg(file="Chim.Dist.Plot.R.jpg",res=300, width=15, height=8, units="in")
g

dev.off()

##checking for chimeras (sequences outside the expected size range) Look at table to figure out spread (min,max for amplicon size)

seqtab2 <- seqtab[,nchar(colnames(seqtab)) %in% seq(110,130)]

seqtab.nochim <- removeBimeraDenovo(seqtab2, multithread=TRUE, verbose=TRUE)

###Produce a table which shows the number or reads at each stage

getN <- function(x) sum(getUniques(x))
track2 <- cbind(out, sapply(dadaF.pseudo, getN), sapply(dadaR.pseudo, getN), sapply(mergers, getN), rowSums(seqtab.nochim))
colnames(track2) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
rownames(track2) <- sample.names


trimtable2 <- as.data.frame(table(nchar(getSequences(seqtab.nochim))))
colnames(trimtable2) <- c("Length.bp", "Frequency")
trimtable2$Frequency <- as.numeric(trimtable2$Frequency)
str(trimtable2)

library(ggplot2)
g <- ggplot(trimtable2, aes(x = Length.bp, y = Frequency)) +
  geom_bar(stat="identity")


jpeg(file="No.Chim.Dist.Plot.R.jpg",res=300, width=15, height=8, units="in")
g
dev.off()


###combining Dada2 results

Tutorial.seqtab1 <- readRDS("~/../Chapter_2_Plate_2_sequences/seqtab1.rds")
Tutorial.seqtab2 <- readRDS("~/../Chapter_2_Plate_4_sequences/12S/seqtab2.rds")

Tutorial.dada.all <- mergeSequenceTables(Tutorial.seqtab1, Tutorial.seqtab2)

###Look at the length distribution of the sequences

trimtable <- as.data.frame(table(nchar(getSequences(Tutorial.dada.all))))
colnames(trimtable) <- c("Length.bp", "Frequency")
trimtable$Frequency <- as.numeric(trimtable$Frequency)
str(trimtable)

g <- ggplot(trimtable, aes(x = Length.bp, y = Frequency)) +
  geom_bar(stat="identity")

jpeg(file="Chim.Dist.Plot.R2.jpg",res=300, width=15, height=8, units="in")
g

dev.off()


##checking for chimeras (sequences outside the expected size range) Look at table to figure out spread (min,max for amplicon size)

seqtab2 <- Tutorial.dada.all[,nchar(colnames(seqtab)) %in% seq(110,130)]

seqtab.nochim <- removeBimeraDenovo(seqtab2, multithread=TRUE, verbose=TRUE)

Tutorial.dada.all.nochimera.forsav <- as.data.frame(seqtab.nochim)
write.csv(Tutorial.dada.all.nochimera.forsav, "~/../Chapter_2_Plate_4_sequences/12S/seqtab.nochim.csv", row.names = T)

####assign taxonomy #####

tutorial.nochimera.tax <- assignTaxonomy(seqtab.nochim, "~/../ReferenceDB/fish_eDNA_12S_combined.fasta", multithread=TRUE)

##Save cleaned sequences

Tutorial.all.nochimera.tax.forsav <- as.data.frame(tutorial.nochimera.tax)
write.csv(Tutorial.all.nochimera.tax.forsav, "~/../Chapter_2_Plate_4_sequences/12S/tutorial.nochimera.tax.csv", row.names = T)

####################################################################################################
#####phyloseq#####
####################################################################################################

# Have added the files from here

####Import sequence table into phyloseq - using edited blast plus DaDa2 taxonomy

##asv table

ASV16s <- fread("all12s.nochim_ASVs_counts.tsv", header = T)
head(ASV16s)
ASV16s<- as.data.frame(ASV16s)
rownames(ASV16s) <- ASV16s$V1
head(ASV16s)
ASV16s <- subset(ASV16s, select = -c(`V1`) )
head(ASV16s)

###tax

tax_merge2 <- fread("Blast12S_merge_Dada2_tax_edited2.csv", header = T)
names(tax_merge2)

tax_merge_only_tax <- subset(tax_merge2, select = c("Row.names",
                                                    "Kingdom",
                                                    "Phylum",
                                                    "Class",
                                                    "Order",
                                                    "Family",
                                                    "Genus",
                                                    "Species", 
                                                    "Species.Com.Name"))
tax_merge_only_tax <- as.data.frame(tax_merge_only_tax)

head(tax_merge_only_tax)
names(tax_merge_only_tax)


###get back into phyloseq
ps2_otu_table = otu_table(ASV16s, taxa_are_rows = T)
head(ps2_otu_table)

row.names(tax_merge_only_tax) = tax_merge_only_tax$Row.names
tax_merge_only_tax <- subset(tax_merge_only_tax, select = -c(1) )

ps_tax = tax_table(tax_merge_only_tax)
head(ps_tax)
#adding ID's back
rownames(ps_tax)=rownames(tax_merge_only_tax)
colnames(ps_tax)=colnames(tax_merge_only_tax)
head(ps_tax)


###get sample data

sample.data <-fread("12S_Sample_Data.csv")

sample.data$Lake <- as.factor(sample.data$Lake)
levels(sample.data$Lake)

sample.data$Site <- as.factor(sample.data$Site)
levels(sample.data$Site)

str(sample.data)
sample.data$Age <- as.factor(sample.data$Age)

sample.data$Trout <- factor(sample.data$Trout, levels = c("Before Trout", "After Trout"))

sample.data$Site <- factor(sample.data$Site, levels = c("Depocentre",
                                                        "SC Depocentre", 
                                                        "Littoral",
                                                        "Blank",
                                                        "Water" ))


sample.data$Depth.F <- as.factor(sample.data$Depth.2)
levels(sample.data$Depth.F)


sample.data$Depth.F <- factor(sample.data$Depth.F, levels = c("2.5",
                                                              "3.5", 
                                                              "5.5",
                                                              "5.5 rpt",
                                                              "7",
                                                              "7.5",
                                                              "9",
                                                              "9 rpt", 
                                                              "9.5",
                                                              "10.5",
                                                              "11",
                                                              "11.5",
                                                              "12",
                                                              "12.5",
                                                              "13",
                                                              "13.5",
                                                              "13.5 rpt",
                                                              "14",
                                                              "14.5",
                                                              "15",
                                                              "15.5",
                                                              "16.5",
                                                              "17", 
                                                              "17.5",
                                                              "18.5",
                                                              "19",
                                                              "19.5",
                                                              "20.5",
                                                              "21",
                                                              "22",
                                                              "22.5",
                                                              "23", 
                                                              "24",
                                                              "24 rpt",
                                                              "25", 
                                                              "26.5",
                                                              "27",
                                                              "28.5",
                                                              "29.5", 
                                                              "30",
                                                              "30.5",
                                                              "31",
                                                              "31.5",
                                                              "32",
                                                              "32.5",
                                                              "34.5",
                                                              "35",
                                                              "36",
                                                              "36 rpt",
                                                              "37",
                                                              "38",
                                                              "39",
                                                              "40",
                                                              "41",
                                                              "42",
                                                              "44",
                                                              "44.5",
                                                              "45",
                                                              "45.5",
                                                              "46",
                                                              "46.5",
                                                              "47.5",
                                                              "47.5 rpt",
                                                              "48.5",
                                                              "49",
                                                              "50",
                                                              "54",
                                                              "55",
                                                              "60",
                                                              "64",
                                                              "65",
                                                              "70",
                                                              "70.5",
                                                              "75",
                                                              "78",
                                                              "80",
                                                              "85.5",
                                                              "88",
                                                              "90",
                                                              "91",
                                                              "94",   
                                                              "95",
                                                              "100",
                                                              "104",
                                                              "104.5",
                                                              "106",
                                                              "108",
                                                              "110",
                                                              "112",
                                                              "114.5",
                                                              "120",
                                                              "130",
                                                              "140",
                                                              "155" ,
                                                              "163",
                                                              "168.5",
                                                              "172",
                                                              "173.5",
                                                              "183",
                                                              "212",
                                                              "219.5",
                                                              "228",
                                                              "237.5",
                                                              "250",
                                                              "411",
                                                              "449",
                                                              "Blank",
                                                              "Blank 1", 
                                                              "Blank 2", 
                                                              "Water"))


ps_map = sample_data(sample.data)
head(ps_map)
rownames(ps_map)=sample.data$'SampleID'
head(ps_map)

##new phyloseq

phyloseq <- phyloseq(ps2_otu_table, ps_map, ps_tax)

phyloseq


####Write ASV tax table

phylo_paringa = subset_samples(phyloseq, Lake == "Paringa")
phylo_paringa1 = prune_taxa(taxa_sums(phylo_paringa) > 0, phylo_paringa)

ASV <- as.data.frame(otu_table(phylo_paringa1))
ASV <- t(ASV)
head(ASV)

TAX <- as.data.frame(tax_table(phylo_paringa1))
head(TAX)

TAX.ASV <- merge(ASV,TAX, by="row.names")

saveRDS(phylo_paringa1, "Chpt2RVphyloseq.rds")

phylo_paringa1 = readRDS(file = "Chpt2RVphyloseq.rds")

phylo_paringa1.depths <- subset_samples(phylo_paringa1, Depth.F %in% c("2.5",
                                                                       "3.5",
                                                                       "5.5 rpt",
                                                                       "7",
                                                                       "7.5",
                                                                       "9 rpt", 
                                                                       "9.5",
                                                                       "10.5",
                                                                       "11",
                                                                       "11.5",
                                                                       "12",
                                                                       "12.5",
                                                                       "13",
                                                                       "13.5 rpt",
                                                                       "14",
                                                                       "14.5",
                                                                       "15",
                                                                       "15.5",
                                                                       "16.5",
                                                                       "17", 
                                                                       "17.5",
                                                                       "18.5",
                                                                       "19",
                                                                       "19.5",
                                                                       "20.5",
                                                                       "21",
                                                                       "22",
                                                                       "22.5",
                                                                       "23", 
                                                                       "24 rpt",
                                                                       "25", 
                                                                       "26.5",
                                                                       "27",
                                                                       "28.5",
                                                                       "29.5", 
                                                                       "30",
                                                                       "30.5",
                                                                       "31",
                                                                       "31.5",
                                                                       "32",
                                                                       "32.5",
                                                                       "34.5",
                                                                       "35",
                                                                       "36 rpt",
                                                                       "37",
                                                                       "38",
                                                                       "39",
                                                                       "40",
                                                                       "41",
                                                                       "42",
                                                                       "44",
                                                                       "44.5",
                                                                       "45",
                                                                       "45.5",
                                                                       "46",
                                                                       "46.5",
                                                                       "47.5 rpt",
                                                                       "48.5",
                                                                       "49",
                                                                       "50",
                                                                       "54",
                                                                       "55",
                                                                       "60",
                                                                       "64",
                                                                       "65",
                                                                       "70",
                                                                       "70.5",
                                                                       "75",
                                                                       "78",
                                                                       "80",
                                                                       "85.5",
                                                                       "88",
                                                                       "90",
                                                                       "91",
                                                                       "94",   
                                                                       "95",
                                                                       "100",
                                                                       "Blank",
                                                                       "Blank 1", 
                                                                       "Blank 2", 
                                                                       "Water"))
####

####subtraction of contamination

Controls = subset_samples(phylo_paringa1.depths, Site %in% c("Blank"))
Controls = filter_taxa(Controls, function(x) sum(x) > 0, TRUE)
sample_sums(Controls)

plot_bar(Controls, fill = "Genus")


###removing ASV contamination - ASV total removal

Controls.ASVs <- taxa_names(Controls)

allTaxa = taxa_names(phylo_paringa1.depths)
allTaxa.filtered <- allTaxa[!(allTaxa %in% Controls.ASVs)]
phylo_paringa1.depths.nocont = prune_taxa(allTaxa.filtered, phylo_paringa1.depths)

##compare sample sums with and without ASV removal

org.ss <- as.data.frame(sample_sums(phylo_paringa1.depths))
org.ss$Names <- rownames(org.ss)
new.ss <- as.data.frame(sample_sums(phylo_paringa1.depths.nocont))
new.ss$Names <- rownames(new.ss)
control.rem.sums <- dplyr::left_join(org.ss, new.ss, by="Names")
colnames(control.rem.sums) <- c("Original", "Names", "New")
control.rem.sums.final <- control.rem.sums %>%
  mutate(perc = New/Original*100)
control.rem.sums.final

##subtract AVS instead
###find the max value for each ASV and store as a vector

Extraction_neg = subset_samples(phylo_paringa1.depths, Site == "Blank")

Extraction_neg_max <- apply(as.data.frame(as.matrix(t(otu_table(Extraction_neg)))), 1, max)
Extraction_neg_max_vec <- as.vector(Extraction_neg_max)
Extraction_neg_max_vec

Extraction_neg_sums <- colSums(otu_table(Controls))
Extraction_neg_sums_vec <-as.vector(Extraction_neg_sums)
Extraction_neg_sums_vec

Extraction = as(otu_table(phylo_paringa1.depths), "matrix")
Extractiondf = as.data.frame(Extraction)

Extractiondf[,1:length(Extractiondf)] <- sweep(Extractiondf[,1:length(Extractiondf)],2,Extraction_neg_sums_vec)
Extractiondf <- replace(Extractiondf, Extractiondf < 0, 0)

##make new phyloseq

phylo_paringa1.depths.subtract.ps <- phyloseq(otu_table(Extractiondf, taxa_are_rows=TRUE),
                                              sample_data(phylo_paringa1.depths), tax_table(phylo_paringa1.depths))


##compare phyloseq objects

org.ss <- as.data.frame(sample_sums(phylo_paringa1.depths))
org.ss$Names <- rownames(org.ss)
new.ss <- as.data.frame(sample_sums(phylo_paringa1.depths.subtract.ps))
new.ss$Names <- rownames(new.ss)
control.rem.maxsub <- dplyr::left_join(org.ss, new.ss, by="Names")
colnames(control.rem.maxsub) <- c("Original", "Names", "New")
control.rem.maxsub.final <- control.rem.maxsub %>%
  mutate(perc = New/Original*100)
control.rem.maxsub.final


phylo_paringa1.5 = subset_samples(phylo_paringa1.depths.subtract.ps, Site != "Blank")

# added additional samples from other run
Paringa.all12S.Phyloseq = readRDS("Paringa.all12S.Phyloseq.rds")

#######
newtroutdata <- subset_samples(Paringa.all12S.Phyloseq, Depth.F %in% c("2.5",
                                                                       "3.5",
                                                                       "5.5 rpt",
                                                                       "7",
                                                                       "7.5",
                                                                       "9 rpt", 
                                                                       "9.5",
                                                                       "10.5",
                                                                       "11",
                                                                       "11.5",
                                                                       "12",
                                                                       "12.5",
                                                                       "13",
                                                                       "13.5 rpt",
                                                                       "14",
                                                                       "14.5",
                                                                       "15",
                                                                       "15.5",
                                                                       "16.5",
                                                                       "17", 
                                                                       "17.5",
                                                                       "18.5",
                                                                       "19",
                                                                       "19.5",
                                                                       "20.5",
                                                                       "21",
                                                                       "22",
                                                                       "22.5",
                                                                       "23", 
                                                                       "24 rpt",
                                                                       "25", 
                                                                       "26.5",
                                                                       "27",
                                                                       "28.5",
                                                                       "29.5", 
                                                                       "30",
                                                                       "30.5",
                                                                       "31",
                                                                       "31.5",
                                                                       "32",
                                                                       "32.5",
                                                                       "34.5",
                                                                       "35",
                                                                       "36 rpt",
                                                                       "37",
                                                                       "38",
                                                                       "39",
                                                                       "40",
                                                                       "41",
                                                                       "42",
                                                                       "44",
                                                                       "44.5",
                                                                       "45",
                                                                       "45.5",
                                                                       "46",
                                                                       "46.5",
                                                                       "47.5 rpt",
                                                                       "48.5",
                                                                       "49",
                                                                       "50",
                                                                       "54",
                                                                       "55",
                                                                       "60",
                                                                       "64",
                                                                       "65",
                                                                       "70",
                                                                       "70.5",
                                                                       "75",
                                                                       "78",
                                                                       "80",
                                                                       "85.5",
                                                                       "88",
                                                                       "90",
                                                                       "91 rpt",
                                                                       "94",   
                                                                       "95",
                                                                       "100",
                                                                       "Blank",
                                                                       "Blank 1", 
                                                                       "Blank 2", 
                                                                       "Water"))


##rarefy

sample_sums(newtroutdata)
set.seed(100)
phylo_paringa1.5r5000 = rarefy_even_depth(newtroutdata, sample.size = 5000, 
                                          replace = FALSE, trimOTUs = TRUE, verbose = TRUE)


###select only fish

phylo_paringa1.5r5000FISH = subset_taxa(phylo_paringa1.5r5000, Class %in% "Actinopteri")
phylo_paringa1.5r5000FISH  = prune_taxa(taxa_sums(phylo_paringa1.5r5000FISH) > 0, phylo_paringa1.5r5000FISH)
phylo_paringa1.5r5000FISH  <- prune_samples(sample_sums(phylo_paringa1.5r5000FISH ) > 0,phylo_paringa1.5r5000FISH )
sample_sums(phylo_paringa1.5r5000FISH)
write.csv(otu_table(phylo_paringa1.5r5000FISH), file = "fishreads.csv")

phylo_paringa1.5r5000FISHNEW = subset_samples(phylo_paringa1.5r5000FISH, Age != "2015" & Age != "2017")

saveRDS(phylo_paringa1.5r5000FISHNEW, file = "12SFish.rds")

#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~
#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~
#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~






















### Statistical Analysis for: Insights into the ecological impact of trout introduction in an 
### oligotrophic lake using sedimentary environmental DNA 


# The datasets:

~#~#~#~#~#~ FISH ~#~#~#~#~#~#~#~
  
  Fish = readRDS("FishwithSalmo18SFinal.rds")

Fish12S = readRDS("12SFish.rds")

#~#~#~#~# EUKARYOTES #~#~#~#~#~#

# Final 18S Euks
Euk18S = readRDS("Euks18SFinal.rds")
ntaxa(Euk18S)
# 3123
nsamples(Euk18S)
#37

# Final COI Long Euks
EukCOIL = readRDS("COILongEuksFinal.rds")
ntaxa(EukCOIL)
# 3948
nsamples(EukCOIL)
#28

# Final COI Short Euks
EukCOIS = readRDS("UPDATED25000RARE_COIShortEuksFinal.rds")
ntaxa(EukCOIS)
# 6367 UPDATED: 7521
nsamples(EukCOIS)
#28 UPDATED: 27

#~#~#~#~#~#~#~#~# PREY #~#~#~#~#~#~#~#~##~#~

# Zooplankton 18S
zoo18S = readRDS("Invert18SFinal.rds")
ntaxa(zoo18S)
#239
nsamples(zoo18S)
#37

# Zooplankton Short COI
zooshort = readRDS("UPDATED_25000RARE_COIShortInvertsFinal.rds")
ntaxa(zooshort)
#690 UPDATED: 869
nsamples(zooshort)
#28 UPDATED: 27

#Zooplankton Long COI
zoolong = readRDS("COILongInvertsFinal.rds")
ntaxa(zoolong)
#536
nsamples(zoolong)
#28

#~#~#~#~#~#~#~ ALGAE #~#~#~#~#~#~#~#~#

Algae18S = readRDS("Algae18SFinal.rds")
ntaxa(Algae18S)
#40
nsamples(Algae18S)
#37

# Algae Long COI
AlgaeLong = readRDS("COILongAlgaeFinal.rds")
ntaxa(AlgaeLong)
# 83
nsamples(AlgaeLong)
#28

# Algae Short COI
AlgaeShort = readRDS("UPDATED_25000RARE_COIShortAlgaeFinal.rds")
ntaxa(AlgaeShort)
# 125 UPDATED: 142
nsamples(AlgaeShort)
#28 UPDATED: 27



################### ORDERED BY DEPTH FOR PLOTS #########################

Euk18SOrder <- Euk18S %>% ps_arrange(Depth)
EukCOILOrder <- EukCOIL %>% ps_arrange(Depth)
EukCOISOrder <- EukCOIS %>% ps_arrange(Depth)
zoo18SOrder <- zoo18S %>% ps_arrange(Depth)
zooshortOrder <- zooshort %>% ps_arrange(Depth)
zoolongOrder <- zoolong %>% ps_arrange(Depth)
AlgaeLongOrder <- AlgaeLong %>% ps_arrange(Depth)
AlgaeShortOrder <- AlgaeShort %>% ps_arrange(Depth)
Algae18SOrder <- Algae18S %>% ps_arrange(Depth)


# adding sample interval to COI data to make GAM

#write.csv(sample_data(zooshortOrder), file = "zooShort_sample.csv")
sampzooshort = read.csv(file="zooShort_sample_interval.csv")
row.names(sampzooshort)=sampzooshort$"SampleID"
sampzooshort1 = sample_data(sampzooshort)
zooshortfinal <- merge_phyloseq(otu_table(zooshortOrder), tax_table(zooshortOrder), sampzooshort1)

#write.csv(sample_data(zoolongOrder), file = "zooLong_sample.csv")
sampzoolong = read.csv(file="zooLong_sample_interval.csv")
row.names(sampzoolong)=sampzoolong$"SampleID"
sampzoolong1 = sample_data(sampzoolong)
zoolongfinal <- merge_phyloseq(otu_table(zoolongOrder), tax_table(zoolongOrder), sampzoolong1)

#write.csv(sample_data(EukCOILOrder), file = "EukCOIL_sample.csv")
sampeuklong = read.csv(file="EukCOIL_sample_interval.csv")
row.names(sampeuklong)=sampeuklong$"SampleID"
sampeuklong1 = sample_data(sampeuklong)
euklongfinal <- merge_phyloseq(otu_table(EukCOILOrder), tax_table(EukCOILOrder), sampeuklong1)

#write.csv(sample_data(EukCOISOrder), file = "EukCOIS_sample.csv")
sampeukshort = read.csv(file="EukCOIS_sample_interval.csv")
row.names(sampeukshort)=sampeukshort$"SampleID"
sampeukshort1 = sample_data(sampeukshort)
eukshortfinal <- merge_phyloseq(otu_table(EukCOISOrder), tax_table(EukCOISOrder), sampeukshort1)

#write.csv(sample_data(AlgaeLongOrder), file = "AlgaeLong_sample.csv")
sampalgaelong = read.csv(file="AlgaeLong_sample_interval.csv")
row.names(sampalgaelong)=sampalgaelong$"SampleID"
sampalgaelong1 = sample_data(sampalgaelong)
algaelongfinal <- merge_phyloseq(otu_table(AlgaeLongOrder), tax_table(AlgaeLongOrder), sampalgaelong1)


#write.csv(sample_data(AlgaeShortOrder), file = "AlgaeShort_sample.csv")
sampalgaeshort = read.csv(file="AlgaeShort_sample_interval.csv")
row.names(sampalgaeshort)=sampalgaeshort$"SampleID"
sampalgaeshort1 = sample_data(sampalgaeshort)
algaeshortfinal <- merge_phyloseq(otu_table(AlgaeShortOrder), tax_table(AlgaeShortOrder), sampalgaeshort1)


############### FISH HEATMAPS

Species_glom <-tax_glom(Fish12S, taxrank="Species")
Species_glom = prune_taxa(taxa_sums(Species_glom) > 0, Species_glom)
Species_glom <- prune_samples(sample_sums(Species_glom) > 0,Species_glom)


Species_glom.dataf = psmelt(Species_glom)
Species_glom.dataf$Trout <- factor(Species_glom.dataf$Trout, levels = c("Before Trout", "After Trout"))
Species_glom.dataf$Age <- as.factor(Species_glom.dataf$Age)
levels(Species_glom.dataf$Age)
levels(Species_glom.dataf$Species)

Species_glom.dataf$Species <- factor(Species_glom.dataf$Species, 
                                     levels = c("Oncorhynchus tshawytscha","Gobiomorphus breviceps",
                                                "Anguilla dieffenbachii", "Anguilla australis", 
                                                "Galaxias paucispondylus", "Retropinna retropinna", 
                                                "Galaxias sp.",  "Galaxias fasciatus", 
                                                "Galaxias argenteus", "Salmo trutta"))

Species_glom.dataf$Age <- factor(Species_glom.dataf$Age, levels = c("1731",
                                                                    "1768",
                                                                    "1776",
                                                                    "1810",
                                                                    "1847",
                                                                    "1866",
                                                                    "1870",
                                                                    "1883",
                                                                    "1888",
                                                                    "1891",
                                                                    "1894",
                                                                    "1906",
                                                                    "1912",
                                                                    "1920",
                                                                    "1925",
                                                                    "1929",
                                                                    "1933",
                                                                    "1937",
                                                                    "1946",
                                                                    "1951",
                                                                    "1952",
                                                                    "1954",
                                                                    "1956",
                                                                    "1963",
                                                                    "1964",
                                                                    "1965",
                                                                    "1966",
                                                                    "1972",
                                                                    "1986",
                                                                    "1997"
))


#Species_glom.dataf$logAbundance = sqrt(Species_glom.dataf$Abundance)

ABUN.heatmap <- ggplot(Species_glom.dataf, aes(x = Species, 
                                               y = Age, 
                                               fill = Abundance)) + 
  scale_fill_gradient(low = "white", high = "dodgerblue4", na.value = "white") + 
  geom_tile(colour = "white") + 
  facet_grid(~Trout, scales = "free", space = "free") + coord_flip() + 
  theme(
    panel.border = element_rect(colour = "black", fill = NA),
    panel.background = element_blank(),
    strip.background = element_rect(colour = "black", fill = NA),
    strip.text = element_text(face = "bold"),
    panel.grid.minor = element_blank(),
    #legend.title = element_text(face = "italic"),
    axis.title.x= element_text(face = "bold"),
    axis.title.y = element_text(face = "bold"),
    axis.text.y = element_text(face = "italic"),
    axis.text.x = element_text(size = 10, angle = 45, hjust = 1, vjust = 1)) 

ABUN.heatmap
#1000x300

################################################################################################ 
########18S FIsh
################################################################################################ 


Species_glomfish <-tax_glom(Fish, taxrank="Genus")
Species_glomfish = prune_taxa(taxa_sums(Species_glomfish) > 0, Species_glomfish)
Species_glomfish <- prune_samples(sample_sums(Species_glomfish) > 0,Species_glomfish)

Species_glom.datafish = psmelt(Fish)

Species_glom.datafish$FishBefore <- factor(Species_glom.datafish$FishBefore, levels = c("Before", "After"))
Species_glom.datafish$Age <- as.factor(Species_glom.datafish$Age)
levels(Species_glom.datafish$Age)

Species_glom.datafish$Age <- factor(Species_glom.datafish$Age, levels = c("1768",
                                                                          "1776",
                                                                          "1781",
                                                                          "1793",
                                                                          "1810",
                                                                          "1830",
                                                                          "1847",
                                                                          "1856",
                                                                          "1866",
                                                                          "1870",
                                                                          "1883",
                                                                          "1886",
                                                                          "1888",
                                                                          "1891",
                                                                          "1894",
                                                                          "1906",
                                                                          "1909",
                                                                          "1912",
                                                                          "1914",
                                                                          "1920",
                                                                          "1925",
                                                                          "1929",
                                                                          "1933",
                                                                          "1937",
                                                                          "1946",
                                                                          "1950",
                                                                          "1951",
                                                                          "1952",
                                                                          "1956",
                                                                          "1963",
                                                                          "1964",
                                                                          "1965",
                                                                          "1966",
                                                                          "1967",
                                                                          "1972",
                                                                          "1981",
                                                                          "1986"))

Species_glom.datafish$logAbundance = log(Species_glom.datafish$Abundance)
#max(Species_glom.datafish$logAbundance)

ABUN.heatmap <- ggplot(Species_glom.datafish, aes(x = Genus, 
                                                  y = Age, 
                                                  fill = logAbundance)) + 
  scale_fill_gradient(low = "white", high = "darkslategrey", na.value = "white") + 
  geom_tile(colour = "white") + 
  facet_grid(~FishBefore, scales = "free", space = "free") + coord_flip() + 
  scale_x_discrete(limits = rev(levels(Species_glom.datafish$Genus))) +
  theme(
    panel.border = element_rect(colour = "black", fill = NA),
    panel.background = element_blank(),
    strip.background = element_rect(colour = "black", fill = NA),
    strip.text = element_text(face = "bold"),
    panel.grid.minor = element_blank(),
    #legend.title = element_text(face = "italic"),
    axis.title.x= element_text(face = "bold"),
    axis.title.y = element_text(face = "bold"),
    axis.text.y = element_text(face = "italic"),
    axis.text.x = element_text(size = 10, angle = 45, hjust = 1, vjust = 1)) 

ABUN.heatmap
#1000x500

"#4682b4"


################ CONNIS PLOTS

################## INPUT THE ORDERED SUBSET GROUPS ########################

###########################################################################################
################## Run from here individually for each data subset! ########################
###########################################################################################

downcore.ps = prune_samples(sample_sums(AlgaeLongOrder)>0, AlgaeLongOrder) %>%
  filter_taxa(function(x) sum(x) > 0, TRUE)

downcore.ps.perc.df <- data.frame(otu_table(downcore.ps))

## TAXONOMY to rename OTU labels
downcore.ps.tax.df <- data.frame(tax_table(downcore.ps))
rownames(downcore.ps.tax.df) <- downcore.ps.tax.df$Genus

#Transpose OTU table
downcore.ps.perc.df<-as.data.frame(t(downcore.ps.perc.df))
# make a column of OTU table the taxonomy
downcore.ps.perc.df$otulabels <- downcore.ps.tax.df$Genus
rownames(downcore.ps.perc.df) <- downcore.ps.perc.df$otulabels

downcore.ps.perc.df$`otulabels`=NULL
# transpose back
downcore.ps.perc.df<-as.data.frame(t(downcore.ps.perc.df))

# metadata
downcore.ps.perc.metadata <- data.frame(sample_data(downcore.ps))

rownames(downcore.ps.perc.df) <- downcore.ps.perc.metadata$Age

downcore.ps.perc.df_hel <-  decostand(downcore.ps.perc.df, "hellinger")


###########################################################################################
######################## CONNIS STRAT PLOTS ##################################################
###########################################################################################

# alltogether plot

strat.per <- downcore.ps.perc.df
# names(strat.per) <- paste("gene",
#                           seq(1, (ncol(strat.per))),sep = "")
#names(strat.per)
rowSums(strat.per) #not percentages

strat.clust <- tran(strat.per, method = "proportion") %>% 
  vegdist(method = "bray", binary = FALSE, diag = FALSE, upper = FALSE,
          na.rm = FALSE) %>% 
  chclust(method = "coniss")

bstick(strat.clust, 20)
strat.clust.sig <- 3 #2 or 3

strat.per.red <- tran(strat.per, method = "percent") %>% 
  chooseTaxa(n.occ = 2, max.abun = 4)


strat.per.ex <- names(strat.per.red) %in% names(which(apply(strat.per.red, 2,
                                                            max) < 5))

x <- strat.plot(strat.per.red, yvar = downcore.ps.perc.metadata$Depth,
                y.tks = seq(round_any(min(downcore.ps.perc.metadata$Depth),
                                      5,
                                      f = floor),
                            round_any(max(downcore.ps.perc.metadata$Depth),
                                      5,
                                      f = floor),
                            5),
                x.names = names(strat.per.red),
                mgp =  c(3, 0.9, 0.3), yBottom = 0.1,
                cex.xlabel = 0.75, y.rev = TRUE, xLeft = 0.16, 
                plot.poly = TRUE, plot.bar = TRUE, col.bar = "black",
                col.poly = "darkslategray3", col.poly.line = "black",
                scale.percent = TRUE, xSpace = 0.01, x.pc.lab = TRUE,
                x.pc.omit0 = TRUE, las = 2, exag = strat.per.ex,
                col.exag = "auto", exag.alpha = 0.50, exag.mult = 10,
                srt.xlabel = 45, x.pc.inc = 10, clust = strat.clust)
addClustZone(x, strat.clust, nZone = strat.clust.sig,
             lwd = 1.5, lty = 2, col = "grey25")
secondary_scale_mod(x, yvar = downcore.ps.perc.metadata$Depth,
                    yvar2 = as.numeric(rownames(strat.per.red)),
                    
                    xLeft = 0.13, n = 20, ylabel2 = "Age (CE)")

# 1200x700



################## INPUT HERE THE ORDERED SUBSET GROUPS ########################

downcore.ps = prune_samples(sample_sums(Euk18SOrder)>0, Euk18SOrder) %>%
  filter_taxa(function(x) sum(x) > 0, TRUE)

downcore.ps.perc <- transform_sample_counts(downcore.ps, function(OTU) OTU/sum(OTU)*100)

downcore.ps.perc.df <- data.frame(otu_table(downcore.ps.perc))

## TAXONOMY to rename OTU labels
downcore.ps.tax.df <- data.frame(tax_table(downcore.ps))
rownames(downcore.ps.tax.df) <- downcore.ps.tax.df$Genus

#Transpose OTU table
downcore.ps.perc.df<-as.data.frame(t(downcore.ps.perc.df))
# make a column of OTU table the taxonomy
downcore.ps.perc.df$otulabels <- downcore.ps.tax.df$Genus
rownames(downcore.ps.perc.df) <- downcore.ps.perc.df$otulabels

downcore.ps.perc.df$`otulabels`=NULL
# transpose back
downcore.ps.perc.df<-as.data.frame(t(downcore.ps.perc.df))

# metadata
downcore.ps.perc.metadata <- data.frame(sample_data(downcore.ps))

rownames(downcore.ps.perc.df) <- downcore.ps.perc.metadata$Age


strat.clust <- downcore.ps.perc.df %>% 
  vegdist(method = "bray", binary = FALSE, diag = FALSE, upper = FALSE,
          na.rm = FALSE) %>% 
  chclust(method = "coniss")

#bstick(strat.clust, 20)
downcore.ps.perc.metadata$Age <- as.numeric(as.character(downcore.ps.perc.metadata$Age))
ct <- cutree(strat.clust, 2)
zn <- which(diff(ct) > 0)
zone <- (downcore.ps.perc.metadata$Age[zn] + downcore.ps.perc.metadata$Age[zn + 1])/2
zn.df <- data.frame(zone)

dendro <- as.dendrogram(strat.clust)
ddata <- dendro_data(dendro, type="rectangle")

depth <- downcore.ps.perc.metadata$Age


ddata$segments->yo
yo$xx=NA
yo$xxend=NA

for(i in 1:length(depth)){
  yo$xx[which(yo$x == i)]=depth[i]
  yo$xxend[which(yo$xend == i)]=depth[i]
  da=depth[i+1]-depth[i]
  yo$xx[yo$x>i & yo$x<i+1]<-depth[i]+(yo$x[yo$x>i & yo$x<i+1]-i)*da
  yo$xxend[yo$xend>i & yo$xend<i+1]<-depth[i]+(yo$xend[yo$xend>i & yo$xend<i+1]-i)*da
} # Please dont ask....
ddata$segments$x<-yo$xx
ddata$segments$xend<-yo$xxend

theme_new <- theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), # remove grids
                   panel.background = element_blank(), axis.line = element_line(colour = "black"),
                   strip.text.x = element_text(size=10, angle=45, hjust=0, vjust=0.05), 
                   # strip.text.x = element_blank(), # Taxa names
                   strip.background = element_blank(),
                   strip.text.y = element_text(angle = 0),
                   legend.position="none",panel.border = element_blank(),
                   axis.text.x=element_text(angle=45,hjust=1, size = 11),
                   axis.text.y = element_text(size=11),
                   axis.title = element_text(size=14),
                   plot.title = element_text(size=15, vjust = 40)) +
  theme(panel.background = element_rect(fill="transparent"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_rect(fill="transparent"),
        plot.background = element_rect(fill="transparent", colour=NA),
        legend.background = element_rect(fill="transparent"),
        legend.key = element_rect(fill = NA, colour = NA))



dendro.18S =ggplot(segment(ddata)) +
  #geom_vline(data=zn.df, mapping = aes(xintercept=1906), linetype="dashed", color = "red") + 
  #geom_rect(data= rects, aes(xmin= xstart, xmax=xend, ymin = ystart, ymax=yend), fill = "red", alpha=0.1) +
  geom_segment(aes(x=x, y=y, xend=xend, yend=yend)) +
  coord_flip() +
  ylab("Distance")+
  #scale_x_reverse()+
  theme_new + theme(axis.title.y=element_blank(),
                    axis.text.y=element_blank(),
                    axis.ticks.y=element_blank(),
                    axis.line.y = element_blank())   +
  scale_x_continuous(expand=c(0,1), limits = c(1760, 1990)) + 
  scale_y_continuous(expand=c(0,0), limits = c(0,12)) +
  ggtitle("")


strat.per.red <- downcore.ps.perc.df %>% 
  chooseTaxa(n.occ = 20)


strat.per.red1 <- strat.per.red %>% select("Chytridiomycetes", "Leidyana1_1", "Unclassified_Kingdom_Eukaryota_2", "Unclassified_Kingdom_Eukaryota_3")

#colnames(strat.per.red1) <- c("Chytridiomycetes         ", "Leidyana1_1", "Unclass_Euk_2", "Unclass_Euk_3")
strat.per.red.depth <- rownames_to_column(strat.per.red1) 

colnames(strat.per.red.depth)[1] <- "Age"

strat.per.red.depth <- strat.per.red.depth %>% arrange(Age) 

strat.per.red.depth.melt <- reshape2::melt(strat.per.red.depth, id.vars=c("Age"))

strat.per.red.depth.melt$Age <- as.numeric(strat.per.red.depth.melt$Age)

strat.per.red.depth.melt <- strat.per.red.depth.melt %>% add_row(Age = 1500, variable = "Leidyana1_1", value= 8)


strat.per.red.depth.melt.p1 <-  strat.per.red.depth.melt %>% filter(variable == "Chytridiomycetes")
strat.per.red.depth.melt.p1$variable <- gsub("Chytridiomycetes", "Chytridiomycetes sp.", strat.per.red.depth.melt.p1$variable)

strat.per.red.depth.melt.p2 <-  strat.per.red.depth.melt %>% filter(variable == "Leidyana1_1")
strat.per.red.depth.melt.p2$variable <- gsub("Leidyana1_1", "Leidyana sp.", strat.per.red.depth.melt.p2$variable)

strat.per.red.depth.melt.p3 <-  strat.per.red.depth.melt %>% filter(variable == "Unclassified_Kingdom_Eukaryota_2")
strat.per.red.depth.melt.p3$variable <- gsub("Unclassified_Kingdom_Eukaryota_2", "Unclassified Eukaryote", strat.per.red.depth.melt.p3$variable)

strat.per.red.depth.melt.p4 <-  strat.per.red.depth.melt %>% filter(variable == "Unclassified_Kingdom_Eukaryota_3")
strat.per.red.depth.melt.p4$variable <- gsub("Unclassified_Kingdom_Eukaryota_3", "Unclassified Eukaryote", strat.per.red.depth.melt.p4$variable)




rects <- data.frame(ystart = c(0), yend = c(Inf), xstart = c(1906), xend = c(1870))

myGrob <- grobTree(rectGrob(gp=gpar(fill="red", alpha=0.1)))

myGrob2 <- grobTree(segmentsGrob(x0=0, x1=1, y0=0, y1=0, default.units="npc", gp=gpar(lty="dashed", col="red", lwd=2)))

myGrob3 <- grobTree(segmentsGrob(x0=0, x1=1, y0=0, y1=0, default.units="npc", gp=gpar(lty="dashed", col="black", lwd=3)))


myGrob4 <- grobTree(textGrob("Chytridiomycetes sp.", rot=60))

plot.18S.p1 <- ggplot()+
  #geom_rect(data= rects, aes(xmin= xstart, xmax=xend, ymin = ystart, ymax=yend), fill = "red", alpha=0.1) +
  geom_line(data=strat.per.red.depth.melt.p1, aes(Age, value), colour="#21b6a8") +
  geom_area(data=strat.per.red.depth.melt.p1, aes(Age, value), fill="#21b6a8")+ 
  geom_bar(data=strat.per.red.depth.melt.p1, aes(Age, value), width=0.1, stat="identity")+
  #facet_grid(~variable, scales = "fixed", labeller = label_wrap_gen(multi_line = TRUE)) +
  #geom_vline(data=zn.df, mapping = aes(xintercept=zone), linetype="dashed", color = "black") + 
  # geom_vline(data=zn.df, mapping = aes(xintercept=1906), linetype="dashed", color = "red") + 
  xlab("Age (CE)")+ylab("") +
  coord_flip() +
  theme_new  +
  scale_x_continuous(expand=c(0,1), limits = c(1760, 1990), breaks = c(1775, 1800, 1825, 1850, 1875, 1900, 1925, 1950, 1975, 2000)) + 
  scale_y_continuous(expand=c(0,0)) +
  theme(plot.margin = unit(c(10,0,0,0), "lines")) +
  theme(axis.text.y = element_text(color = c(NA, "black", NA, "black", NA,  "black", NA,  "black", NA, "black"))) +
  ggtitle("a) Euk18S") +
  annotation_custom(myGrob, xmin=1906, xmax=1870, ymin=0, ymax=700) +
  annotation_custom(myGrob4, xmin=2000, xmax=2040, ymin=3, ymax=20) 


plot.18S.p1.g <- ggplotGrob(plot.18S.p1)
plot.18S.p1.g$layout$clip[plot.18S.p1.g$layout$name=="panel"] <- "off"

myGrob4 <- grobTree(textGrob("Leidyana sp.", rot=60))


plot.18S.p2 <- ggplot() +
  #geom_rect(data= rects, aes(xmin= xstart, xmax=xend, ymin = ystart, ymax=yend), fill = "red", alpha=0.1) +
  geom_line(data=strat.per.red.depth.melt.p2, aes(Age, value), colour="#21b6a8") +
  geom_area(data=strat.per.red.depth.melt.p2, aes(Age, value), fill="#21b6a8") + 
  geom_bar(data=strat.per.red.depth.melt.p2, aes(Age, value), width=0.1, stat="identity")+
  #facet_grid(~variable, scales = "fixed", labeller = label_wrap_gen(multi_line = TRUE)) +
  #geom_vline(data=zn.df, mapping = aes(xintercept=zone), linetype="dashed", color = "black") + 
  # geom_vline(data=zn.df, mapping = aes(xintercept=1906), linetype="dashed", color = "red") + 
  xlab("Age")+ylab("Relative Abundance (%)") +
  coord_flip() +
  theme_new  +
  scale_x_continuous(expand=c(0,1), limits = c(1760, 1990), breaks = c(1775, 1800, 1825, 1850, 1875, 1900, 1925, 1950, 1975, 2000)) + 
  scale_y_continuous(expand=c(0,0), breaks = c(0,8)) +
  theme(plot.margin = unit(c(10,0,0,0), "lines")) +
  ggtitle("") + theme(axis.text.y = element_blank(),
                      axis.title.y = element_blank(),
                      axis.title.x = element_text(vjust=0.1)) +
  annotation_custom(myGrob4, xmin=1980, xmax=2040, ymin=1, ymax=20) 


plot.18S.p2.g <- ggplotGrob(plot.18S.p2)
plot.18S.p2.g$layout$clip[plot.18S.p2.g$layout$name=="panel"] <- "off"


myGrob4 <- grobTree(textGrob("Unclassified Eukaryote", rot=60))


plot.18S.p3 <- ggplot()+
  #geom_rect(data= rects, aes(xmin= xstart, xmax=xend, ymin = ystart, ymax=yend), fill = "red", alpha=0.1) +
  geom_line(data=strat.per.red.depth.melt.p3, aes(Age, value), colour="#21b6a8") +
  geom_area(data=strat.per.red.depth.melt.p3, aes(Age, value), fill="#21b6a8") + 
  geom_bar(data=strat.per.red.depth.melt.p3, aes(Age, value), width=0.1, stat="identity")+
  #facet_grid(~variable, scales = "fixed", labeller = label_wrap_gen(multi_line = TRUE)) +
  #geom_vline(data=zn.df, mapping = aes(xintercept=zone), linetype="dashed", color = "black") + 
  # geom_vline(data=zn.df, mapping = aes(xintercept=1906), linetype="dashed", color = "red") + 
  xlab("Age")+ylab("") +
  coord_flip() +
  theme_new  +
  scale_x_continuous(expand=c(0,1), limits = c(1760, 1990), breaks = c(1775, 1800, 1825, 1850, 1875, 1900, 1925, 1950, 1975, 2000)) + 
  scale_y_continuous(expand=c(0,0)) +
  theme(plot.margin = unit(c(10,0,0,0), "lines")) +
  ggtitle("") + theme(axis.text.y = element_blank(),
                      axis.title.y = element_blank()) +
  annotation_custom(myGrob4, xmin=2005, xmax=2040, ymin=4, ymax=27) 


plot.18S.p3.g <- ggplotGrob(plot.18S.p3)
plot.18S.p3.g$layout$clip[plot.18S.p3.g$layout$name=="panel"] <- "off"

myGrob4 <- grobTree(textGrob("Unclassified Eukaryote", rot=60))


plot.18S.p4 <- ggplot() +
  #geom_rect(data= rects, aes(xmin= xstart, xmax=xend, ymin = ystart, ymax=yend), fill = "red", alpha=0.1) +
  geom_line(data=strat.per.red.depth.melt.p4, aes(Age, value), colour="#21b6a8") +
  geom_area(data=strat.per.red.depth.melt.p4, aes(Age, value), fill="#21b6a8")+ 
  geom_bar(data=strat.per.red.depth.melt.p4, aes(Age, value), width=0.1, stat="identity") +
  #facet_grid(~variable, scales = "fixed", labeller = label_wrap_gen(multi_line = TRUE)) +
  #geom_vline(data=zn.df, mapping = aes(xintercept=zone), linetype="dashed", color = "black") + 
  # geom_vline(data=zn.df, mapping = aes(xintercept=1906), linetype="dashed", color = "red") + 
  xlab("Age")+ylab("") +
  coord_flip() +
  theme_new  +
  scale_x_continuous(expand=c(0,1), limits = c(1760, 1990), breaks = c(1775, 1800, 1825, 1850, 1875, 1900, 1925, 1950, 1975, 2000)) + 
  scale_y_continuous(expand=c(0,0), breaks = c(0,10,20)) +
  theme(plot.margin = unit(c(10,0,0,0), "lines")) +
  ggtitle("") + theme(axis.text.y = element_blank(),
                      axis.title.y = element_blank()) +
  annotation_custom(myGrob3, xmin=zn.df$zone, xmax=zn.df$zone, ymin=-149, ymax=60) +
  annotation_custom(myGrob4, xmin=2005, xmax=2040, ymin=7, ymax=30) 


plot.18S.p4.g <- ggplotGrob(plot.18S.p4)
plot.18S.p4.g$layout$clip[plot.18S.p4.g$layout$name=="panel"] <- "off"


coniss.plot <- plot_grid(plot.18S.p1.g, plot.18S.p2.g, plot.18S.p3.g, plot.18S.p4.g, dendro.18S, align = "hv", axis="tblr", ncol=6, rel_widths = c(0.6, 0.2, 0.4, 0.2, 0.2))

ggsave("coniss.plot.18S.trans.jpeg", width=10, height=7.5, units = "in", bg = "transparent")
coniss.plot
dev.off()


################## INPUT HERE THE ORDERED SUBSET GROUPS ########################

downcore.ps = prune_samples(sample_sums(EukCOILOrder)>0, EukCOILOrder) %>%
  filter_taxa(function(x) sum(x) > 0, TRUE)

downcore.ps.perc <- transform_sample_counts(downcore.ps, function(OTU) OTU/sum(OTU)*100)

downcore.ps.perc.df <- data.frame(otu_table(downcore.ps.perc))

## TAXONOMY to rename OTU labels
downcore.ps.tax.df <- data.frame(tax_table(downcore.ps))
rownames(downcore.ps.tax.df) <- downcore.ps.tax.df$Genus

#Transpose OTU table
downcore.ps.perc.df<-as.data.frame(t(downcore.ps.perc.df))
# make a column of OTU table the taxonomy
downcore.ps.perc.df$otulabels <- downcore.ps.tax.df$Genus
rownames(downcore.ps.perc.df) <- downcore.ps.perc.df$otulabels

downcore.ps.perc.df$`otulabels`=NULL
# transpose back
downcore.ps.perc.df<-as.data.frame(t(downcore.ps.perc.df))

# metadata
downcore.ps.perc.metadata <- data.frame(sample_data(downcore.ps))

rownames(downcore.ps.perc.df) <- downcore.ps.perc.metadata$Age


strat.clust <- downcore.ps.perc.df %>% 
  vegdist(method = "bray", binary = FALSE, diag = FALSE, upper = FALSE,
          na.rm = FALSE) %>% 
  chclust(method = "coniss")

#bstick(strat.clust, 20)
downcore.ps.perc.metadata$Age <- as.numeric(as.character(downcore.ps.perc.metadata$Age))
ct <- cutree(strat.clust, 2)
zn <- which(diff(ct) > 0)
zone <- (downcore.ps.perc.metadata$Age[zn] + downcore.ps.perc.metadata$Age[zn + 1])/2
zn.df <- data.frame(zone)

dendro <- as.dendrogram(strat.clust)
ddata <- dendro_data(dendro, type="rectangle")

depth <- downcore.ps.perc.metadata$Age


ddata$segments->yo
yo$xx=NA
yo$xxend=NA

for(i in 1:length(depth)){
  yo$xx[which(yo$x == i)]=depth[i]
  yo$xxend[which(yo$xend == i)]=depth[i]
  da=depth[i+1]-depth[i]
  yo$xx[yo$x>i & yo$x<i+1]<-depth[i]+(yo$x[yo$x>i & yo$x<i+1]-i)*da
  yo$xxend[yo$xend>i & yo$xend<i+1]<-depth[i]+(yo$xend[yo$xend>i & yo$xend<i+1]-i)*da
} # Please dont ask....
ddata$segments$x<-yo$xx
ddata$segments$xend<-yo$xxend


dendro.COIL = ggplot(segment(ddata)) +
  # geom_vline(data=zn.df, mapping = aes(xintercept=1906), linetype="dashed", color = "red") + 
  # geom_rect(data= rects, aes(xmin= xstart, xmax=xend, ymin = ystart, ymax=yend), fill = "red", alpha=0.1) +
  geom_segment(aes(x=x, y=y, xend=xend, yend=yend)) +
  coord_flip() +
  ylab("Distance")+
  #scale_x_reverse()+
  theme_new + theme(axis.title.y=element_blank(),
                    axis.text.y=element_blank(),
                    axis.ticks.y=element_blank(),
                    axis.line.y = element_blank())   +
  scale_x_continuous(expand=c(0,1), limits = c(1760, 1990)) + 
  scale_y_continuous(expand=c(0,0), limits = c(0, 12)) 



strat.per.red <- downcore.ps.perc.df %>% 
  chooseTaxa(n.occ = 1)


#abund.sel <- strat.per.red %>% 
#  tidyr::gather(key,value) %>% 
#   dplyr::group_by(key) %>% 
#   dplyr::summarise(Sum=sum(value)) %>% 
#   arrange(desc(Sum)) %>% 
#   top_n(5,Sum)

strat.per.red1 <- strat.per.red %>% select("Unclassified_Kingdom_Eukaryota", "Unclassified_Kingdom_Eukaryota_2",  "Ceriodaphnia", "Unclassified_Kingdom_Eukaryota_6")


strat.per.red.depth <- rownames_to_column(strat.per.red) 

colnames(strat.per.red.depth)[1] <- "Age"
#downcore.ps.perc.metadata$Age <- as.character(downcore.ps.perc.metadata$Age)
#strat.per.red.depth <- left_join(strat.per.red.depth, downcore.ps.perc.metadata, by = c("rowname" = "Age"))

strat.per.red.depth <- strat.per.red.depth %>% arrange(Age) 

strat.per.red.depth.melt <- reshape2::melt(strat.per.red.depth, id.vars=c("Age"))

strat.per.red.depth.melt$Age <- as.numeric(strat.per.red.depth.melt$Age)

strat.per.red.depth.melt <- strat.per.red.depth.melt %>% add_row(Age = 1500, variable = "Unclassified_Kingdom_Eukaryota_6", value= 8)

strat.per.red.depth.melt.p1 <-  strat.per.red.depth.melt %>% filter(variable == "Unclassified_Kingdom_Eukaryota")
strat.per.red.depth.melt.p1$variable <- gsub("Unclassified_Kingdom_Eukaryota", "Unclassified Eukaryote", strat.per.red.depth.melt.p1$variable)

strat.per.red.depth.melt.p2 <-  strat.per.red.depth.melt %>% filter(variable == "Unclassified_Kingdom_Eukaryota_2")
strat.per.red.depth.melt.p2$variable <- gsub("Unclassified_Kingdom_Eukaryota_2", "Unclassified Eukaryote", strat.per.red.depth.melt.p2$variable)

strat.per.red.depth.melt.p3 <-  strat.per.red.depth.melt %>% filter(variable == "Ceriodaphnia")
strat.per.red.depth.melt.p3$variable <- gsub("Ceriodaphnia", "Ceriodaphnia sp.", strat.per.red.depth.melt.p3$variable)

strat.per.red.depth.melt.p4 <-  strat.per.red.depth.melt %>% filter(variable == "Unclassified_Kingdom_Eukaryota_6")
strat.per.red.depth.melt.p4$variable <- gsub("Unclassified_Kingdom_Eukaryota_6", "Unclassified Eukaryote", strat.per.red.depth.melt.p4$variable)

myGrob4 <- grobTree(textGrob("Unclassified Eukaryote", rot=60))


plot.COIL.p1 <- ggplot()+
  # geom_rect(data= rects, aes(xmin= xstart, xmax=xend, ymin = ystart, ymax=yend), fill = "red", alpha=0.1) +
  geom_line(data=strat.per.red.depth.melt.p1, aes(Age, value), colour="#c94fb1") +
  geom_area(data=strat.per.red.depth.melt.p1, aes(Age, value), fill="#c94fb1")+ 
  geom_bar(data=strat.per.red.depth.melt.p1, aes(Age, value), width=0.1, stat="identity") +
  #facet_grid(~variable, scales = "fixed", labeller = label_wrap_gen(multi_line = TRUE)) +
  # geom_vline(data=zn.df, mapping = aes(xintercept=zone), linetype="dashed", color = "black") + 
  # geom_vline(data=zn.df, mapping = aes(xintercept=1906), linetype="dashed", color = "red") + 
  xlab("")+ ylab("") +
  coord_flip() +
  theme_new  +
  scale_x_continuous(expand=c(0,1), limits = c(1760, 1990), breaks = c(1775, 1800, 1825, 1850, 1875, 1900, 1925, 1950, 1975, 2000)) + 
  scale_y_continuous(expand=c(0,0)) +
  theme(plot.margin = unit(c(1,0,0,0), "lines")) +
  theme(axis.text.y = element_text(color = c(NA, "black", NA, "black", NA,  "black", NA,  "black", NA, "black"))) +
  ggtitle("b) EukCOI-L") +
  #annotation_custom(myGrob, xmin=1906, xmax=1870, ymin=0, ymax=250) +
  #annotation_custom(myGrob2, xmin=1906, xmax=1906, ymin=0, ymax=230) +
  #annotation_custom(myGrob3, xmin=zn.df$zone, xmax=zn.df$zone, ymin=0, ymax=140)  +
  annotation_custom(myGrob4, xmin=2005, xmax=2040, ymin=2, ymax=27) 


plot.COIL.p1.g <- ggplotGrob(plot.COIL.p1)
plot.COIL.p1.g$layout$clip[plot.COIL.p1.g$layout$name=="panel"] <- "off"
grid.draw(plot.COIL.p1.g)


myGrob4 <- grobTree(textGrob("Unclassified Eukaryote", rot=60))


plot.COIL.p2 <- ggplot()+
  #  geom_rect(data= rects, aes(xmin= xstart, xmax=xend, ymin = ystart, ymax=yend), fill = "red", alpha=0.1) +
  geom_line(data=strat.per.red.depth.melt.p2, aes(Age, value), colour="#c94fb1") +
  geom_area(data=strat.per.red.depth.melt.p2, aes(Age, value), fill="#c94fb1")+ 
  geom_bar(data=strat.per.red.depth.melt.p2, aes(Age, value), width=0.1, stat="identity")+
  #facet_grid(~variable, scales = "fixed", labeller = label_wrap_gen(multi_line = TRUE)) +
  # geom_vline(data=zn.df, mapping = aes(xintercept=zone), linetype="dashed", color = "black") + 
  # geom_vline(data=zn.df, mapping = aes(xintercept=1906), linetype="dashed", color = "red") + 
  xlab("Age")+ylab("Relative Abundance (%)") +
  coord_flip() +
  theme_new  +
  scale_x_continuous(expand=c(0,1), limits = c(1760, 1990), breaks = c(1775, 1800, 1825, 1850, 1875, 1900, 1925, 1950, 1975, 2000)) + 
  scale_y_continuous(expand=c(0,0), breaks = c(0,12)) +
  theme(plot.margin = unit(c(1,0,0,0), "lines")) +
  ggtitle("") + theme(axis.text.y = element_blank(),
                      axis.title.y = element_blank()) +
  annotation_custom(myGrob4, xmin=2005, xmax=2040, ymin=10, ymax=30) 


plot.COIL.p2.g <- ggplotGrob(plot.COIL.p2)
plot.COIL.p2.g$layout$clip[plot.COIL.p2.g$layout$name=="panel"] <- "off"

myGrob4 <- grobTree(textGrob("Ceriodaphnia sp.", rot=60))


plot.COIL.p3 <- ggplot()+
  #  geom_rect(data= rects, aes(xmin= xstart, xmax=xend, ymin = ystart, ymax=yend), fill = "red", alpha=0.1) +
  geom_line(data=strat.per.red.depth.melt.p3, aes(Age, value), colour="#c94fb1") +
  geom_area(data=strat.per.red.depth.melt.p3, aes(Age, value), fill="#c94fb1")+ 
  geom_bar(data=strat.per.red.depth.melt.p3, aes(Age, value), width=0.1, stat="identity")+
  #facet_grid(~variable, scales = "fixed", labeller = label_wrap_gen(multi_line = TRUE)) +
  # geom_vline(data=zn.df, mapping = aes(xintercept=zone), linetype="dashed", color = "black") + 
  # geom_vline(data=zn.df, mapping = aes(xintercept=1906), linetype="dashed", color = "red") + 
  xlab("Age")+ylab("") +
  coord_flip() +
  theme_new  +
  scale_x_continuous(expand=c(0,1), limits = c(1760, 1990), breaks = c(1775, 1800, 1825, 1850, 1875, 1900, 1925, 1950, 1975, 2000)) + 
  scale_y_continuous(expand=c(0,0)) +
  theme(plot.margin = unit(c(1,0,0,0), "lines")) +
  ggtitle("") + theme(axis.text.y = element_blank(),
                      axis.title.y = element_blank()) +
  annotation_custom(myGrob4, xmin=1990, xmax=2040, ymin=3, ymax=20) 


plot.COIL.p3.g <- ggplotGrob(plot.COIL.p3)
plot.COIL.p3.g$layout$clip[plot.COIL.p3.g$layout$name=="panel"] <- "off"

myGrob4 <- grobTree(textGrob("Unclassified Eukaryote", rot=60))


plot.COIL.p4 <- ggplot()+
  #  geom_rect(data= rects, aes(xmin= xstart, xmax=xend, ymin = ystart, ymax=yend), fill = "red", alpha=0.1) +
  geom_line(data=strat.per.red.depth.melt.p4, aes(Age, value), colour="#c94fb1") +
  geom_area(data=strat.per.red.depth.melt.p4, aes(Age, value), fill="#c94fb1")+ 
  geom_bar(data=strat.per.red.depth.melt.p4, aes(Age, value), width=0.1, stat="identity")+
  #facet_grid(~variable, scales = "fixed", labeller = label_wrap_gen(multi_line = TRUE)) +
  # geom_vline(data=zn.df, mapping = aes(xintercept=zone), linetype="dashed", color = "black") + 
  # geom_vline(data=zn.df, mapping = aes(xintercept=1906), linetype="dashed", color = "red") + 
  xlab("Age")+ylab("") +
  coord_flip() +
  theme_new  +
  scale_x_continuous(expand=c(0,1), limits = c(1760, 1990), breaks = c(1775, 1800, 1825, 1850, 1875, 1900, 1925, 1950, 1975, 2000)) + 
  scale_y_continuous(expand=c(0,0), breaks = c(0,8)) +
  theme(plot.margin = unit(c(1,0,0,0), "lines")) +
  ggtitle("") + theme(axis.text.y = element_blank(),
                      axis.title.y = element_blank()) +
  annotation_custom(myGrob4, xmin=2005, xmax=2040, ymin=5, ymax=30) +
  annotation_custom(myGrob3, xmin=zn.df$zone, xmax=zn.df$zone, ymin=-103, ymax=40)


plot.COIL.p4.g <- ggplotGrob(plot.COIL.p4)
plot.COIL.p4.g$layout$clip[plot.COIL.p4.g$layout$name=="panel"] <- "off"


coniss.plot <- plot_grid(plot.COIL.p1.g, plot.COIL.p2.g, plot.COIL.p3.g, plot.COIL.p4.g, dendro.COIL, align = "hv", axis="tblr", ncol=6, rel_widths = c(0.6, 0.2, 0.4, 0.2, 0.2))

ggsave("coniss.plot.COIL.trans.jpeg", width=10, height=7.5, units = "in", bg = "transparent")
coniss.plot
dev.off()



################## INPUT HERE THE ORDERED SUBSET GROUPS ########################

downcore.ps = prune_samples(sample_sums(EukCOISOrder)>0, EukCOISOrder) %>%
  filter_taxa(function(x) sum(x) > 0, TRUE)

downcore.ps.perc <- transform_sample_counts(downcore.ps, function(OTU) OTU/sum(OTU)*100)

downcore.ps.perc.df <- data.frame(otu_table(downcore.ps.perc))

## TAXONOMY to rename OTU labels
downcore.ps.tax.df <- data.frame(tax_table(downcore.ps))
rownames(downcore.ps.tax.df) <- downcore.ps.tax.df$Genus

#Transpose OTU table
downcore.ps.perc.df<-as.data.frame(t(downcore.ps.perc.df))
# make a column of OTU table the taxonomy
downcore.ps.perc.df$otulabels <- downcore.ps.tax.df$Genus
rownames(downcore.ps.perc.df) <- downcore.ps.perc.df$otulabels

downcore.ps.perc.df$`otulabels`=NULL
# transpose back
downcore.ps.perc.df<-as.data.frame(t(downcore.ps.perc.df))

# metadata
downcore.ps.perc.metadata <- data.frame(sample_data(downcore.ps))

rownames(downcore.ps.perc.df) <- downcore.ps.perc.metadata$Age


strat.clust <- downcore.ps.perc.df %>% 
  vegdist(method = "bray", binary = FALSE, diag = FALSE, upper = FALSE,
          na.rm = FALSE) %>% 
  chclust(method = "coniss")

#bstick(strat.clust, 20)
downcore.ps.perc.metadata$Age <- as.numeric(as.character(downcore.ps.perc.metadata$Age))
ct <- cutree(strat.clust, 2)
zn <- which(diff(ct) > 0)
zone <- (downcore.ps.perc.metadata$Age[zn] + downcore.ps.perc.metadata$Age[zn + 1])/2
zn.df <- data.frame(zone)

dendro <- as.dendrogram(strat.clust)
ddata <- dendro_data(dendro, type="rectangle")

depth <- downcore.ps.perc.metadata$Age


ddata$segments->yo
yo$xx=NA
yo$xxend=NA

for(i in 1:length(depth)){
  yo$xx[which(yo$x == i)]=depth[i]
  yo$xxend[which(yo$xend == i)]=depth[i]
  da=depth[i+1]-depth[i]
  yo$xx[yo$x>i & yo$x<i+1]<-depth[i]+(yo$x[yo$x>i & yo$x<i+1]-i)*da
  yo$xxend[yo$xend>i & yo$xend<i+1]<-depth[i]+(yo$xend[yo$xend>i & yo$xend<i+1]-i)*da
} # Please dont ask....
ddata$segments$x<-yo$xx
ddata$segments$xend<-yo$xxend



dendro.COIS = ggplot(segment(ddata)) +
  #geom_vline(data=zn.df, mapping = aes(xintercept=1906), linetype="dashed", color = "red") + 
  #geom_rect(data= rects, aes(xmin= xstart, xmax=xend, ymin = ystart, ymax=yend), fill = "red", alpha=0.1) +
  geom_segment(aes(x=x, y=y, xend=xend, yend=yend)) +
  coord_flip() +
  ylab("Distance")+
  #scale_x_reverse()+
  theme_new + theme(axis.title.y=element_blank(),
                    axis.text.y=element_blank(),
                    axis.ticks.y=element_blank(),
                    axis.line.y = element_blank())   +
  scale_x_continuous(expand=c(0,1), limits = c(1760, 1990)) + 
  scale_y_continuous(expand=c(0,0), limits = c(0, 12)) 


strat.per.red <- downcore.ps.perc.df %>% 
  chooseTaxa(n.occ = 15)


strat.per.red1 <- strat.per.red %>% select("Unclassified_Phylum_Chlorophyta_2", "Unclassified_Kingdom_Eukaryota_14", "Unclassified_Kingdom_Eukaryota_101", "Mychonastes")

strat.per.red.depth <- rownames_to_column(strat.per.red1) 

colnames(strat.per.red.depth)[1] <- "Age"
#downcore.ps.perc.metadata$Age <- as.character(downcore.ps.perc.metadata$Age)
#strat.per.red.depth <- left_join(strat.per.red.depth, downcore.ps.perc.metadata, by = c("rowname" = "Age"))

strat.per.red.depth <- strat.per.red.depth %>% arrange(Age) 

strat.per.red.depth.melt <- reshape2::melt(strat.per.red.depth, id.vars=c("Age"))

strat.per.red.depth.melt$Age <- as.numeric(strat.per.red.depth.melt$Age)

strat.per.red.depth.melt <- strat.per.red.depth.melt %>% add_row(Age = 1500, variable = "Unclassified_Kingdom_Eukaryota_14", value= 10)
strat.per.red.depth.melt <- strat.per.red.depth.melt %>% add_row(Age = 1500, variable = "Unclassified_Kingdom_Eukaryota_101", value= 8)
strat.per.red.depth.melt <- strat.per.red.depth.melt %>% add_row(Age = 1500, variable = "Mychonastes", value= 8)

strat.per.red.depth.melt.p1 <-  strat.per.red.depth.melt %>% filter(variable == "Unclassified_Phylum_Chlorophyta_2")
strat.per.red.depth.melt.p1$variable <- gsub("Unclassified_Phylum_Chlorophyta_2", "Unclassified Chlorophyte", strat.per.red.depth.melt.p1$variable)

strat.per.red.depth.melt.p2 <-  strat.per.red.depth.melt %>% filter(variable == "Unclassified_Kingdom_Eukaryota_14")
strat.per.red.depth.melt.p2$variable <- gsub("Unclassified_Kingdom_Eukaryota_14", "Unclassified Eukaryote", strat.per.red.depth.melt.p2$variable)

strat.per.red.depth.melt.p3 <-  strat.per.red.depth.melt %>% filter(variable == "Unclassified_Kingdom_Eukaryota_101")
strat.per.red.depth.melt.p3$variable <- gsub("Unclassified_Kingdom_Eukaryota_101", "Unclassified Eukaryote", strat.per.red.depth.melt.p3$variable)


strat.per.red.depth.melt.p4 <-  strat.per.red.depth.melt %>% filter(variable == "Mychonastes")
strat.per.red.depth.melt.p4$variable <- gsub("Mychonastes", "Mychonastes sp.", strat.per.red.depth.melt.p4$variable)

myGrob4 <- grobTree(textGrob("Unclassified Chlorophyte", rot=60))

plot.COIS.p1 <- ggplot()+
  #  geom_rect(data= rects, aes(xmin= xstart, xmax=xend, ymin = ystart, ymax=yend), fill = "red", alpha=0.1) +
  geom_line(data=strat.per.red.depth.melt.p1, aes(Age, value), colour="#FFE25C") +
  geom_area(data=strat.per.red.depth.melt.p1, aes(Age, value), fill="#FFE25C")+ 
  geom_bar(data=strat.per.red.depth.melt.p1, aes(Age, value), width=0.1, stat="identity")+
  #facet_grid(~variable, scales = "fixed", labeller = label_wrap_gen(multi_line = TRUE)) +
  #geom_vline(data=zn.df, mapping = aes(xintercept=zone), linetype="dashed", color = "black") + 
  #geom_vline(data=zn.df, mapping = aes(xintercept=1906), linetype="dashed", color = "red") + 
  xlab("")+ ylab("") +
  coord_flip() +
  theme_new  +
  scale_x_continuous(expand=c(0,1), limits = c(1760, 1990), breaks = c(1775, 1800, 1825, 1850, 1875, 1900, 1925, 1950, 1975, 2000)) + 
  scale_y_continuous(expand=c(0,0), breaks = c(0,10,20)) +
  theme(plot.margin = unit(c(1,0,0,0), "lines")) +
  theme(axis.text.y = element_text(color = c(NA, "black", NA, "black", NA,  "black", NA,  "black", NA, "black"))) +
  ggtitle("c) EukCOI-S") +
  #annotation_custom(myGrob, xmin=1906, xmax=1870, ymin=0, ymax=200) +
  #annotation_custom(myGrob2, xmin=1906, xmax=1906, ymin=0, ymax=200) +
  #annotation_custom(myGrob3, xmin=zn.df$zone, xmax=zn.df$zone, ymin=0, ymax=200) +
  annotation_custom(myGrob4, xmin=2008, xmax=2040, ymin=4, ymax=30) 

plot.COIS.p1.g <- ggplotGrob(plot.COIS.p1)
plot.COIS.p1.g$layout$clip[plot.COIS.p1.g$layout$name=="panel"] <- "off"
grid.draw(plot.COIS.p1.g)

myGrob4 <- grobTree(textGrob("Unclassified Eukaryote", rot=60))


plot.COIS.p2 <- ggplot()+
  # geom_rect(data= rects, aes(xmin= xstart, xmax=xend, ymin = ystart, ymax=yend), fill = "red", alpha=0.1) +
  geom_line(data=strat.per.red.depth.melt.p2, aes(Age, value), colour="#FFE25C") +
  geom_area(data=strat.per.red.depth.melt.p2, aes(Age, value), fill="#FFE25C")+ 
  geom_bar(data=strat.per.red.depth.melt.p2, aes(Age, value), width=0.1, stat="identity")+
  #facet_grid(~variable, scales = "fixed", labeller = label_wrap_gen(multi_line = TRUE)) +
  #geom_vline(data=zn.df, mapping = aes(xintercept=zone), linetype="dashed", color = "black") + 
  #geom_vline(data=zn.df, mapping = aes(xintercept=1906), linetype="dashed", color = "red") + 
  xlab("Age")+ylab("") +
  coord_flip() +
  theme_new  +
  scale_x_continuous(expand=c(0,1), limits = c(1760, 1990), breaks = c(1775, 1800, 1825, 1850, 1875, 1900, 1925, 1950, 1975, 2000)) + 
  scale_y_continuous(expand=c(0,0), breaks = c(0,10)) +
  theme(plot.margin = unit(c(1,0,0,0), "lines")) +
  ggtitle("") + theme(axis.text.y = element_blank(),
                      axis.title.y = element_blank()) +
  annotation_custom(myGrob4, xmin=2005, xmax=2040, ymin=5, ymax=30) 


plot.COIS.p2.g <- ggplotGrob(plot.COIS.p2)
plot.COIS.p2.g$layout$clip[plot.COIS.p2.g$layout$name=="panel"] <- "off"

myGrob4 <- grobTree(textGrob("Unclassified Eukaryote", rot=60))


plot.COIS.p3 <- ggplot()+
  #geom_rect(data= rects, aes(xmin= xstart, xmax=xend, ymin = ystart, ymax=yend), fill = "red", alpha=0.1) +
  geom_line(data=strat.per.red.depth.melt.p3, aes(Age, value), colour="#FFE25C") +
  geom_area(data=strat.per.red.depth.melt.p3, aes(Age, value), fill="#FFE25C")+ 
  geom_bar(data=strat.per.red.depth.melt.p3, aes(Age, value), width=0.1, stat="identity")+
  #facet_grid(~variable, scales = "fixed", labeller = label_wrap_gen(multi_line = TRUE)) +
  #geom_vline(data=zn.df, mapping = aes(xintercept=zone), linetype="dashed", color = "black") + 
  #geom_vline(data=zn.df, mapping = aes(xintercept=1906), linetype="dashed", color = "red") + 
  xlab("Age")+ylab("Relative Abundance (%)") +
  coord_flip() +
  theme_new  +
  scale_x_continuous(expand=c(0,1), limits = c(1760, 1990), breaks = c(1775, 1800, 1825, 1850, 1875, 1900, 1925, 1950, 1975, 2000)) + 
  scale_y_continuous(expand=c(0,0), breaks = c(0,8)) +
  theme(plot.margin = unit(c(1,0,0,0), "lines")) +
  ggtitle("") + theme(axis.text.y = element_blank(),
                      axis.title.y = element_blank(),
                      axis.title.x = element_text(hjust = 0.7, vjust=0.1)) +
  annotation_custom(myGrob4, xmin=2005, xmax=2040, ymin=5, ymax=30) 


plot.COIS.p3.g <- ggplotGrob(plot.COIS.p3)
plot.COIS.p3.g$layout$clip[plot.COIS.p3.g$layout$name=="panel"] <- "off"

myGrob4 <- grobTree(textGrob("Mychonastes sp.", rot=60))


plot.COIS.p4 <- ggplot()+
  #  geom_rect(data= rects, aes(xmin= xstart, xmax=xend, ymin = ystart, ymax=yend), fill = "red", alpha=0.1) +
  geom_line(data=strat.per.red.depth.melt.p4, aes(Age, value), colour="#FFE25C") +
  geom_area(data=strat.per.red.depth.melt.p4, aes(Age, value), fill="#FFE25C")+ 
  geom_bar(data=strat.per.red.depth.melt.p4, aes(Age, value), width=0.1, stat="identity")+
  #facet_grid(~variable, scales = "fixed", labeller = label_wrap_gen(multi_line = TRUE)) +
  #geom_vline(data=zn.df, mapping = aes(xintercept=zone), linetype="dashed", color = "black") + 
  #geom_vline(data=zn.df, mapping = aes(xintercept=1906), linetype="dashed", color = "red") + 
  xlab("Age")+ylab("") +
  coord_flip() +
  theme_new  +
  scale_x_continuous(expand=c(0,1), limits = c(1760, 1990), breaks = c(1775, 1800, 1825, 1850, 1875, 1900, 1925, 1950, 1975, 2000)) + 
  scale_y_continuous(expand=c(0,0), breaks = c(0,8)) +
  theme(plot.margin = unit(c(1,0,0,0), "lines")) +
  ggtitle("") + theme(axis.text.y = element_blank(),
                      axis.title.y = element_blank()) +
  annotation_custom(myGrob4, xmin=1990, xmax=2040, ymin=2, ymax=25) +
  annotation_custom(myGrob2, xmin=1906, xmax=1906, ymin=-405, ymax=120) +
  annotation_custom(myGrob3, xmin=zn.df$zone, xmax=zn.df$zone, ymin=-73, ymax=120)


plot.COIS.p4.g <- ggplotGrob(plot.COIS.p4)
plot.COIS.p4.g$layout$clip[plot.COIS.p4.g$layout$name=="panel"] <- "off"


coniss.plot <- plot_grid(plot.COIS.p1.g, plot.COIS.p2.g, plot.COIS.p3.g, plot.COIS.p4.g, dendro.COIS, align = "hv", axis="tblr", ncol=6, rel_widths = c(0.6, 0.2, 0.4, 0.2, 0.2))

ggsave("coniss.plot.COIS.trans.jpeg", width=10, height=7.5, units = "in", bg = "transparent")
coniss.plot
dev.off()

coniss.plot <- plot_grid(plot.18S.p1.g, NULL, plot.18S.p2.g, NULL, plot.18S.p3.g, NULL, plot.18S.p4.g, NULL, dendro.18S, NULL, 
                         plot.COIL.p1.g, NULL, plot.COIL.p2.g, NULL, plot.COIL.p3.g, NULL, plot.COIL.p4.g,NULL,  dendro.COIL, NULL, 
                         plot.COIS.p1.g, NULL, plot.COIS.p2.g,NULL,  plot.COIS.p3.g, NULL, plot.COIS.p4.g, NULL, dendro.COIS, align = "hv", axis="tblr", ncol=29, 
                         rel_widths = c(0.75, -0.1, 0.29, -0.1, 0.375, -0.1, 0.375, -0.1, 0.4, -0.1,
                                        0.5, -0.1, 0.32, -0.1, 0.5, -0.1, 0.29, -0.1, 0.4,-0.08,
                                        0.425, -0.1, 0.30, -0.1, 0.29, -0.1, 0.29, -0.1, 0.4))


ggsave("coniss.plot.euks.jpeg", width=15, height=9, units = "in", bg = "transparent")
coniss.plot
dev.off()

#CONISS Prey


################## INPUT HERE THE ORDERED SUBSET GROUPS ########################

downcore.ps = prune_samples(sample_sums(zoo18SOrder)>0, zoo18SOrder) %>%
  filter_taxa(function(x) sum(x) > 0, TRUE)

downcore.ps.perc <- transform_sample_counts(downcore.ps, function(OTU) OTU/sum(OTU)*100)

downcore.ps.perc.df <- data.frame(otu_table(downcore.ps.perc))

## TAXONOMY to rename OTU labels
downcore.ps.tax.df <- data.frame(tax_table(downcore.ps))
rownames(downcore.ps.tax.df) <- downcore.ps.tax.df$Genus

#Transpose OTU table
downcore.ps.perc.df<-as.data.frame(t(downcore.ps.perc.df))
# make a column of OTU table the taxonomy
downcore.ps.perc.df$otulabels <- downcore.ps.tax.df$Genus
rownames(downcore.ps.perc.df) <- downcore.ps.perc.df$otulabels

downcore.ps.perc.df$`otulabels`=NULL
# transpose back
downcore.ps.perc.df<-as.data.frame(t(downcore.ps.perc.df))

# metadata
downcore.ps.perc.metadata <- data.frame(sample_data(downcore.ps))

rownames(downcore.ps.perc.df) <- downcore.ps.perc.metadata$Age


strat.clust <- downcore.ps.perc.df %>% 
  vegdist(method = "bray", binary = FALSE, diag = FALSE, upper = FALSE,
          na.rm = FALSE) %>% 
  chclust(method = "coniss")

#bstick(strat.clust, 20)
downcore.ps.perc.metadata$Age <- as.numeric(as.character(downcore.ps.perc.metadata$Age))
ct <- cutree(strat.clust, 3)
zn <- which(diff(ct) > 0)
zone <- (downcore.ps.perc.metadata$Age[zn] + downcore.ps.perc.metadata$Age[zn + 1])/2
zn.df <- data.frame(zone)

dendro <- as.dendrogram(strat.clust)
ddata <- dendro_data(dendro, type="rectangle")

depth <- downcore.ps.perc.metadata$Age


ddata$segments->yo
yo$xx=NA
yo$xxend=NA

for(i in 1:length(depth)){
  yo$xx[which(yo$x == i)]=depth[i]
  yo$xxend[which(yo$xend == i)]=depth[i]
  da=depth[i+1]-depth[i]
  yo$xx[yo$x>i & yo$x<i+1]<-depth[i]+(yo$x[yo$x>i & yo$x<i+1]-i)*da
  yo$xxend[yo$xend>i & yo$xend<i+1]<-depth[i]+(yo$xend[yo$xend>i & yo$xend<i+1]-i)*da
} # Please dont ask....
ddata$segments$x<-yo$xx
ddata$segments$xend<-yo$xxend


theme_new <- theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), # remove grids
                   panel.background = element_blank(), axis.line = element_line(colour = "black"),
                   strip.text.x = element_text(size=10, angle=45, hjust=0, vjust=0.05), 
                   # strip.text.x = element_blank(), # Taxa names
                   strip.background = element_blank(),
                   strip.text.y = element_text(angle = 0),
                   legend.position="none",panel.border = element_blank(),
                   axis.text.x=element_text(angle=45,hjust=1, size = 11),
                   axis.text.y = element_text(size=11),
                   axis.title = element_text(size=14),
                   plot.title = element_text(size=15, vjust = 55)) +
  theme(panel.background = element_rect(fill="transparent"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_rect(fill="transparent"),
        plot.background = element_rect(fill="transparent", colour=NA),
        legend.background = element_rect(fill="transparent"),
        legend.key = element_rect(fill = NA, colour = NA))



dendro.18S =ggplot(segment(ddata)) +
  geom_segment(aes(x=x, y=y, xend=xend, yend=yend)) +
  coord_flip() +
  ylab("Distance")+
  #scale_x_reverse()+
  theme_new + theme(axis.title.y=element_blank(),
                    axis.text.y=element_blank(),
                    axis.ticks.y=element_blank(),
                    axis.line.y = element_blank())   +
  scale_x_continuous(expand=c(0,1), limits = c(1760, 1990)) + 
  scale_y_continuous(expand=c(0,0), limits = c(0,15)) +
  ggtitle("")


strat.per.red <- downcore.ps.perc.df %>% 
  chooseTaxa(n.occ = 10)


strat.per.red1 <- strat.per.red %>% select("Elliptio", "Unclassified_Family_Spirostreptida")

colnames(strat.per.red1) <- c("Elliptio sp", "Unclassified_Family_Spirostreptida")
strat.per.red.depth <- rownames_to_column(strat.per.red1) 

colnames(strat.per.red.depth)[1] <- "Age"

strat.per.red.depth <- strat.per.red.depth %>% arrange(Age) 

strat.per.red.depth.melt <- reshape2::melt(strat.per.red.depth, id.vars=c("Age"))

strat.per.red.depth.melt$Age <- as.numeric(strat.per.red.depth.melt$Age)

#strat.per.red.depth.melt <- strat.per.red.depth.melt %>% add_row(Age = 1500, variable = "Leidyana1_1", value= 8)


strat.per.red.depth.melt.p1 <-  strat.per.red.depth.melt %>% filter(variable == "Elliptio sp")
strat.per.red.depth.melt.p1$variable <- gsub("Elliptio sp", "Elliptio sp. (likely Echyridella sp.)", strat.per.red.depth.melt.p1$variable)



strat.per.red.depth.melt.p2 <-  strat.per.red.depth.melt %>% filter(variable == "Unclassified_Family_Spirostreptida")
strat.per.red.depth.melt.p2$variable <- gsub("Unclassified_Family_Spirostreptida", "Unclassified Spirostreptida", strat.per.red.depth.melt.p2$variable)



rects <- data.frame(ystart = c(0), yend = c(Inf), xstart = c(1906), xend = c(1870))

myGrob <- grobTree(rectGrob(gp=gpar(fill="red", alpha=0.1)))

myGrob2 <- grobTree(segmentsGrob(x0=0, x1=1, y0=0, y1=0, default.units="npc", gp=gpar(lty="dashed", col="red", lwd=2)))

myGrob3 <- grobTree(segmentsGrob(x0=0, x1=1, y0=0, y1=0, default.units="npc", gp=gpar(lty="dashed", col="black", lwd=3)))

myGrob5 <- grobTree(segmentsGrob(x0=0, x1=1, y0=0, y1=0, default.units="npc", gp=gpar(lty="dashed", col="grey70", lwd=3)))


myGrob4 <- grobTree(textGrob("Elliptio sp. (likely Echyridella sp.)", rot=60))

plot.18S.p1 <- ggplot()+
  geom_line(data=strat.per.red.depth.melt.p1, aes(Age, value), colour="#21b6a8") +
  geom_area(data=strat.per.red.depth.melt.p1, aes(Age, value), fill="#21b6a8")+ 
  xlab("Age (CE)")+ylab("") +
  coord_flip() +
  theme_new  +
  scale_x_continuous(expand=c(0,1), limits = c(1760, 1990), breaks = c(1775, 1800, 1825, 1850, 1875, 1900, 1925, 1950, 1975, 2000)) + 
  scale_y_continuous(expand=c(0,0)) +
  theme(plot.margin = unit(c(13,0,0,0), "lines")) +
  theme(axis.text.y = element_text(color = c(NA, "black", NA, "black", NA,  "black", NA,  "black", NA, "black"))) +
  ggtitle("a) Prey18S") +
  annotation_custom(myGrob4, xmin=2003, xmax=2070, ymin=10, ymax=50) 


plot.18S.p1.g <- ggplotGrob(plot.18S.p1)
plot.18S.p1.g$layout$clip[plot.18S.p1.g$layout$name=="panel"] <- "off"

myGrob4 <- grobTree(textGrob("Unclassified Spirostreptida", rot=60))


plot.18S.p2 <- ggplot()+
  geom_line(data=strat.per.red.depth.melt.p2, aes(Age, value), colour="#21b6a8") +
  geom_area(data=strat.per.red.depth.melt.p2, aes(Age, value), fill="#21b6a8")+ 
  geom_bar(data=strat.per.red.depth.melt.p2, aes(Age, value), width=0.1, stat="identity")+
  xlab("Age")+ylab("Relative Abundance (%)") +
  coord_flip() +
  theme_new  +
  scale_x_continuous(expand=c(0,1), limits = c(1760, 1990), breaks = c(1775, 1800, 1825, 1850, 1875, 1900, 1925, 1950, 1975, 2000)) + 
  scale_y_continuous(expand=c(0,0), breaks = c(0,10,20)) +
  theme(plot.margin = unit(c(1,0,0,0), "lines")) +
  ggtitle("") + 
  theme(axis.text.y = element_blank(),
        axis.title.y = element_blank(),
        axis.title.x = element_text(hjust=1.15)) +
  annotation_custom(myGrob3, xmin=zn.df$zone[1], xmax=zn.df$zone[1], ymin=-66, ymax=53) +
  annotation_custom(myGrob5, xmin=zn.df$zone[2], xmax=zn.df$zone[2], ymin=-66, ymax=53) +
  annotation_custom(myGrob4, xmin=1998, xmax=2060, ymin=6, ymax=40) 


plot.18S.p2.g <- ggplotGrob(plot.18S.p2)
plot.18S.p2.g$layout$clip[plot.18S.p2.g$layout$name=="panel"] <- "off"



################# INPUT HERE THE ORDERED SUBSET GROUPS ########################

downcore.ps = prune_samples(sample_sums(zoolongOrder)>0, zoolongOrder) %>%
  filter_taxa(function(x) sum(x) > 0, TRUE)

downcore.ps.perc <- transform_sample_counts(downcore.ps, function(OTU) OTU/sum(OTU)*100)

downcore.ps.perc.df <- data.frame(otu_table(downcore.ps.perc))

## TAXONOMY to rename OTU labels
downcore.ps.tax.df <- data.frame(tax_table(downcore.ps))
rownames(downcore.ps.tax.df) <- downcore.ps.tax.df$Genus

#Transpose OTU table
downcore.ps.perc.df<-as.data.frame(t(downcore.ps.perc.df))
# make a column of OTU table the taxonomy
downcore.ps.perc.df$otulabels <- downcore.ps.tax.df$Genus
rownames(downcore.ps.perc.df) <- downcore.ps.perc.df$otulabels

downcore.ps.perc.df$`otulabels`=NULL
# transpose back
downcore.ps.perc.df<-as.data.frame(t(downcore.ps.perc.df))

# metadata
downcore.ps.perc.metadata <- data.frame(sample_data(downcore.ps))

rownames(downcore.ps.perc.df) <- downcore.ps.perc.metadata$Age


strat.clust <- downcore.ps.perc.df %>% 
  vegdist(method = "bray", binary = FALSE, diag = FALSE, upper = FALSE,
          na.rm = FALSE) %>% 
  chclust(method = "coniss")

#bstick(strat.clust, 20)
downcore.ps.perc.metadata$Age <- as.numeric(as.character(downcore.ps.perc.metadata$Age))
ct <- cutree(strat.clust, 2)
zn <- which(diff(ct) > 0)
zone <- (downcore.ps.perc.metadata$Age[zn] + downcore.ps.perc.metadata$Age[zn + 1])/2
zn.df <- data.frame(zone)

dendro <- as.dendrogram(strat.clust)
ddata <- dendro_data(dendro, type="rectangle")

depth <- downcore.ps.perc.metadata$Age


ddata$segments->yo
yo$xx=NA
yo$xxend=NA

for(i in 1:length(depth)){
  yo$xx[which(yo$x == i)]=depth[i]
  yo$xxend[which(yo$xend == i)]=depth[i]
  da=depth[i+1]-depth[i]
  yo$xx[yo$x>i & yo$x<i+1]<-depth[i]+(yo$x[yo$x>i & yo$x<i+1]-i)*da
  yo$xxend[yo$xend>i & yo$xend<i+1]<-depth[i]+(yo$xend[yo$xend>i & yo$xend<i+1]-i)*da
} # Please dont ask....
ddata$segments$x<-yo$xx
ddata$segments$xend<-yo$xxend


dendro.COIL=ggplot(segment(ddata)) +
  geom_segment(aes(x=x, y=y, xend=xend, yend=yend)) +
  coord_flip() +
  ylab("Distance")+
  #scale_x_reverse()+
  theme_new + theme(axis.title.y=element_blank(),
                    axis.text.y=element_blank(),
                    axis.ticks.y=element_blank(),
                    axis.line.y = element_blank())   +
  scale_x_continuous(expand=c(0,1), limits = c(1760, 1990)) + 
  scale_y_continuous(expand=c(0,0), limits = c(0, 15)) 


strat.per.red <- downcore.ps.perc.df %>% 
  chooseTaxa(n.occ = 2)


strat.per.red1 <- strat.per.red %>% select("Ceriodaphnia", "Unclassified_Phylum_Arthropoda_1")


strat.per.red.depth <- rownames_to_column(strat.per.red) 

colnames(strat.per.red.depth)[1] <- "Age"

strat.per.red.depth <- strat.per.red.depth %>% arrange(Age) 

strat.per.red.depth.melt <- reshape2::melt(strat.per.red.depth, id.vars=c("Age"))

strat.per.red.depth.melt$Age <- as.numeric(strat.per.red.depth.melt$Age)

strat.per.red.depth.melt <- strat.per.red.depth.melt %>% add_row(Age = 1500, variable = "Unclassified_Phylum_Arthropoda_1", value= 10)

strat.per.red.depth.melt.p1 <-  strat.per.red.depth.melt %>% filter(variable == "Ceriodaphnia")
strat.per.red.depth.melt.p1$variable <- gsub("Ceriodaphnia", "Ceriodaphnia sp.", strat.per.red.depth.melt.p1$variable)

strat.per.red.depth.melt.p2 <-  strat.per.red.depth.melt %>% filter(variable == "Unclassified_Phylum_Arthropoda_1")
strat.per.red.depth.melt.p2$variable <- gsub("Unclassified_Phylum_Arthropoda_1", "Unclassified Arthropoda", strat.per.red.depth.melt.p2$variable)

myGrob4 <- grobTree(textGrob("Ceriodaphnia sp.", rot=60))


plot.COIL.p1 <- ggplot()+
  geom_line(data=strat.per.red.depth.melt.p1, aes(Age, value), colour="#c94fb1") +
  geom_area(data=strat.per.red.depth.melt.p1, aes(Age, value), fill="#c94fb1")+ 
  geom_bar(data=strat.per.red.depth.melt.p1, aes(Age, value), width=0.1, stat="identity")+
  xlab("")+ylab("") +
  coord_flip() +
  theme_new  +
  scale_x_continuous(expand=c(0,1), limits = c(1760, 1990), breaks = c(1775, 1800, 1825, 1850, 1875, 1900, 1925, 1950, 1975, 2000)) + 
  scale_y_continuous(expand=c(0,0)) +
  theme(plot.margin = unit(c(1,0,0,0), "lines")) +
  theme(axis.text.y = element_text(color = c(NA, "black", NA, "black", NA,  "black", NA,  "black", NA, "black"))) +
  ggtitle("b) PreyCOI-L") +
  annotation_custom(myGrob4, xmin=2003, xmax=2030, ymin=2, ymax=27) 


plot.COIL.p1.g <- ggplotGrob(plot.COIL.p1)
plot.COIL.p1.g$layout$clip[plot.COIL.p1.g$layout$name=="panel"] <- "off"



myGrob4 <- grobTree(textGrob("Unclassified Arthropoda", rot=60))


plot.COIL.p2 <- ggplot()+
  geom_line(data=strat.per.red.depth.melt.p2, aes(Age, value), colour="#c94fb1") +
  geom_area(data=strat.per.red.depth.melt.p2, aes(Age, value), fill="#c94fb1")+ 
  geom_bar(data=strat.per.red.depth.melt.p2, aes(Age, value), width=0.1, stat="identity")+
  xlab("Age")+ylab("Relative Abundance (%)") +
  coord_flip() +
  theme_new  +
  scale_x_continuous(expand=c(0,1), limits = c(1760, 1990), breaks = c(1775, 1800, 1825, 1850, 1875, 1900, 1925, 1950, 1975, 2000)) + 
  scale_y_continuous(expand=c(0,0), breaks = c(0, 10)) +
  theme(plot.margin = unit(c(1,0,0,0), "lines")) +
  ggtitle("") + theme(axis.text.y = element_blank(),
                      axis.title.y = element_blank(),
                      axis.title.x = element_text(hjust=1.15)) +
  annotation_custom(myGrob4, xmin=2008, xmax=2040, ymin=4, ymax=20) +
  annotation_custom(myGrob3, xmin=zn.df$zone[1], xmax=zn.df$zone[1], ymin=-64, ymax=43) +
  annotation_custom(myGrob5, xmin=zn.df$zone[2], xmax=zn.df$zone[2], ymin=-64, ymax=43)

plot.COIL.p2.g <- ggplotGrob(plot.COIL.p2)
plot.COIL.p2.g$layout$clip[plot.COIL.p2.g$layout$name=="panel"] <- "off"


################## INPUT HERE THE ORDERED SUBSET GROUPS ########################

downcore.ps = prune_samples(sample_sums(zooshortOrder)>0, zooshortOrder) %>%
  filter_taxa(function(x) sum(x) > 0, TRUE)

downcore.ps.perc <- transform_sample_counts(downcore.ps, function(OTU) OTU/sum(OTU)*100)

downcore.ps.perc.df <- data.frame(otu_table(downcore.ps.perc))

## TAXONOMY to rename OTU labels
downcore.ps.tax.df <- data.frame(tax_table(downcore.ps))
rownames(downcore.ps.tax.df) <- downcore.ps.tax.df$Genus

#Transpose OTU table
downcore.ps.perc.df<-as.data.frame(t(downcore.ps.perc.df))
# make a column of OTU table the taxonomy
downcore.ps.perc.df$otulabels <- downcore.ps.tax.df$Genus
rownames(downcore.ps.perc.df) <- downcore.ps.perc.df$otulabels

downcore.ps.perc.df$`otulabels`=NULL
# transpose back
downcore.ps.perc.df<-as.data.frame(t(downcore.ps.perc.df))

# metadata
downcore.ps.perc.metadata <- data.frame(sample_data(downcore.ps))

rownames(downcore.ps.perc.df) <- downcore.ps.perc.metadata$Age


strat.clust <- downcore.ps.perc.df %>% 
  vegdist(method = "bray", binary = FALSE, diag = FALSE, upper = FALSE,
          na.rm = FALSE) %>% 
  chclust(method = "coniss")

#bstick(strat.clust, 20)
downcore.ps.perc.metadata$Age <- as.numeric(as.character(downcore.ps.perc.metadata$Age))
ct <- cutree(strat.clust, 2)
zn <- which(diff(ct) > 0)
zone <- (downcore.ps.perc.metadata$Age[zn] + downcore.ps.perc.metadata$Age[zn + 1])/2
zn.df <- data.frame(zone)

dendro <- as.dendrogram(strat.clust)
ddata <- dendro_data(dendro, type="rectangle")

depth <- downcore.ps.perc.metadata$Age


ddata$segments->yo
yo$xx=NA
yo$xxend=NA

for(i in 1:length(depth)){
  yo$xx[which(yo$x == i)]=depth[i]
  yo$xxend[which(yo$xend == i)]=depth[i]
  da=depth[i+1]-depth[i]
  yo$xx[yo$x>i & yo$x<i+1]<-depth[i]+(yo$x[yo$x>i & yo$x<i+1]-i)*da
  yo$xxend[yo$xend>i & yo$xend<i+1]<-depth[i]+(yo$xend[yo$xend>i & yo$xend<i+1]-i)*da
} # Please dont ask....
ddata$segments$x<-yo$xx
ddata$segments$xend<-yo$xxend


dendro.COIS = ggplot(segment(ddata)) +
  geom_segment(aes(x=x, y=y, xend=xend, yend=yend)) +
  coord_flip() +
  ylab("Distance")+
  #scale_x_reverse()+
  theme_new + theme(axis.title.y=element_blank(),
                    axis.text.y=element_blank(),
                    axis.ticks.y=element_blank(),
                    axis.line.y = element_blank())   +
  scale_x_continuous(expand=c(0,1), limits = c(1760, 1990)) + 
  scale_y_continuous(expand=c(0,0), limits = c(0, 15)) 


strat.per.red <- downcore.ps.perc.df %>% 
  chooseTaxa(n.occ = 3)


strat.per.red1 <- strat.per.red %>% select("Unclassified_Phylum_Arthropoda_6", "Unclassified_Phylum_Arthropoda_10")

strat.per.red.depth <- rownames_to_column(strat.per.red1) 

colnames(strat.per.red.depth)[1] <- "Age"

strat.per.red.depth <- strat.per.red.depth %>% arrange(Age) 

strat.per.red.depth.melt <- reshape2::melt(strat.per.red.depth, id.vars=c("Age"))

strat.per.red.depth.melt$Age <- as.numeric(strat.per.red.depth.melt$Age)


strat.per.red.depth.melt.p1 <-  strat.per.red.depth.melt %>% filter(variable == "Unclassified_Phylum_Arthropoda_6")
strat.per.red.depth.melt.p1$variable <- gsub("Unclassified_Phylum_Arthropoda_6", "Unclassified Arthropoda", strat.per.red.depth.melt.p1$variable)

strat.per.red.depth.melt.p2 <-  strat.per.red.depth.melt %>% filter(variable == "Unclassified_Phylum_Arthropoda_10")
strat.per.red.depth.melt.p2$variable <- gsub("Unclassified_Phylum_Arthropoda_10", "Unclassified Arthropoda", strat.per.red.depth.melt.p2$variable)

myGrob4 <- grobTree(textGrob("Unclassified Arthropoda", rot=60))


plot.COIS.p1 <- ggplot()+
  geom_line(data=strat.per.red.depth.melt.p1, aes(Age, value), colour="#FFE25C") +
  geom_area(data=strat.per.red.depth.melt.p1, aes(Age, value), fill="#FFE25C")+ 
  geom_bar(data=strat.per.red.depth.melt.p1, aes(Age, value), width=0.1, stat="identity")+
  xlab("")+ylab("") +
  coord_flip() +
  theme_new  +
  scale_x_continuous(expand=c(0,1), limits = c(1760, 1990), breaks = c(1775, 1800, 1825, 1850, 1875, 1900, 1925, 1950, 1975, 2000)) + 
  scale_y_continuous(expand=c(0,0), breaks = c(0,10,20)) +
  theme(plot.margin = unit(c(1,0,0,0), "lines")) +
  theme(axis.text.y = element_text(color = c(NA, "black", NA, "black", NA,  "black", NA,  "black", NA, "black"))) +
  ggtitle("c) PreyCOI-S") +
  annotation_custom(myGrob4, xmin=2008, xmax=2040, ymin=3, ymax=30) 


plot.COIS.p1.g <- ggplotGrob(plot.COIS.p1)
plot.COIS.p1.g$layout$clip[plot.COIS.p1.g$layout$name=="panel"] <- "off"


myGrob4 <- grobTree(textGrob("Unclassified Arthropoda", rot=60))


plot.COIS.p2 <- ggplot()+
  geom_line(data=strat.per.red.depth.melt.p2, aes(Age, value), colour="#FFE25C") +
  geom_area(data=strat.per.red.depth.melt.p2, aes(Age, value), fill="#FFE25C")+ 
  geom_bar(data=strat.per.red.depth.melt.p2, aes(Age, value), width=0.1, stat="identity")+
  xlab("Age")+ylab("Rel. Abundance (%)") +
  coord_flip() +
  theme_new  +
  scale_x_continuous(expand=c(0,1), limits = c(1760, 1990), breaks = c(1775, 1800, 1825, 1850, 1875, 1900, 1925, 1950, 1975, 2000)) + 
  scale_y_continuous(expand=c(0,0), breaks = c(0, 10, 20)) +
  theme(plot.margin = unit(c(1,0,0,0), "lines")) +
  ggtitle("") + theme(axis.text.y = element_blank(),
                      axis.title.y = element_blank(),
                      axis.title.x = element_text(hjust=0.7)) +
  annotation_custom(myGrob4, xmin=2008, xmax=2040, ymin=5, ymax=30) +
  annotation_custom(myGrob2, xmin=1906, xmax=1906, ymin=-330, ymax=120) +
  annotation_custom(myGrob3, xmin=zn.df$zone[1], xmax=zn.df$zone[1], ymin=-24, ymax=66) +
  annotation_custom(myGrob, xmin=1906, xmax=1870, ymin=-330, ymax=700)

plot.COIS.p2.g <- ggplotGrob(plot.COIS.p2)
plot.COIS.p2.g$layout$clip[plot.COIS.p2.g$layout$name=="panel"] <- "off"



coniss.plot <- plot_grid(plot.18S.p1.g, NULL, plot.18S.p2.g, NULL, dendro.18S, NULL,
                         plot.COIL.p1.g, NULL, plot.COIL.p2.g, NULL, dendro.COIL, NULL,
                         plot.COIS.p1.g, NULL, plot.COIS.p2.g, NULL, dendro.COIS, align = "hv", axis="tblr", ncol=17, 
                         rel_widths = c(0.5, -0.1, 0.3, -0.1, 0.4, -0.1, 
                                        0.7, -0.1, 0.25, -0.1, 0.4, -0.08,
                                        0.25, -0.1, 0.3, -0.1, 0.4))


ggsave("coniss.plot.invert.jpeg", width=12, height=9, units = "in", bg = "transparent")
coniss.plot
dev.off()


################## INPUT HERE THE ORDERED SUBSET GROUPS ########################

downcore.ps = prune_samples(sample_sums(Algae18SOrder)>0, Algae18SOrder) %>%
  filter_taxa(function(x) sum(x) > 0, TRUE)

downcore.ps.perc <- transform_sample_counts(downcore.ps, function(OTU) OTU/sum(OTU)*100)

downcore.ps.perc.df <- data.frame(otu_table(downcore.ps.perc))

## TAXONOMY to rename OTU labels
downcore.ps.tax.df <- data.frame(tax_table(downcore.ps))
rownames(downcore.ps.tax.df) <- downcore.ps.tax.df$Genus

#Transpose OTU table
downcore.ps.perc.df<-as.data.frame(t(downcore.ps.perc.df))
# make a column of OTU table the taxonomy
downcore.ps.perc.df$otulabels <- downcore.ps.tax.df$Genus
rownames(downcore.ps.perc.df) <- downcore.ps.perc.df$otulabels

downcore.ps.perc.df$`otulabels`=NULL
# transpose back
downcore.ps.perc.df<-as.data.frame(t(downcore.ps.perc.df))

# metadata
downcore.ps.perc.metadata <- data.frame(sample_data(downcore.ps))

rownames(downcore.ps.perc.df) <- downcore.ps.perc.metadata$Age


strat.clust <- downcore.ps.perc.df %>% 
  vegdist(method = "bray", binary = FALSE, diag = FALSE, upper = FALSE,
          na.rm = FALSE) %>% 
  chclust(method = "coniss")

#bstick(strat.clust, 20)
downcore.ps.perc.metadata$Age <- as.numeric(as.character(downcore.ps.perc.metadata$Age))
ct <- cutree(strat.clust, 2)
zn <- which(diff(ct) > 0)
zone <- (downcore.ps.perc.metadata$Age[zn] + downcore.ps.perc.metadata$Age[zn + 1])/2
zn.df <- data.frame(zone)

dendro <- as.dendrogram(strat.clust)
ddata <- dendro_data(dendro, type="rectangle")

depth <- downcore.ps.perc.metadata$Age


ddata$segments->yo
yo$xx=NA
yo$xxend=NA

for(i in 1:length(depth)){
  yo$xx[which(yo$x == i)]=depth[i]
  yo$xxend[which(yo$xend == i)]=depth[i]
  da=depth[i+1]-depth[i]
  yo$xx[yo$x>i & yo$x<i+1]<-depth[i]+(yo$x[yo$x>i & yo$x<i+1]-i)*da
  yo$xxend[yo$xend>i & yo$xend<i+1]<-depth[i]+(yo$xend[yo$xend>i & yo$xend<i+1]-i)*da
} # Please dont ask....
ddata$segments$x<-yo$xx
ddata$segments$xend<-yo$xxend


theme_new <- theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), # remove grids
                   panel.background = element_blank(), axis.line = element_line(colour = "black"),
                   strip.text.x = element_text(size=10, angle=45, hjust=0, vjust=0.05), 
                   # strip.text.x = element_blank(), # Taxa names
                   strip.background = element_blank(),
                   strip.text.y = element_text(angle = 0),
                   legend.position="none",panel.border = element_blank(),
                   axis.text.x=element_text(angle=45,hjust=1, size = 11),
                   axis.text.y = element_text(size=11),
                   axis.title = element_text(size=14),
                   plot.title = element_text(size=15, vjust = 45)) +
  theme(panel.background = element_rect(fill="transparent"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_rect(fill="transparent"),
        plot.background = element_rect(fill="transparent", colour=NA),
        legend.background = element_rect(fill="transparent"),
        legend.key = element_rect(fill = NA, colour = NA))



dendro.18S =ggplot(segment(ddata)) +
  geom_segment(aes(x=x, y=y, xend=xend, yend=yend)) +
  coord_flip() +
  ylab("Distance")+
  #scale_x_reverse()+
  theme_new + theme(axis.title.y=element_blank(),
                    axis.text.y=element_blank(),
                    axis.ticks.y=element_blank(),
                    axis.line.y = element_blank())   +
  scale_x_continuous(expand=c(0,1), limits = c(1760, 1990)) + 
  scale_y_continuous(expand=c(0,0), limits = c(0,15)) +
  ggtitle("")


strat.per.red <- downcore.ps.perc.df %>% 
  chooseTaxa(n.occ = 5)


strat.per.red1 <- strat.per.red %>% select("Choricystis", "Mychonastes", "Desmodesmus", "DinophyceaeX")


colnames(strat.per.red1) <- c("Choricystis", "Mychonastes", "Desmodesmus", "DinophyceaeX")

strat.per.red.depth <- rownames_to_column(strat.per.red1) 

colnames(strat.per.red.depth)[1] <- "Age"

strat.per.red.depth <- strat.per.red.depth %>% arrange(Age) 

strat.per.red.depth.melt <- reshape2::melt(strat.per.red.depth, id.vars=c("Age"))

strat.per.red.depth.melt$Age <- as.numeric(strat.per.red.depth.melt$Age)

#strat.per.red.depth.melt <- strat.per.red.depth.melt %>% add_row(Age = 1500, variable = "Leidyana1_1", value= 8)


strat.per.red.depth.melt.p1 <-  strat.per.red.depth.melt %>% filter(variable == "Choricystis")
strat.per.red.depth.melt.p1$variable <- gsub("Choricystis", "Choricystis sp.", strat.per.red.depth.melt.p1$variable)

strat.per.red.depth.melt.p2 <-  strat.per.red.depth.melt %>% filter(variable == "Mychonastes")
strat.per.red.depth.melt.p2$variable <- gsub("Mychonastes", "Mychonastes sp.", strat.per.red.depth.melt.p2$variable)

strat.per.red.depth.melt.p3 <-  strat.per.red.depth.melt %>% filter(variable == "Desmodesmus")
strat.per.red.depth.melt.p3$variable <- gsub("Desmodesmus", "Desmodesmus sp.", strat.per.red.depth.melt.p3$variable)

strat.per.red.depth.melt.p4 <-  strat.per.red.depth.melt %>% filter(variable == "DinophyceaeX")
strat.per.red.depth.melt.p4$variable <- gsub("DinophyceaeX", "Dinophyceae sp.", strat.per.red.depth.melt.p4$variable)


rects <- data.frame(ystart = c(0), yend = c(Inf), xstart = c(1906), xend = c(1870))

myGrob <- grobTree(rectGrob(gp=gpar(fill="red", alpha=0.1)))

myGrob2 <- grobTree(segmentsGrob(x0=0, x1=1, y0=0, y1=0, default.units="npc", gp=gpar(lty="dashed", col="red", lwd=2)))

myGrob3 <- grobTree(segmentsGrob(x0=0, x1=1, y0=0, y1=0, default.units="npc", gp=gpar(lty="dashed", col="black", lwd=3)))

myGrob5 <- grobTree(segmentsGrob(x0=0, x1=1, y0=0, y1=0, default.units="npc", gp=gpar(lty="dashed", col="grey70", lwd=3)))


myGrob4 <- grobTree(textGrob("Choricystis sp.", rot=60))

plot.18S.p1 <- ggplot()+
  geom_line(data=strat.per.red.depth.melt.p1, aes(Age, value), colour="#21b6a8") +
  geom_area(data=strat.per.red.depth.melt.p1, aes(Age, value), fill="#21b6a8")+ 
  xlab("Age (CE)")+ylab("") +
  coord_flip() +
  theme_new  +
  scale_x_continuous(expand=c(0,1), limits = c(1760, 1990), breaks = c(1775, 1800, 1825, 1850, 1875, 1900, 1925, 1950, 1975, 2000)) + 
  scale_y_continuous(expand=c(0,0)) +
  theme(plot.margin = unit(c(13,0,0,0), "lines")) +
  theme(axis.text.y = element_text(color = c(NA, "black", NA, "black", NA,  "black", NA,  "black", NA, "black"))) +
  ggtitle("a) Algae18S") +
  annotation_custom(myGrob4, xmin=1990, xmax=2040, ymin=10, ymax=50) 


plot.18S.p1.g <- ggplotGrob(plot.18S.p1)
plot.18S.p1.g$layout$clip[plot.18S.p1.g$layout$name=="panel"] <- "off"


myGrob4 <- grobTree(textGrob("Mychonastes sp.", rot=60))

plot.18S.p2 <- ggplot()+
  geom_line(data=strat.per.red.depth.melt.p2, aes(Age, value), colour="#21b6a8") +
  geom_area(data=strat.per.red.depth.melt.p2, aes(Age, value), fill="#21b6a8")+ 
  xlab("Age (CE)")+ylab("Relative Abundance (%)") +
  coord_flip() +
  theme_new  +
  scale_x_continuous(expand=c(0,1), limits = c(1760, 1990), breaks = c(1775, 1800, 1825, 1850, 1875, 1900, 1925, 1950, 1975, 2000)) + 
  scale_y_continuous(expand=c(0,0)) +
  theme(plot.margin = unit(c(13,0,0,0), "lines")) +
  theme(axis.text.y = element_text(color = c(NA, "black", NA, "black", NA,  "black", NA,  "black", NA, "black"))) +
  theme(axis.text.y = element_blank(),
        axis.title.y = element_blank()) +
  ggtitle("") +
  annotation_custom(myGrob4, xmin=1999, xmax=2040, ymin=10, ymax=50) 


plot.18S.p2.g <- ggplotGrob(plot.18S.p2)
plot.18S.p2.g$layout$clip[plot.18S.p2.g$layout$name=="panel"] <- "off"


myGrob4 <- grobTree(textGrob("Desmodesmus sp.", rot=60))

plot.18S.p3 <- ggplot()+
  geom_line(data=strat.per.red.depth.melt.p3, aes(Age, value), colour="#21b6a8") +
  geom_area(data=strat.per.red.depth.melt.p3, aes(Age, value), fill="#21b6a8")+ 
  xlab("Age (CE)")+ylab("") +
  coord_flip() +
  theme_new  +
  scale_x_continuous(expand=c(0,1), limits = c(1760, 1990), breaks = c(1775, 1800, 1825, 1850, 1875, 1900, 1925, 1950, 1975, 2000)) + 
  scale_y_continuous(expand=c(0,0), breaks = c(0,20, 40)) +
  theme(plot.margin = unit(c(13,0,0,0), "lines")) +
  theme(axis.text.y = element_text(color = c(NA, "black", NA, "black", NA,  "black", NA,  "black", NA, "black"))) +
  theme(axis.text.y = element_blank(),
        axis.title.y = element_blank()) +
  ggtitle("") +
  annotation_custom(myGrob4, xmin=2003, xmax=2040, ymin=35, ymax=60) 


plot.18S.p3.g <- ggplotGrob(plot.18S.p3)
plot.18S.p3.g$layout$clip[plot.18S.p3.g$layout$name=="panel"] <- "off"


myGrob4 <- grobTree(textGrob("Dinophyceae sp.", rot=60))


plot.18S.p4 <- ggplot()+
  geom_line(data=strat.per.red.depth.melt.p4, aes(Age, value), colour="#21b6a8") +
  geom_area(data=strat.per.red.depth.melt.p4, aes(Age, value), fill="#21b6a8")+ 
  geom_bar(data=strat.per.red.depth.melt.p4, aes(Age, value), width=0.1, stat="identity")+
  xlab("Age")+ylab("") +
  coord_flip() +
  theme_new  +
  scale_x_continuous(expand=c(0,1), limits = c(1760, 1990), breaks = c(1775, 1800, 1825, 1850, 1875, 1900, 1925, 1950, 1975, 2000)) + 
  scale_y_continuous(expand=c(0,0), breaks = c(0,20)) +
  theme(plot.margin = unit(c(1,0,0,0), "lines")) +
  ggtitle("") + 
  theme(axis.text.y = element_blank(),
        axis.title.y = element_blank(),
        axis.title.x = element_text(hjust=1.15)) +
  annotation_custom(myGrob3, xmin=zn.df$zone[1], xmax=zn.df$zone[1], ymin=-480, ymax=155) +
  #annotation_custom(myGrob3, xmin=zn.df$zone[2], xmax=zn.df$zone[2], ymin=-66, ymax=53) +
  annotation_custom(myGrob4, xmin=1998, xmax=2040, ymin=26, ymax=60) 


plot.18S.p4.g <- ggplotGrob(plot.18S.p4)
plot.18S.p4.g$layout$clip[plot.18S.p4.g$layout$name=="panel"] <- "off"


################## INPUT HERE THE ORDERED SUBSET GROUPS ########################

downcore.ps = prune_samples(sample_sums(AlgaeLongOrder)>0, AlgaeLongOrder) %>%
  filter_taxa(function(x) sum(x) > 0, TRUE)

downcore.ps.perc <- transform_sample_counts(downcore.ps, function(OTU) OTU/sum(OTU)*100)

downcore.ps.perc.df <- data.frame(otu_table(downcore.ps.perc))

## TAXONOMY to rename OTU labels
downcore.ps.tax.df <- data.frame(tax_table(downcore.ps))
rownames(downcore.ps.tax.df) <- downcore.ps.tax.df$Genus

#Transpose OTU table
downcore.ps.perc.df<-as.data.frame(t(downcore.ps.perc.df))
# make a column of OTU table the taxonomy
downcore.ps.perc.df$otulabels <- downcore.ps.tax.df$Genus
rownames(downcore.ps.perc.df) <- downcore.ps.perc.df$otulabels

downcore.ps.perc.df$`otulabels`=NULL
# transpose back
downcore.ps.perc.df<-as.data.frame(t(downcore.ps.perc.df))

# metadata
downcore.ps.perc.metadata <- data.frame(sample_data(downcore.ps))

rownames(downcore.ps.perc.df) <- downcore.ps.perc.metadata$Age


strat.clust <- downcore.ps.perc.df %>% 
  vegdist(method = "bray", binary = FALSE, diag = FALSE, upper = FALSE,
          na.rm = FALSE) %>% 
  chclust(method = "coniss")

#bstick(strat.clust, 20)
downcore.ps.perc.metadata$Age <- as.numeric(as.character(downcore.ps.perc.metadata$Age))
ct <- cutree(strat.clust, 3)
zn <- which(diff(ct) > 0)
zone <- (downcore.ps.perc.metadata$Age[zn] + downcore.ps.perc.metadata$Age[zn + 1])/2
zn.df <- data.frame(zone)

dendro <- as.dendrogram(strat.clust)
ddata <- dendro_data(dendro, type="rectangle")

depth <- downcore.ps.perc.metadata$Age


ddata$segments->yo
yo$xx=NA
yo$xxend=NA

for(i in 1:length(depth)){
  yo$xx[which(yo$x == i)]=depth[i]
  yo$xxend[which(yo$xend == i)]=depth[i]
  da=depth[i+1]-depth[i]
  yo$xx[yo$x>i & yo$x<i+1]<-depth[i]+(yo$x[yo$x>i & yo$x<i+1]-i)*da
  yo$xxend[yo$xend>i & yo$xend<i+1]<-depth[i]+(yo$xend[yo$xend>i & yo$xend<i+1]-i)*da
} # Please dont ask....
ddata$segments$x<-yo$xx
ddata$segments$xend<-yo$xxend


dendro.COIL=ggplot(segment(ddata)) +
  geom_segment(aes(x=x, y=y, xend=xend, yend=yend)) +
  coord_flip() +
  ylab("Distance")+
  #scale_x_reverse()+
  theme_new + theme(axis.title.y=element_blank(),
                    axis.text.y=element_blank(),
                    axis.ticks.y=element_blank(),
                    axis.line.y = element_blank())   +
  scale_x_continuous(expand=c(0,1), limits = c(1760, 1990)) + 
  scale_y_continuous(expand=c(0,0), limits = c(0, 15)) 


strat.per.red <- downcore.ps.perc.df %>% 
  chooseTaxa(n.occ = 5)


strat.per.red1 <- strat.per.red %>% select("Unclassified_Phylum_Chlorophyta", "Unclassified_Phylum_Chlorophyta_1",  "Unclassified_Order_Sphaeropleales", "Nannochloropsis")


strat.per.red.depth <- rownames_to_column(strat.per.red) 

colnames(strat.per.red.depth)[1] <- "Age"

strat.per.red.depth <- strat.per.red.depth %>% arrange(Age) 

strat.per.red.depth.melt <- reshape2::melt(strat.per.red.depth, id.vars=c("Age"))

strat.per.red.depth.melt$Age <- as.numeric(strat.per.red.depth.melt$Age)


strat.per.red.depth.melt.p1 <-  strat.per.red.depth.melt %>% filter(variable == "Unclassified_Phylum_Chlorophyta")
strat.per.red.depth.melt.p1$variable <- gsub("Unclassified_Phylum_Chlorophyta", "Unclassified Chlorophyta", strat.per.red.depth.melt.p1$variable)

strat.per.red.depth.melt.p2 <-  strat.per.red.depth.melt %>% filter(variable == "Unclassified_Phylum_Chlorophyta_1")
strat.per.red.depth.melt.p2$variable <- gsub("Unclassified_Phylum_Chlorophyta_1", "Unclassified Chlorophyta", strat.per.red.depth.melt.p2$variable)

strat.per.red.depth.melt.p3 <-  strat.per.red.depth.melt %>% filter(variable == "Unclassified_Order_Sphaeropleales")
strat.per.red.depth.melt.p3$variable <- gsub("Unclassified_Order_Sphaeropleales", "Unclassified Sphaeropleales", strat.per.red.depth.melt.p3$variable)

strat.per.red.depth.melt.p4 <-  strat.per.red.depth.melt %>% filter(variable == "Nannochloropsis")
strat.per.red.depth.melt.p4$variable <- gsub("Nannochloropsis", "Nannochloropsis sp.", strat.per.red.depth.melt.p4$variable)




myGrob4 <- grobTree(textGrob("Unclassified Chlorophyta", rot=60))


plot.COIL.p1 <- ggplot()+
  geom_line(data=strat.per.red.depth.melt.p1, aes(Age, value), colour="#c94fb1") +
  geom_area(data=strat.per.red.depth.melt.p1, aes(Age, value), fill="#c94fb1")+ 
  geom_bar(data=strat.per.red.depth.melt.p1, aes(Age, value), width=0.1, stat="identity")+
  xlab("")+ylab("") +
  coord_flip() +
  theme_new  +
  scale_x_continuous(expand=c(0,1), limits = c(1760, 1990), breaks = c(1775, 1800, 1825, 1850, 1875, 1900, 1925, 1950, 1975, 2000)) + 
  scale_y_continuous(expand=c(0,0)) +
  theme(plot.margin = unit(c(1,0,0,0), "lines")) +
  theme(axis.text.y = element_text(color = c(NA, "black", NA, "black", NA,  "black", NA,  "black", NA, "black"))) +
  ggtitle("b) AlgaeCOI-L") +
  annotation_custom(myGrob4, xmin=1998, xmax=2060, ymin=35, ymax=70) 


plot.COIL.p1.g <- ggplotGrob(plot.COIL.p1)
plot.COIL.p1.g$layout$clip[plot.COIL.p1.g$layout$name=="panel"] <- "off"


myGrob4 <- grobTree(textGrob("Unclassified Chlorophyta", rot=60))

plot.COIL.p2 <- ggplot()+
  geom_line(data=strat.per.red.depth.melt.p2, aes(Age, value), colour="#c94fb1") +
  geom_area(data=strat.per.red.depth.melt.p2, aes(Age, value), fill="#c94fb1")+ 
  geom_bar(data=strat.per.red.depth.melt.p2, aes(Age, value), width=0.1, stat="identity")+
  xlab("")+ylab("Relative Abundance (%)") +
  coord_flip() +
  theme_new  +
  scale_x_continuous(expand=c(0,1), limits = c(1760, 1990), breaks = c(1775, 1800, 1825, 1850, 1875, 1900, 1925, 1950, 1975, 2000)) + 
  scale_y_continuous(expand=c(0,0)) +
  theme(plot.margin = unit(c(1,0,0,0), "lines")) +
  theme(axis.text.y = element_text(color = c(NA, "black", NA, "black", NA,  "black", NA,  "black", NA, "black"))) +
  theme(axis.text.y = element_blank(),
        axis.title.y = element_blank()) +
  ggtitle("") +
  annotation_custom(myGrob4, xmin=1998, xmax=2060, ymin=35, ymax=70) 


plot.COIL.p2.g <- ggplotGrob(plot.COIL.p2)
plot.COIL.p2.g$layout$clip[plot.COIL.p2.g$layout$name=="panel"] <- "off"



myGrob4 <- grobTree(textGrob("Unclassified Sphaeropleales", rot=60))

plot.COIL.p3 <- ggplot()+
  geom_line(data=strat.per.red.depth.melt.p3, aes(Age, value), colour="#c94fb1") +
  geom_area(data=strat.per.red.depth.melt.p3, aes(Age, value), fill="#c94fb1")+ 
  geom_bar(data=strat.per.red.depth.melt.p3, aes(Age, value), width=0.1, stat="identity")+
  xlab("")+ylab("") +
  coord_flip() +
  theme_new  +
  scale_x_continuous(expand=c(0,1), limits = c(1760, 1990), breaks = c(1775, 1800, 1825, 1850, 1875, 1900, 1925, 1950, 1975, 2000)) + 
  scale_y_continuous(expand=c(0,0), breaks = c(0, 20, 40)) +
  theme(plot.margin = unit(c(1,0,0,0), "lines")) +
  theme(axis.text.y = element_text(color = c(NA, "black", NA, "black", NA,  "black", NA,  "black", NA, "black"))) +
  theme(axis.text.y = element_blank(),
        axis.title.y = element_blank()) +
  ggtitle("") +
  annotation_custom(myGrob4, xmin=2005, xmax=2060, ymin=60, ymax=90) 


plot.COIL.p3.g <- ggplotGrob(plot.COIL.p3)
plot.COIL.p3.g$layout$clip[plot.COIL.p3.g$layout$name=="panel"] <- "off"


myGrob4 <- grobTree(textGrob("Nannochloropsis sp.", rot=60))


plot.COIL.p4 <- ggplot()+
  geom_line(data=strat.per.red.depth.melt.p4, aes(Age, value), colour="#c94fb1") +
  geom_area(data=strat.per.red.depth.melt.p4, aes(Age, value), fill="#c94fb1")+ 
  geom_bar(data=strat.per.red.depth.melt.p4, aes(Age, value), width=0.1, stat="identity")+
  xlab("Age")+ylab("") +
  coord_flip() +
  theme_new  +
  scale_x_continuous(expand=c(0,1), limits = c(1760, 1990), breaks = c(1775, 1800, 1825, 1850, 1875, 1900, 1925, 1950, 1975, 2000)) + 
  scale_y_continuous(expand=c(0,0), breaks = c(0, 20, 40)) +
  theme(plot.margin = unit(c(1,0,0,0), "lines")) +
  ggtitle("") + theme(axis.text.y = element_blank(),
                      axis.title.y = element_blank(),
                      axis.title.x = element_text(hjust=1.15)) +
  annotation_custom(myGrob4, xmin=1998, xmax=2050, ymin=35, ymax=70) +
  annotation_custom(myGrob5, xmin=zn.df$zone[1], xmax=zn.df$zone[1], ymin=-450, ymax=170) +
  annotation_custom(myGrob3, xmin=zn.df$zone[2], xmax=zn.df$zone[2], ymin=-450, ymax=170)

plot.COIL.p4.g <- ggplotGrob(plot.COIL.p4)
plot.COIL.p4.g$layout$clip[plot.COIL.p4.g$layout$name=="panel"] <- "off"




################## INPUT HERE THE ORDERED SUBSET GROUPS ########################

downcore.ps = prune_samples(sample_sums(AlgaeShortOrder)>0, AlgaeShortOrder) %>%
  filter_taxa(function(x) sum(x) > 0, TRUE)

downcore.ps.perc <- transform_sample_counts(downcore.ps, function(OTU) OTU/sum(OTU)*100)

downcore.ps.perc.df <- data.frame(otu_table(downcore.ps.perc))

## TAXONOMY to rename OTU labels
downcore.ps.tax.df <- data.frame(tax_table(downcore.ps))
rownames(downcore.ps.tax.df) <- downcore.ps.tax.df$Genus

#Transpose OTU table
downcore.ps.perc.df<-as.data.frame(t(downcore.ps.perc.df))
# make a column of OTU table the taxonomy
downcore.ps.perc.df$otulabels <- downcore.ps.tax.df$Genus
rownames(downcore.ps.perc.df) <- downcore.ps.perc.df$otulabels

downcore.ps.perc.df$`otulabels`=NULL
# transpose back
downcore.ps.perc.df<-as.data.frame(t(downcore.ps.perc.df))

# metadata
downcore.ps.perc.metadata <- data.frame(sample_data(downcore.ps))

rownames(downcore.ps.perc.df) <- downcore.ps.perc.metadata$Age


strat.clust <- downcore.ps.perc.df %>% 
  vegdist(method = "bray", binary = FALSE, diag = FALSE, upper = FALSE,
          na.rm = FALSE) %>% 
  chclust(method = "coniss")

#bstick(strat.clust, 20)
downcore.ps.perc.metadata$Age <- as.numeric(as.character(downcore.ps.perc.metadata$Age))
ct <- cutree(strat.clust, 3)
zn <- which(diff(ct) > 0)
zone <- (downcore.ps.perc.metadata$Age[zn] + downcore.ps.perc.metadata$Age[zn + 1])/2
zn.df <- data.frame(zone)

dendro <- as.dendrogram(strat.clust)
ddata <- dendro_data(dendro, type="rectangle")

depth <- downcore.ps.perc.metadata$Age


ddata$segments->yo
yo$xx=NA
yo$xxend=NA

for(i in 1:length(depth)){
  yo$xx[which(yo$x == i)]=depth[i]
  yo$xxend[which(yo$xend == i)]=depth[i]
  da=depth[i+1]-depth[i]
  yo$xx[yo$x>i & yo$x<i+1]<-depth[i]+(yo$x[yo$x>i & yo$x<i+1]-i)*da
  yo$xxend[yo$xend>i & yo$xend<i+1]<-depth[i]+(yo$xend[yo$xend>i & yo$xend<i+1]-i)*da
} # Please dont ask....
ddata$segments$x<-yo$xx
ddata$segments$xend<-yo$xxend


dendro.COIS = ggplot(segment(ddata)) +
  geom_segment(aes(x=x, y=y, xend=xend, yend=yend)) +
  coord_flip() +
  ylab("Distance")+
  #scale_x_reverse()+
  theme_new + theme(axis.title.y=element_blank(),
                    axis.text.y=element_blank(),
                    axis.ticks.y=element_blank(),
                    axis.line.y = element_blank())   +
  scale_x_continuous(expand=c(0,1), limits = c(1760, 1990)) + 
  scale_y_continuous(expand=c(0,0), limits = c(0, 15)) 


strat.per.red <- downcore.ps.perc.df %>% 
  chooseTaxa(n.occ = 5)


strat.per.red1 <- strat.per.red %>% select("Unclassified_Phylum_Chlorophyta_2", "Mychonastes", "Unclassified_Order_Sphaeropleales_1")

strat.per.red.depth <- rownames_to_column(strat.per.red1) 

colnames(strat.per.red.depth)[1] <- "Age"

strat.per.red.depth <- strat.per.red.depth %>% arrange(Age) 

strat.per.red.depth.melt <- reshape2::melt(strat.per.red.depth, id.vars=c("Age"))

strat.per.red.depth.melt$Age <- as.numeric(strat.per.red.depth.melt$Age)


strat.per.red.depth.melt.p1 <-  strat.per.red.depth.melt %>% filter(variable == "Unclassified_Phylum_Chlorophyta_2")
strat.per.red.depth.melt.p1$variable <- gsub("Unclassified_Phylum_Chlorophyta_2", "Unclassified Chlorophyta", strat.per.red.depth.melt.p1$variable)

strat.per.red.depth.melt.p2 <-  strat.per.red.depth.melt %>% filter(variable == "Mychonastes")
strat.per.red.depth.melt.p2$variable <- gsub("Mychonastes", "Mychonastes sp.", strat.per.red.depth.melt.p2$variable)


strat.per.red.depth.melt.p2 <-  strat.per.red.depth.melt %>% filter(variable == "Unclassified_Order_Sphaeropleales_1")
strat.per.red.depth.melt.p2$variable <- gsub("Unclassified_Order_Sphaeropleales_1", "Unclassified Sphaeropleales", strat.per.red.depth.melt.p2$variable)


myGrob4 <- grobTree(textGrob("Unclassified Chlorophyta", rot=60))


plot.COIS.p1 <- ggplot()+
  geom_line(data=strat.per.red.depth.melt.p1, aes(Age, value), colour="#FFE25C") +
  geom_area(data=strat.per.red.depth.melt.p1, aes(Age, value), fill="#FFE25C")+ 
  geom_bar(data=strat.per.red.depth.melt.p1, aes(Age, value), width=0.1, stat="identity")+
  xlab("")+ylab("") +
  coord_flip() +
  theme_new  +
  scale_x_continuous(expand=c(0,1), limits = c(1760, 1990), breaks = c(1775, 1800, 1825, 1850, 1875, 1900, 1925, 1950, 1975, 2000)) + 
  scale_y_continuous(expand=c(0,0), breaks = c(0, 25, 50, 75)) +
  theme(plot.margin = unit(c(1,0,0,0), "lines")) +
  theme(axis.text.y = element_text(color = c(NA, "black", NA, "black", NA,  "black", NA,  "black", NA, "black"))) +
  ggtitle("c) AlgaeCOI-S") +
  annotation_custom(myGrob4, xmin=2006, xmax=2050, ymin=22, ymax=70) 


plot.COIS.p1.g <- ggplotGrob(plot.COIS.p1)
plot.COIS.p1.g$layout$clip[plot.COIS.p1.g$layout$name=="panel"] <- "off"


myGrob4 <- grobTree(textGrob("Mychonastes sp.", rot=60))


plot.COIS.p2 <- ggplot()+
  geom_line(data=strat.per.red.depth.melt.p2, aes(Age, value), colour="#FFE25C") +
  geom_area(data=strat.per.red.depth.melt.p2, aes(Age, value), fill="#FFE25C")+ 
  geom_bar(data=strat.per.red.depth.melt.p2, aes(Age, value), width=0.1, stat="identity")+
  xlab("")+ylab("Rel. Abundance (%)") +
  coord_flip() +
  theme_new  +
  scale_x_continuous(expand=c(0,1), limits = c(1760, 1990), breaks = c(1775, 1800, 1825, 1850, 1875, 1900, 1925, 1950, 1975, 2000)) + 
  scale_y_continuous(expand=c(0,0), breaks = c(0,15)) +
  theme(plot.margin = unit(c(1,0,0,0), "lines")) +
  theme(axis.text.y = element_text(color = c(NA, "black", NA, "black", NA,  "black", NA,  "black", NA, "black"))) +
  theme(axis.text.y = element_blank(),
        axis.title.y = element_blank()) +
  ggtitle("") +
  annotation_custom(myGrob4, xmin=2005, xmax=2030, ymin=25, ymax=45) 


plot.COIS.p2.g <- ggplotGrob(plot.COIS.p2)
plot.COIS.p2.g$layout$clip[plot.COIS.p2.g$layout$name=="panel"] <- "off"




myGrob4 <- grobTree(textGrob("Unclassified Sphaeropleales", rot=60))


plot.COIS.p3 <- ggplot()+
  geom_line(data=strat.per.red.depth.melt.p3, aes(Age, value), colour="#FFE25C") +
  geom_area(data=strat.per.red.depth.melt.p3, aes(Age, value), fill="#FFE25C")+ 
  geom_bar(data=strat.per.red.depth.melt.p3, aes(Age, value), width=0.1, stat="identity")+
  xlab("Age")+ ylab("") +
  coord_flip() +
  theme_new  +
  scale_x_continuous(expand=c(0,1), limits = c(1760, 1990), breaks = c(1775, 1800, 1825, 1850, 1875, 1900, 1925, 1950, 1975, 2000)) + 
  scale_y_continuous(expand=c(0,0), breaks = c(0, 25, 50)) +
  theme(plot.margin = unit(c(1,0,0,0), "lines")) +
  ggtitle("") + theme(axis.text.y = element_blank(),
                      axis.title.y = element_blank(),
                      axis.title.x = element_text(hjust=0.7)) +
  annotation_custom(myGrob4, xmin=2015, xmax=2050, ymin=18, ymax=70) +
  annotation_custom(myGrob2, xmin=1906, xmax=1906, ymin=-1006, ymax=150) +
  annotation_custom(myGrob3, xmin=zn.df$zone[1], xmax=zn.df$zone[1], ymin=-162, ymax=115) +
  annotation_custom(myGrob5, xmin=zn.df$zone[2], xmax=zn.df$zone[2], ymin=-162, ymax=115) +
  annotation_custom(myGrob, xmin=1906, xmax=1870, ymin=-1006, ymax=700)

plot.COIS.p3.g <- ggplotGrob(plot.COIS.p3)
plot.COIS.p3.g$layout$clip[plot.COIS.p3.g$layout$name=="panel"] <- "off"




coniss.plot <- plot_grid(plot.18S.p1.g, NULL, plot.18S.p2.g, NULL, plot.18S.p3.g, NULL, plot.18S.p4.g, NULL, dendro.18S, NULL, 
                         plot.COIL.p1.g, NULL, plot.COIL.p2.g, NULL, plot.COIL.p3.g, NULL, plot.COIL.p4.g, NULL, dendro.COIL, NULL, 
                         plot.COIS.p1.g, NULL, plot.COIS.p2.g, NULL, plot.COIS.p3.g, NULL, dendro.COIS, align = "hv", axis="tblr", ncol=27, 
                         rel_widths = c(0.625, -0.1, 0.7, -0.1, 0.4, -0.1, 0.35, -0.1, 0.47, -0.1, 
                                        0.625, -0.1, 0.5, -0.1, 0.4, -0.1, 0.4, -0.1, 0.47, -0.1,
                                        0.625, -0.1, 0.32, -0.1, 0.5, -0.1, 0.47))


ggsave("coniss.plot.algae.jpeg", width=15, height=9, units = "in", bg = "transparent")
coniss.plot
dev.off()



##########################################################################################
################# Supplementary plots #################################################
###########################################################################################


downcore.ps_final <- prcurve(downcore.ps.perc.df_hel, method = "ca", plotit = FALSE, vary =TRUE, maxit = 50)
#using CA as DCA axis <2 (warrants unimodal method)
#downcore.ps_final

Bact_scores <- scores(downcore.ps_final, display = "curve")
Bact_scores <- as.data.frame(Bact_scores)
Bact_scores$Age <- rownames(Bact_scores)
Bact_scores$SampleInterval <- downcore.ps.perc.metadata$SampleInterval
Bact_scores %<>%
  mutate(var = "prc")%>%
  mutate(Age = as.numeric(Age))

Bact_scores$Era <- "Trout"
Bact_scores$Era[Bact_scores$Age <= 1900] <- "NoTrout"

PRC.plot <- ggplot()+
  geom_point(data=Bact_scores, aes(y=PrC, x=Age, colour=Era), size=3)+
  geom_line(data=Bact_scores, aes(y=PrC, x=Age))+
  scale_color_manual(values=c("Trout" = "#189100","NoTrout" = "#0066dc")) +
  labs(y="PrC scores")+
  coord_flip()+
  #scale_x_reverse()+
  theme_classic()+
  # ylim(0, 5) +
  theme(legend.position = "")
PRC.plot

########################### GAM PLOT ########################


k.length <-  nrow(Bact_scores)
#K check

prc_gamm_check <- gam(PrC ~s(Age, k=3), data = Bact_scores, method = "REML") 
#k needs to be adjusted so that the k-index is > 1 (see what prints out from the following line)

gam.check(prc_gamm_check)
rsd <- residuals(prc_gamm_check)
plot(Bact_scores$Age,rsd);

qq.gam(prc_gamm_check,rep=100)
#Once you find the appropriate k, run the actual gamm while adding the depth (or age) correlation
gam_bact <-gamm(PrC ~s(Age, k = 3, bs = "tp", m = 1), data = Bact_scores, correlation = corCAR1(form = ~ Age))

gam_bact_d<- Deriv(gam_bact, n=200)

gam_bact_dat <- with(Bact_scores,
                     data.frame(Age = seq(min(Age), max(Age),
                                          length = 200)))

bact_prc_p1 <- predict(gam_bact$gam, newdata = gam_bact_dat)
bact_prc_CI <- confint(gam_bact_d, alpha = 0.01)
bact_prc_S = signifD(bact_prc_p1, gam_bact_d$Age$deriv, bact_prc_CI$Age$upper, bact_prc_CI$Age$lower, eval = 0)


gam_bact_dat <- gam_bact_dat %>% mutate(gam = bact_prc_p1, incr = bact_prc_S$incr, decr = bact_prc_S$decr, upper = bact_prc_CI$epth$upper, lower = bact_prc_CI$Age$lower)

# trial GAM plot
gam.plot <- ggplot()+
  geom_point(data = Bact_scores, aes(y = PrC, x= Age, color = Era))+
  geom_path(data = gam_bact_dat, aes(x=Age, y=gam))+
  geom_path(data = gam_bact_dat, aes(x=Age, y=incr), colour = "red", size = 1)+
  geom_path(data = gam_bact_dat, aes(x=Age, y=decr), colour = "blue", size = 1)+
  #geom_hline(yintercept=0)+
  #scale_x_reverse()+
  #scale_y_reverse()+
  coord_flip()+
  scale_color_manual(values=c("Trout" = "#189100", "NoTrout" = "#0066dc")) +
  labs(y = "Prc Score", x = "Age") + #, title = paste(downcore.list[i], #sep = ".")) +
  theme_classic()+
  theme(legend.position = "")

gam.plot


######################### Plot the final GAM ################################

# As geom_vline, add in the CONISS shift date/s manually

ggplot(data = Bact_scores,
       aes(x = Age, y = PrC)) +
  
  geom_line(data = gam_bact_dat,
            aes(x = Age, y = fit),
            inherit.aes = FALSE,
            size = 1, show.legend = FALSE) +
  geom_line(data = gam_bact_dat,
            aes(x = Age, y = incr),
            inherit.aes = FALSE,
            colour = "tomato1", size = 1) +
  geom_line(data = gam_bact_dat,
            aes(x = Age, y = decr),
            inherit.aes = FALSE,
            colour = "steelblue1", size = 1) +
  geom_ribbon(data = gam_bact_dat,
              aes(x = Age, y = NULL,
                  ymin = lower, ymax = upper),
              alpha = 0.3, 
              inherit.aes = FALSE, fill = "grey") +
  geom_vline(xintercept = c(1929, 1969), linetype = "dashed",
             colour = "lightgrey") +
  geom_point(aes(fill = Era), size = 2, colour = "white", shape = 21,
             show.legend = FALSE) +
  scale_x_continuous(limits = c(1775, 2000), breaks = seq(1775, 2000, 25),
                     expand = c(0, 0)) +
  scale_y_continuous(limits = c(0, 25), breaks = seq(0,25, 5),
                     expand = c(0, 0)) +
  labs(x = "Age (CE)", y = "PrC") +
  theme_classic()

#600x400


####################################################################################################
############################# BAR PLOTS FOR SUPPLEMENTARY ##########################################
####################################################################################################

y12 <- tax_glom(AlgaeShort, taxrank = 'Class') # agglomerate taxa
y13 <- transform_sample_counts(y12, function(x) x/sum(x)) #get abundance in %
y14 <- psmelt(y13) # create dataframe from phyloseq object
y14$Class <- as.character(y14$Class) #convert to character
y14$Class[y14$Abundance < 0.01] <- "Species < 1% abund." #rename genera with < 1% abundance

p1 = ggplot(data=y14, aes(x=Sample, y=Abundance, fill=Class)) 
p1$data$Sample <- factor(p1$data$Age)

p1 + geom_bar(aes(), stat="identity",position="Stack") + 
  xlab("Median age") + ylab("Proportional abundance") + coord_flip() + 
  #guides(fill=guide_legend(ncol=2)) +
  scale_fill_manual(values = scale38)  + 
  theme(panel.background = element_rect(fill="transparent"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.background = element_rect(fill="transparent", colour=NA),
        axis.title = element_text(size = 12, color = "black",
                                  face = "italic"),
        legend.title = element_text(face = "italic", color="black")) +
  guides(color = guide_legend(override.aes = list(size = 0.3)))
