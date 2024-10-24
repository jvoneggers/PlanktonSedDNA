## 8. Create Phyloseq
require(phyloseq)
require(tidyverse)

### a. read in metadata tables
metadata<-read.csv(list.files(pattern = "*_sample_data_dates.csv", full.names = TRUE), row.names=1, header=T)

### b. read in V7 files
rarefied_esv_table <- read.csv("1_LoadData/2024-06-19_18S_V7_ZOTUexact_tax80_Rarefied_ESV_table.csv", row.names=1, header=T)
taxonomy <- read.csv("1_LoadData/2024-06-19_18S_V7_ZOTUexact_tax80_tax_table_notnorm.csv", row.names=1, header=T)

### c. Put into Phyloseq
#put in "Unassigned" in blanks in taxa tab
tax_tab<-as.data.frame(taxonomy)
tax_tab[is.na(tax_tab)==TRUE]<-"Unassigned"
tax_tab <-as.matrix(tax_tab)

#covert taxa table to phyloseq sub-object
tax_tab <- tax_table(tax_tab)
rm(taxonomy)

#convert metadata into a phyloseq object
row.names(metadata)<-metadata$sample_id
samp_dat <- sample_data(metadata)

#convert ESV table into a phyloseq object
esv_tab_norm <- otu_table(rarefied_esv_table, taxa_are_rows = T)
rm(rarefied_esv_table)

# make phyloseq objects
ps_V7_norm <- phyloseq(esv_tab_norm, samp_dat, tax_tab)

#remove phyloseq sub-objects (not metadata because we want that separately)
rm(esv_tab_norm)
rm(tax_tab)
rm(samp_dat)

#transform sample data
#ps_V7_norm_tr <- transform_sample_counts(ps_V7_norm, function(x) x / sum(x))


### d. read in V9 files
rarefied_esv_table <- read.csv("1_LoadData/2024-06-19_18S_V9_ZOTUexact_tax80_Rarefied_ESV_table.csv", row.names=1, header=T)
taxonomy <- read.csv("1_LoadData/2024-06-19_18S_V9_ZOTUexact_tax80_tax_table_notnorm.csv", row.names=1, header=T)

### e. put into Phyloseq
#put in "Unassigned" in blanks in taxa tab
tax_tab<-as.data.frame(taxonomy)
tax_tab[is.na(tax_tab)==TRUE]<-"Unassigned"
tax_tab <-as.matrix(tax_tab)

#covert taxa table to phyloseq sub-object
tax_tab <- tax_table(tax_tab)
rm(taxonomy)

#convert metadata into a phyloseq object
row.names(metadata)<-metadata$sample_id
samp_dat <- sample_data(metadata)

#convert ESV table into a phyloseq object
esv_tab_norm <- otu_table(rarefied_esv_table, taxa_are_rows = T)
rm(rarefied_esv_table)

# make phyloseq objects
ps_V9_norm <- phyloseq(esv_tab_norm, samp_dat, tax_tab)

#remove phyloseq sub-objects (not metadata because we want that separately)
rm(esv_tab_norm)
rm(tax_tab)
rm(samp_dat)

#transform sample data
#ps_V9_norm_tr <- transform_sample_counts(ps_V9_norm, function(x) x / sum(x))