###############################################################################################
### This code takes the output .clstr files from cd-hit-est and assigns taxonomic classification 
### to sequences
### It also converts data to phyloseq object
###############################################################################################

##Install packages
#install.packages("rlist")
#install.packages("taxize")
#install.packages("TAI")
#install.packages("gtools")
#install.packages("rstatix")
#install.packages("rio")
#install.packages("filesstrings")
#install.packages("stringr")
#install.packages("dplyr")
#install.packages("tidyr")
#install.packages("stringi")
#install.packages("ggplot2")
#install.packages("data.table")
#install.packages("phia")
#install.packages("multcomp")
#install.packages("base")

##Load required packages
library(filesstrings)
library(stringr)
library(dplyr)
library(tidyr)
library(stringi)
library(rlist)
library(taxize)
library(myTAI)
library(gtools)
library(ggplot2)
library(data.table)
library(phia)
library(multcomp)
library(rstatix)
library(base)
library(rio)

## Ensure all .clstr files are in the same folder and this is the working directory
all.files <- list.files(pattern = "*.clstr")

## The .clstr file is first changed into a dataframe that is easier to manipulate.
## Then upstream taxonomy is assigned to sequences in the .clstr file using NCBI so that we can see 
## whether the sample reads have been mapped to the correct species, given not all genes provide species level resolution. 
## If clustering occurs with the same species for a given gene region, then the mapping was correct at the species level, otherwise the next common 
## taxonmic ranking is recorded and assigned to the sample sequence read.

##Change .clstr file to dataframe
for (i in all.files){
  sample <- tools::file_path_sans_ext(i)
  #Load clstr file output from cd-hit
  clstr <- read.csv(i, sep = "\t", row.names = NULL, header = FALSE, stringsAsFactors = FALSE)
  head(clstr)
  #Format into a dataframe
  clstr2 <- clstr
  n = nrow(clstr)
  x = 0
  numbers_only <- function(x) !grepl("\\D", x)
  for (row in c(1:n)) {
    if (numbers_only(clstr2[row,1]) == TRUE) {
      clstr2[row,1] <- x}
    else (NULL)
    x <- clstr2[row,1]
  }
  head(clstr2)
  clstr.sums <- data.frame(dplyr::count(clstr2,V1))
  head(clstr.sums)
  switch <- clstr.sums[1,2]
  clstr3 <- cbind(clstr2[1], clstr)
  head(clstr3)
  clstr4 <- clstr2[-which(clstr2$V2 == ""), ]
  head(clstr4)
  clstr5 <- clstr4
  clstr5[] <- lapply(clstr5, gsub, pattern='>', replacement='')
  clstr5.2 <- data.frame(str_split_fixed(clstr5$V2, "nt, ", 2))
  clstr5.3 <- data.frame(str_split_fixed(clstr5.2$X2, "... ", 2))
  clstr6 <- cbind(clstr5[1],clstr5.2[1],clstr5.3[1:2])
  head(clstr6)
  clstr7 <- data.frame(str_split_fixed(clstr6$V1, " ", 2))
  clstr8 <- cbind(clstr6[-1], clstr7[1:2])
  
  #check to see it has worked- should be a lot easier to read now
  head(clstr8)
  View(clstr8)
  
  #rename columns
  colnames(clstr8) <- c("nt","gene","stat","cluster", "cluster_no")
  
  #create a column that just has the % similarity value
  numextract <- function(string){ 
    str_extract(string, "\\-*\\d+\\.*\\d*")
  } 
  clstr8$perc_cluster <- numextract(clstr8$stat)
  
  #check to see it has worked
  head(clstr8)
  
  #Now all columns are separated so the data can be explored properly
  View(clstr8)
  
  #subset by the chloroplast genes that were targeted in the hybridization capture
  maingenes <- c("accD", "atpB", "atpF", "atpH", "atpI", "matK", "ndhC", "ndhF", "ndhK", "petA", "petD", "psbA", "psbD", "psbE","psbH", "psbK","psbZ", "rbcL", "rpl16", "rpoC1")
  
  #search for main genes within the names of the sequences and export gene match to new column
  clstr8$Gene <- (str_extract(clstr8$gene, paste(maingenes, collapse = "|")))
  
  #keep only samples for 20 genes of interest
  clstr8sub <-  clstr8 %>% filter(str_detect(clstr8$Gene, paste(maingenes, collapse = "|")))
  
  #subset by sample reads
  clstr9 <- filter(clstr8sub, grepl('Soil', gene))
  View(clstr8sub)
  
  #only retain clusters that contain your sample reads
  final <- clstr8sub %>% filter(clstr8sub$cluster_no %in% clstr9$cluster_no)
  
  #create a new column to denote which rows are your samples
  final2 <- final %>% mutate(sample = case_when(str_detect(final$gene, regex("Soil")) ~ "soil sample"))
  View(final2)
  
  #get rid of sample name prefix so that species can be looked at
  patx='Soil_'
  final2$gene <- str_remove(final2$gene, paste(patx))

  #extract sample species name to new column
  pat= '[A-Z][a-z]*_[a-z]*'
  final2$species <- str_extract(final2$gene, paste(pat))
  
  #check everything has worked correctly thus far
  View(final2)
  #keeping only columns of interest
  df <- final2[c(2,5,6,7,8,9)]
  head(df)
  View(df)
  
  #ensure representative cluster sequences are indicated (these are the longest sequences that cd-hit identifies then aligns other sequences to)
  df$perc_cluster[is.na(df$perc_cluster)] <- "Representative_cluster_sequence"
  
  #getting species column to not have an underscore for better matching
  df$species <- gsub('_', ' ', df$species)
  
  # getting rid of samples that only clustered with themselves
  as.numeric(df$cluster_no)
  
  df2 <- df %>% distinct(species,gene,cluster_no, .keep_all= TRUE)
  #View(df2)
  df3 <- df2 %>% group_by(cluster_no) %>% filter(n()>1)
  View(df3)
  
  
  ##output taxonomic classifications from species column using taxize and ncbi database
  if (nrow(df3)>0){
    taxclass <- as.data.frame(df3$species)
    outputlist <- apply(taxclass, 1,  function(x) taxonomy(organism = x , db = "ncbi", output = "classification" ))
    taxdata = data.frame()
    for(x in 1:length(outputlist)){
      tryCatch({
        phylum=filter(outputlist[[x]], rank =="phylum")$name
        class=filter(outputlist[[x]], rank =="class")$name
        order=filter(outputlist[[x]], rank =="order")$name
        family=filter(outputlist[[x]], rank =="family")$name
        genus=filter(outputlist[[x]], rank =="genus")$name
        species=filter(outputlist[[x]], rank =="species")$name
        row <- data.frame(cbind(phylum=phylum,class=class,order=order,family=family,genus=genus,species=species))
        taxdata <- bind_rows(taxdata, row)    
      }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
    }}
  View(taxdata)
  #merge original dataframe with taxonomy dataframe
  m1 <- merge(df3, taxdata, by = "species", all = TRUE)
  View(m1)
  
  #complete taxonomy data frame by adding a species column
  m1$species1 <- gsub(".* ","", m1$species)
  
  #getrid of duplicates
  m2 <- m1 %>% distinct(species,gene,cluster_no, .keep_all= TRUE)
  View(m2)

  #save and export results so far
  write.csv(m2, file = paste({sample}, ".csv", sep = "_all"))
  
}

#All .clstr files have now been processed and upstream taxonomy assigned
#The following code assigns taxonomy to the species that were not found in the NCBI database either through not being deposited or incorrect naming

file <- read.csv("NA_samples.csv") # this file is a list of the species where upstream taxonomy could not be found and can be updated as needed 

all.files <- list.files(pattern = "*_all.csv")

for (i in all.files){
  sample <- tools::file_path_sans_ext(i)
  m2 <- read.csv(i)
  m3 <- filter(m2[is.na(m2$phylum)==0, ])
  no <- filter(m2,is.na(phylum))
  no$phylum <- file$phylum[match(no$species,file$species)]
  no$class<- file$class[match(no$species,file$species)]
  no$order <- file$order[match(no$species,file$species)]
  no$family <- file$family[match(no$species,file$species)]
  no$genus <- file$genus[match(no$species,file$species)]
  no2 <- rbind(no,m3)
  write.csv(no2, paste0(sample, "_noNA.csv")) 
}

#Now that upstream taxonomy has been assigned to every species in the cluster the next step is to take your samples and 
#classify the highest taxonomic ranking the samples classify at

all.files <- list.files(pattern = "*_all_noNA.csv")
for (i in all.files){
  sample <- tools::file_path_sans_ext(i)
  m2 <- read.csv(i)
  
  #output another column denoting lowest common taxonomic level for clusters
  #n_distinct == 1 is when they are all the same whereas != 1 is when they are different
  
  m3 <- m2 %>% group_by(cluster_no) %>%
    mutate(rank = case_when(n_distinct(class, na.rm = TRUE) != 1  ~ as.character(phylum), n_distinct(order, na.rm = TRUE) != 1  ~ as.character(class), n_distinct(family, na.rm = TRUE) != 1  ~ as.character(order), n_distinct(genus, na.rm = TRUE) != 1 ~ as.character(family), 
                            n_distinct(species, na.rm = TRUE) != 1 ~ as.character(genus), n_distinct(species, na.rm = TRUE) == 1 ~ as.character(species))) 
  View(m3)
  
  
  # retain only samples from the main dataframe
  m4 <- subset(m3, m3$sample == "soil sample") 
  
  final <- m4
  
  #check final table
  View(final)
  
  #save and expoert results to csv file
  write.csv(final, file = paste({sample}, ".csv", sep = "_final"))
  
}

#The final data files are then assigned a depth column (for where they were collected in the core) and combined into one single file combining all replicates
#This part of the code will change depending how the files are named

#install required packages
#install.packages("gridExtra")
#install.packages("cowplot")

#load required packages
library(ggplot2)
library(gridExtra)
library(cowplot)

##clear global environment
rm(list = ls())

#collate all final output.csv files from cd-hit-est analysis
list<- list.files(pattern = ".soilref_all_noNA_final.csv") #adapt to file name

## this is to create a depth column for all the samples (sample filenames are their depth)
for (i in list){
  d = read.csv(i)
  dep = gsub("output_TIN3_", "", i) # adapt to sample name
  depth = gsub("*[a-z].soilref_all_noNA_final.csv", "", dep)
  d$depth <- print(depth)
  filename = gsub(".soilref_all_noNA_final.csv", "", i)
  assign(filename,d)   
}

#combine all files for each sediment core
allfiles <- do.call(rbind, mget(ls(pattern="[0-9].*"))) 

as.data.frame(allfiles)
row.names(allfiles) <- NULL
View(allfiles)

##specify name of file here for all samples in single core
write.csv(allfiles, "allfilesTIN3.csv")


##Assigning a new column called "level" to assign Class, Order, Family Genus and species to the ranking column

#load packages
library(dplyr)
library(tidyverse)
library(RColorBrewer)

samples <- list.files(pattern = "allfilesTIN3.csv") #change to name of own file

for (i in samples){
  pat = ('.csv')
  name <- str_remove(i, pat)
  filename <- paste(name)
  file <- read.csv(i, stringsAsFactors=F)
  file <- na.omit(file)#get rid of any lingering NA's
  head(file)
  
  ##Get rid of rank classification to "Streptophyta" as this is the phylum, which will be the same for all plants
  file <- subset(file, rank != "Streptophyta")

  file1 <- file %>% mutate(level = ifelse(str_detect(file$rank, file$class), "class", "")

  file2 <- subset(file1, level != "class")
  
  file2$level = ifelse(str_detect(file2$rank, file2$order), "order", "")
  
  file3 <- subset(file2, level != "order")
  
  file3$level = ifelse(str_detect(file3$rank, file3$family), "family", "")
  
  file4 <- subset(file3, level != "family")
  
  file4$level = ifelse(str_detect(file4$rank, file4$species1), "species", "")
  
  file5 <- subset(file4, level != "species")
  
  file5$level = ifelse(str_detect(file5$rank, file5$genus), "genus", "")
  
  m1 <- subset(file1, level == "class")
  m2 <- subset(file2, level == "order")
  m3 <- subset(file3, level == "family")
  m4 <- subset(file4, level == "species")
  m5 <- subset(file5, level == "genus")
  
  fileX <- rbind(m1,m2,m3,m4,m5)
  View(fileX)
  

  write.csv(fileX, paste0(filename, "_withlevel.csv")) #overwrites original file
}

#-----------------------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------------------
#-----------------------------The data is then ready to be analyzed and plotted------------------------------------------

#---------------------------------------------------------------------------------------------------------
############################### Converting files to phyloseq ###########################
#load phyloseq package
library("phyloseq")

#load all files from previous work - change names to soil cores analysed
data.TIN1 <- read.csv("allfilesTIN1_withlevel.csv")
data.TIN2<- read.csv("allfilesTIN2_withlevel.csv")
data.TI2 <- read.csv("allfilesTI2_withlevel.csv")
data.TIN3 <- read.csv("allfilesTIN3_withlevel.csv")

#keeping only species hits that were recovered for more than 4 genes - based on Foster et al. 2021 https://doi.org/10.3389/fevo.2021.735744
data2 <- data.TIN1 %>% group_by(depth, species) %>% distinct(Gene) %>% filter(n()>4) #change for each core
View(data2)

data <- merge(data.TIN1, data2, by = c("species", "depth", "Gene"), all = FALSE) #change for each core
#View(data)
head(data)

#convert to data frame 
play <- as.data.frame(data)[c("rank", "level")]

#get unique taxa within each depth/sample to give a number to - referring to the taxa column as the unique taxa i.e. what classifications did we get for each samples
play2 <- play %>% distinct(rank, .keep_all= TRUE)

#generate  numbers for unique taxa (using OTU for compatibility with phyloseq)
rownames(play2) <- paste0("OTU", 1:nrow(play2))

#get rownames as OTU numbers (not really OTU - just for compatibility with phyloseq)
play2$otu <- rownames(play2)

#merge OTU data frame and original retaining the OTU numbers 
m1 <- merge(data, play2, by = c("rank", "level"), all = TRUE)
#View(m1) #check all is okay
head(m1)

#creating otu data i.e. just OTU numbers and sample and preparing it for physeq conversion
otu.data <- m1[c("otu", "depth")] # subset data
otu.data2 <- otu.data %>% group_by(depth, otu) %>% count() 
otu.data2 <- as.data.frame(otu.data2)
otu.data2_tryagain <- pivot_wider(otu.data2, names_from = depth, values_from = n) # get into correct table format
View(otu.data2_tryagain)

otu.data2_tryagain <- as.data.frame(otu.data2_tryagain)

otu.data2_tryagain[is.na(otu.data2_tryagain)] <- 0 # change NA to zero so that phyloseq functions work

rownames(otu.data2_tryagain) <- otu.data2_tryagain$otu #rownames as OTU for matching to TAX file

FINAL_sam_data <- otu.data2_tryagain

FINAL_sam_data <- as.data.frame(otu.data2_tryagain)[,-1] # remove OTU column as is row rownames

View(FINAL_sam_data)

FINAL_sam_data <- as.matrix(FINAL_sam_data)

#preparing tax_table
FINAL_tax_table <-play2[-3]

FINAL_tax_table <- as.matrix(FINAL_tax_table)

# create OTU physeq object
OTU = otu_table(FINAL_sam_data, taxa_are_rows = TRUE)

#create TAX phyloseq object
TAX = tax_table(FINAL_tax_table)

#check both
OTU
TAX

#create final phyloseq file
TIN1_phy <- phyloseq(OTU,TAX)

##Creating sample data for core
##Create a .csv file with column headers sampleID, Depth (fill data)
#reading in sample data 
sampledata <- read.csv("sample_data_TIN1.csv")
sampledata$Depth <- as.factor(sampledata$Depth)

rownames(sampledata) <- sampledata$sampleID

samplesdata <- sampledata[, -c(1)]

SAMP <- as.data.frame(samplesdata)

samp_TIN1 <- sample_data(SAMP)

# final phyloseq file
physeq_final_TIN1 <- merge_phyloseq(TIN1_phy, samp_TIN1)

##Repeat for each core