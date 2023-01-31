# Targeted-capture-of-eDNA-in-sediment-cores
The following study used environmental DNA buried in sediment cores to reconstruct historical changes in coastal plant communities. 

## Project overview
Sediment cores (~1m long and 7.5cm wide) were collected from a wetland site in South Australia which today is dominated by mangroves. The sediment cores were sub-sectioned and samples were sent for either dating, stable isotope analysis or eDNA analysis. The code and descriptions here refer to the eDNA analysis method and bioinformatic processing. 

## Targeted Capture of eDNA
Targeted capture of eDNA was conducted on sediment core samples using RNA baits/probes designed to capture plant communities in the Alismatales order (see bait design - https://doi.org/10.1101/2021.09.06.456727). A description of this process and comparisons to tradtitional metabarcoding can be found in Foster et al. 2020 https://doi.org/10.1071/MF19175. 
Briefly, DNA was extracted using the DNeasy Powerlyzer Soil Kit (QIAGEN), library prep was conducted using NEBNext Ultra II Library preparation kit (New England Biolabs®) and targeted capture was performed using myBaits® Targeted NGS Manual Version 4.01. Full methods can be found in the manuscript https://doi.org/10.3389/fevo.2021.735744

## Bioinformatic processing
Raw sequences reads were received from the Illumina HiSeq X Ten and demultiplexed using Illumina Bcl2fastq v2.18.0. See Soil_core_processing_code_commandline.R file for next steps, then Soil_core_processing_code_run_with_R.R

Briefly, reads are run through the 'Paleomix' pipeline https://paleomix.readthedocs.io/en/stable/introduction.html See ".yaml" file for settings. This conducts read trimming and filtering then read mapping.
A custom reference database was created using NCBI and sequences from Foster et al. 2021  https://doi.org/10.1002/ece3.8816 these sequences were separated by the different gene regions so that mapping was region specific. See "Reference_library_targetgenes.fasta" file
Resulting bam files are then processed using bcftools and consensus sequences are created for the mapped reads.
This resulted in multiple consensus sequences across multiple genes and species. However, given not all gene regions discriminate at the species level, an additional clustering step was conducted to ensure sample reads were mapped to the correct taxa and taxonomic ranking.
Clustering was performed at 95% similarity using cd-hit-est and both the consensus sample sequences and the reference database. The custom R code "Soil_core_processing_code_run_with_R.R" was then run to take the output of cd-hit-est, assign upstream taxonomy to sequences, then determine whether samples consensus sequences clustered at the species, genus, family or order level and assign a taxonomic ranking.

![Picture 1](https://user-images.githubusercontent.com/122473452/215876420-1a55ba45-203a-4b8d-9679-cf81c5615389.png)

If the sample consensus sequence clustered broadly with other sequences, it meant either, the gene in question did not have high enough interspecific variation or the sample consensus sequence did not contain enough genetic data to generate informative sequences i.e. read depth did not meet the assigned threshold to call a base (depth <50 bp) and instead missing data values were inserted (N’s).
This approach helps to prevent false positives and is conservative in calling species presence.
The 20 Chloropast gene regions used in this study increased the chances of recovering DNA for all species present and improved species assignment because chloroplast genes have different levels of discrimination for different plant species. 

## Scripts

### Soil_core_processing_code_commandline.R 
For command line processing; Paleomix, bcftools consensus call and cd-hit-est clustering

### Soil_core_processing_code_run_with_R.R
For R processing; converting .clstr files from above script to dataframe for manipulation; determining upstream taxonomy; assigning taxonomic ranking; converting data to a phyloseq object for plotting

### .yaml
Example of .yaml input and conditions used with Paleomix




