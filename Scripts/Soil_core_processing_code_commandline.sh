################################################################################################################
######## The below code takes raw sequenced reads from Illumina Hiseq X ten after undergoing targeted capture using 
######## specially designed RNA baits to capture plant communities in the Alismatales order 
######## (see bait design - https://doi.org/10.1101/2021.09.06.456727)
######## 
####
##
#

### This code is run through a HPC using the command line program BASH
### The first part of the code uses a program called paleomix https://paleomix.readthedocs.io/en/stable/introduction.html
### Refer to the above documentation for how to install and run the program
### Instructions and an example of the .yaml file are provided to adapt to your data
### Create a separate .yaml file for each sample and save in a folder together

### Download Adapter Removal https://adapterremoval.readthedocs.io/en/stable/installation.html
### Create barcode list as desribed in Adapter Removal above
### Use adapter removal to demutiplex raw reads 

# AdapterRemoval --file1 HybCap46_2_p4_p5_S5_L008_R1_001.fastq.gz --file2 HybCap46_2_p4_p5_S5_L008_R2_001.fastq.gz #your own read file names
# --basename Nic02 #name of folder where reads are stores
# --adapter1 AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTNNNNNNNNGTGTAGATCTCGGTGGTCGCCGTATCATT #adjust to own adapters
# --adapter2 AGATCGGAAGAGCACACGTCTGAACTCCAGTCACNNNNNNNNATCTCGTATGCCGTCTTCTGCTTG 
# --barcode-list Demuxfile_poll.txt #use own barcode list
# --barcode-mm 1 --demultiplex-only


#-----------------------------------------------------------------------------------------------------------------------
### The below codes takes the demultiplexed raw fastq files and runs them through the program paleomix where files are
### filtered and trimmed then mapped using bwa to the reference database.
### The resulting bam output files are then sorted and indexed using samtools and variants are called usng bcftools mpileup
### reads are then filtered for coverage  < 50 reads, normalised and a consensus sequence is called.
### Ambiguity in base calls (N) are removed from the consensus sequences and reference database and these are combined
### The combined consensus sequences and reference data base are then clustered using cd-hit-est and the output is then
### exported to a folder and manipulated in R


## Set up .sh file to fun through HPC or command line interface using .yaml files set up from paleomix above
## Download/install paleomix using the instructions provided on their website
## Download samtools, bcftools and bioawk from conda and cd-hit-est from https://sites.google.com/view/cd-hit?pli=1
## Ensure .yaml files for each sample are in the same folder as your reference database

#!/bin/bash
#SBATCH -p batch        # partition (this is the queue your job will be added to)
#SBATCH -N 1            # number of nodes (due to the nature of sequential processing, here uses single node)
#SBATCH -n 1            # number of cores
#SBATCH --time=24:00:00 # time allocation, which has the format (D-HH:MM)
#SBATCH --mem=64GB      # memory pool for all cores (here set to 32 GB)

# sample=`ls *.yaml`
# 
# for i in ${sample}; do
#   sample=${i%%.yaml*}
#   paleomix bam_pipeline run ${sample}.yaml #running paleomix command
# 
# done
# 
# sample=`ls *.bam` #select output .bam files from paleomix pipeline output
# reference_file="Reference_library_targetgenes" #reference database
# 
# for i in ${sample}; do
# sample=${i%%.bam*}
# echo ${sample}
# samtools sort ${sample}.bam -o ${sample}_sorted.bam
# samtools index ${sample}_sorted.bam
# 
# bcftools mpileup -B -Ou -Q 30 -q 30 -f ${reference_file}.fasta ${sample}_sorted.bam | bcftools call -Ob --ploidy 1 -m -o ${sample}_sorted.vcf.gz
# 
# bcftools index ${sample}_sorted.vcf.gz
# 
# bedtools genomecov -ibam ${sample}_sorted.bam -bga | awk '($4 < 50)' > ${sample}_flt.bed
# 
# bcftools norm -f ${reference_file}.fasta -m +any -Oz ${sample}_sorted.vcf.gz -o ${sample}_sorted_norm.vcf.gz
# 
# bcftools index ${sample}_sorted_norm.vcf.gz
# 
# cat ${reference_file}.fasta | bcftools consensus ${sample}_sorted_norm.vcf.gz -H 1pIu --mask ${sample}_flt.bed -o ${sample}.fa
# done
# 
# sample=`ls *.fa`
# sed '/^[^>]/s/[N]//g' ${reference_file}.fasta > ${reference_file}_noN.fasta
# 
# for i in ${sample}; do
# 
# sample=${i%%.fa*}
# #Add'Soil' to all sample sequence names
# sed 's/^>\(.*\)$/>Soil_\1/'  ${sample}.fa > ${sample}_soil.fa
# 
# ## get rid of N in reference and fasta file
# 
# sed '/^[^>]/s/[N]//g' ${sample}_soil.fa > ${sample}_soil_trimN.fa
# 
# #trim sequences <100bp
# 
# bioawk -c fastx '{ if(length($seq) > 100) { print ">"$name; print $seq }}' ${sample}_soil_trimN.fa > ${sample}_soil_trimN_filtered.fa
# 
# cat *.fasta > ${file_name}.fasta # to combine output of sample with references for clustering
# 
# cat ${sample}_soil_trimN_filtered.fa ${reference_file}_noN.fasta > ${sample}_cd-hit.fasta
# 
# cd-hit-est -i ${sample}_cd-hit.fasta -o output_${sample} -c 0.95 -n 9 -d 200 -g 1 -r 1 -M 1600 -aS 1 -aL 0.1
# 
# done
# 
# mkdir final_fasta_files
# mv *.fa final_fasta_files
# mv *.fasta final_fasta_files
# 
# mkdir final_cdhit_output
# mv *.clstr final_cdhit_output
#-------------------------------------------------------------------------------------------------------------------------------------------

