library(GAMBLR)
library(tidyverse)
library(data.table)

genome_build <- snakemake@wildcards[["genome_build"]]
in_maf <- snakemake@input[["maf"]]
outfile <- snakemake@output[["maf"]]

# Load the input maf file for this sample
maf <- fread_maf(in_maf)

coding_types <- c(
  "Frame_Shift_Del",
  "Frame_Shift_Ins",
  "In_Frame_Del",
  "In_Frame_Ins",
  "Missense_Mutation",
  "Nonsense_Mutation",
  "Nonstop_Mutation",
  "Silent",
  "Splice_Region",
  "Splice_Site",
  "Translation_Start_Site"
)

# Load the correct aSHM regions file for the current genome_build
if (str_detect(genome_build, "grch37|hg19|hs37d5")){
  ashm_regions <- grch37_ashm_regions[1:3]
} else {
  ashm_regions <- hg38_ashm_regions[1:3]
}

# Remove the chr prefix if necessary
if (!str_detect("chr", maf$Chromosome[1])){
  ashm_regions <- mutate(ashm_regions, 
                         chr_name = str_remove(chr_name, "chr"))
}

# Rename the columns and set the keys for foverlaps
colnames(ashm_regions) <- c("Chromosome", "Start_Position", "End_Position") 
ashm_regions <- data.table(ashm_regions)
setkey(ashm_regions, Chromosome, Start_Position, End_Position)

# Subset the maf to aSHM regions using foverlaps
subset_maf <- foverlaps(maf, ashm_regions) %>% 
  filter(!is.na(Start_Position) | 
           Variant_Classification %in% coding_types) %>% 
  select(-Start_Position, -End_Position) %>% 
  select(Start_Position = i.Start_Position, 
         End_Position = i.End_Position, 
         everything()) %>% 
  select(all_of(colnames(maf)))

# Write to file
write_tsv(subset_maf, outfile)

