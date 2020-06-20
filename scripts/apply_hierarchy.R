#### RELABEL HAPLOTYPE GROUPS BASED ON A HIERARCHY ####
#########    AUTHOR: STEPHANIE COFFMAN    #############

# This script takes in a feather file or csv file containing
# haplotype group assignments for SNP windows across chromosomes
# and renames the haplotype groups based on a hierarchy.
# The hierarchy is an ordered list of all samples that are contained
# in the haplotype data. Within each haplotype group, the sample
# that is highest in the hierarchy will become the name of that
# haplotype group. Renaming the haplotype groups in this manner 
# provides a consistent naming scheme across all SNP windows on 
# all chromosomes. 


library(data.table)
# library(feather) # only needed if you want to export results to feather format
library(dplyr)

# set working directory to the folder that contains the output folders from the plugin
setwd("../data/") #\maize_expvp_haplotypes\data

# bring in hierarchy
hier = fread("./hierarchy/sample_hierarchy.csv")

# bring in haplotype group assignments
myhapsfile = "../output/hap_results.csv"
hap_data = fread(myhapsfile, sep = ",", header = T)
# hap_data = read_feather(myhapsfile) # use this instead of fread() if you have a feather file

colnames(hap_data)[which(colnames(hap_data) == "hmp_sample")] <- "sample"

# confirm that ever sample in the haplotype results is present in the hierarchy file
hier_samples = sort(hier$sample)
hap_samples = sort(unique(hap_data$sample))
stopifnot(hier_samples == hap_samples)


# join hap assignments and hierarchy
hap_data = left_join(hap_data, hier, by = "sample")

# by haplotype group within each haplotype block
#   which sample has the smallest hierarchy integer (highest in the hierarchy)
#   rename the hapgroup based on that sample
#   change the color for that hapgroup to match the new name

hap_data$hapgrp_name = paste(hap_data$chr, hap_data$block, hap_data$hapname, sep = "_")
hap_data_split = split(hap_data, c(hap_data$hapgrp_name))

apply_hier = function(x){
    index_val = which.min(x$order) # highest ranking sample in the hierarchy
    new_hapname = x$sample[index_val]
    new_hapgrp_hex = x$Hex_Code[index_val]
    x$hapgroup = new_hapname
    x$Hex_Code = new_hapgrp_hex
    return(x)
}


hap_data_split = lapply(hap_data_split, apply_hier)
hap_data_out = rbindlist(hap_data_split) 

fwrite(hap_data_out[,c("chr","block","sample","startpos","endpos","hapgroup","Hex_Code"),], 
       "../output/hap_results_new_hier.csv", sep = ",", row.names = F)
# write_feather(hap_data_out[,c("chr","block","sample","startpos","endpos","hapgroup","Hex_Code"),], 
#               "../output/hap_results_new_hier.feather")

