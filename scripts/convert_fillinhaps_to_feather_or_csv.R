#### CONVERT FILLINFINDHAPLOTYPES PLUGIN HAPLOTYPE RESULTS TO A SINGLE FILE ####
###################    AUTHOR: STEPHANIE COFFMAN    ############################

     ############## OUTPUT FOLDER STRUCTURE FROM THE PLUGIN ##############
     #   Each output folder from the plugin contains several text files  #
     #   -- one per SNP window tested. Each of these text files contains #
     #   the consensus SNP calls by haplotype group in that SNP window.  #
     #   The "extOut" folder contains the names of the samples that were #
     #   assigned to each haplotype group in each SNP window. So         #
     #   together, each SNP window has two results files that we can     #
     #   combine, bringing together window position information, the     #
     #   haplotype groups named in that window and which of the samples  #
     #   belong to each haplotype group                                  #



library(data.table)
library(feather)
library(dplyr)


# 1) BRING IN THE CONSENSUS SNPS BY HAPLOTYPE BLOCK TO GET POSITION INFO AND GROUP COUNTS

# make sure your working directory (where you are running this from) is set to the
# folder that contains all of the output folders from the plugin results
# get subfolder names - example subfolder name: mintaxa1_mxdiv03_mxhet05_c1
# modify the strings to look for based on your subfolder names 
thedirs = list.dirs()[which(!(grepl("extOut", list.dirs())) & grepl("mintaxa1_mxdiv03_mxhet05", list.dirs()))]

# get list of files in each subfolder
thefiles = lapply(thedirs, function(x){
    paste(x, "/", list.files(x, pattern = "txt"), sep = "")
})

hapblocks_files = unlist(thefiles)

hapblocks = list()


for (i in 1:length(hapblocks_files)){
    
    hapblocks[[i]] <- fread(hapblocks_files[[i]], header=T, verbose = F) 
    
    chr = as.integer(gsub("/.*", "", gsub("\\./mintaxa1_mxdiv03_mxhet05_c", "", hapblocks_files[[i]]))) # get the chrom from the file name
    hapblocks[[i]]$chr = chr # add the chrom as a column
    
    block = as.integer(gsub("\\..*", "", gsub(paste(".*gc", chr, "s", sep = ""), "", hapblocks_files[[i]])))
    hapblocks[[i]]$block = block # add the hapblock number as a column
    hapblocks[[i]] = as.matrix(hapblocks[[i]])
}


print("hapblocks_files loaded")


# make a data table to summarize the hapblock information
# start pos, stop pos, nloci, nhapgroups, names of the hapgroups

get_block_info = function(x){
    chr = unique(x[,"chr"])
    block = unique(x[,"block"])
    startpos = gsub(" ", "", x[1,"pos"])
    endpos = gsub(" ", "", x[nrow(x),"pos"])
    nloci = nrow(x)
    nhaps = length(which(grepl(":", colnames(x))))
    
    out = cbind(chr, block, startpos, endpos, nloci, nhaps)
    return(out)
}

hapblock_info = lapply(hapblocks, get_block_info)
hapblock_info = do.call(rbind, hapblock_info)
hapblock_info = data.table(apply(hapblock_info, 2, as.integer))


write_feather(hapblock_info, "hapblock_info_for_fillin_summary.feather")
# fwrite(hapblock_info, "haplotype_info_for_fillin_summary.csv", row.names = F, sep = ",")



# 2) BRING IN THE SAMPLE INFORMATION SO WE KNOW WHAT SAMPLES WERE ASSIGNED TO WHICH HAPLOTYPE GROUPS

# make the extOut subfolder names; modify strings as needed
thedirs = list.dirs()[which(grepl("extOut", list.dirs()) & grepl("mintaxa1_mxdiv03_mxhet05", list.dirs()))]

# get list of files in each extOut subfolder
thefiles = lapply(thedirs, function(x){
    paste(x, "/", list.files(x, pattern = "txt"), sep = "")
})

groupassigns_files = unlist(thefiles)

hap_results = list()

for (i in 1:length(groupassigns_files)){
    hap_results[[i]] <- fread(groupassigns_files[[i]], sep = "\t", header=F, fill = T, verbose = F) 
    
    chr = as.integer(gsub("/.*", "", gsub("\\./mintaxa1_mxdiv03_mxhet05_c", "", groupassigns_files[[i]]))) # get the chrom from the file name
    hap_results[[i]]$chr = chr # add the chrom as a column
    
    block = as.integer(gsub("\\..*", "", gsub(paste(".*gc", chr, "s", sep = ""), "", groupassigns_files[[i]])))
    hap_results[[i]]$block = block # add the hapblock number as a column
    hap_results[[i]] = as.matrix(hap_results[[i]])
    hap_results[[i]][which(hap_results[[i]] == "")] <- NA
}


print("hapgroups_files loaded")


# melt the data for each data table within the list and add in the start and stop position
hap_results = lapply(hap_results, function(x){
    melt(data.table(x), id.vars = c("chr","block","V1") , na.rm = T, value.name = "hmp_sample")
})

hap_results = do.call(rbind, hap_results)
colnames(hap_results)[which(colnames(hap_results) == "V1")] <- "hapname"
hap_results$hmp_sample = gsub(" ", "", hap_results$hmp_sample) # these ar our sample names
hap_results$chr = as.integer(hap_results$chr) # this is our chromosome
hap_results$block = as.integer(hap_results$block) # this is the snp window

# add our position information for each haplotype block / SNP window
hap_results = left_join(hap_results, hapblock_info[,c("chr","block","startpos","endpos")], by = c("chr","block"))

# hmp_samples with "V2" in the variable column are ones that the plugin used to name the haplotype groups
# we will define the haplotype group name using this information and call it the "origin inbred"
origin_inbreds_list = hap_results[which(hap_results$variable == "V2"),] 
colnames(origin_inbreds_list)[which(colnames(origin_inbreds_list) == "hmp_sample")] <- "origin_inbred"

# bring in origin inbred to hap_results
hap_results = left_join(hap_results, origin_inbreds_list[,c("chr","block","hapname","origin_inbred")], by = c("chr","block","hapname"))


write_feather(hap_results, "hap_results_for_fillin_summary.feather")
# fwrite(hap_results, "hap_results_for_fillin_summary.csv", sep = ",", row.names = F)


