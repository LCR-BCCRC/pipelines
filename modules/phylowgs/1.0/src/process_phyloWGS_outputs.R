#'
#' processing phyloWGS outputs pipeline that takes output json files and preprocessing output files
#' SAMPLE_ID.mutass.zip file must be unzipped before runing the script

# E example: how to run
# mkdir -p output
# Rscript ./process.R --samplename SAMPLE_ID -j SAMPLE_ID.summ.json -t unziped.mutass/ -s ssm_data.txt -c cnv_data.txt -a SAMPLE_ID--matched_slms-3.final_deblacklisted_augmented.maf -b SAMPLE_ID_matched_slms-3.final_deblacklisted_augmented.maf -m SAMPLE_ID.muts.json -o out

##################################################
# load required libraries
##################################################

# library("optparse")
library("rjson")
library("tidyverse")
library("ggrepel")
library("data.table")


##########################
#### Snakemake Input #####
##########################

samplename <- snakemake@wildcards[["patient_id"]]
json_file <- snakemake@input[["summ"]]
trees_out <- snakemake@input[["mutass"]]
ssm_file <- snakemake@input[["ssms"]]
cnv_file <- snakemake@input[["cnvs"]]
mafs <- unlist(strsplit(snakemake@params[["maf_list"]], ","))
mut_file <- snakemake@input[["muts"]]
driver_genes <- snakemake@params[["drivers"]]
sample_order <- unlist(strsplit(snakemake@params[["sample_order"]], ","))
genome_build <- snakemake@wildcards[["genome_build"]]

# Define the chr_prefix parameter based on the genome_build
chr_prefixed <- str_detect(genome_build, "hg")


##################################################
# Process input files
###################################################
# Parse the input file and obtain the required data for this run
result1 <- fromJSON(file = json_file)
result_mut <- fromJSON(file = mut_file)
ssm_pre <- read.table(file = ssm_file, header = TRUE)
cnv_pre <- read.delim(file = cnv_file, header = TRUE)[, c("cnv", "a", "d")]



##################################################
# define output files
##################################################
out_json_to_Rtable <- snakemake@output[["tree_summary"]]
ssm_to_trees <- snakemake@output[["maf"]]
cnv_to_trees <- snakemake@output[["cnvs"]]
cellular_prevalence_plot <- file.path(snakemake@output[["plots"]], paste0(samplename, "_cellular_prevalence.pdf"))
CCF_plot <- file.path(snakemake@output[["plots"]], paste0(samplename, "_cancer_cell_fraction_.pdf"))
VAF_plot <- file.path(snakemake@output[["plots"]], paste0(samplename, "_vaf_.pdf"))
VAF_coding_plot <- file.path(snakemake@output[["plots"]], paste0(samplename, "_vaf_coding.pdf"))
tree_plot <- file.path(snakemake@output[["plots"]], paste0(samplename, "_tree.pdf"))
CCF_table <- snakemake@output[["CCF"]]

if (!dir.exists(snakemake@output[["plots"]])) {
  dir.create(snakemake@output[["plots"]])
}

# out_json_to_Rtable= file.path(output_dir, paste("out_res_",samplename,"_json_converted_toR.table", sep = ""))
# ssm_to_trees= file.path(output_dir, paste("out_res_",samplename,"_assigned_ssms_to_best_tree_maf_format.table", sep = ""))
# cnv_to_trees= file.path(output_dir, paste("out_res_",samplename,"_assigned_cnvs_to_best_tree_maf_format.table", sep = ""))
# cellular_prevalence_plot= file.path(output_dir, paste("cellular_prevalence_",samplename,".pdf", sep = ""))
# CCF_plot= file.path(output_dir, paste("cancer_cell_fraction_",samplename,".pdf", sep = ""))
# VAF_plot= file.path(output_dir, paste("vaf_",samplename,".pdf", sep = ""))
# VAF_coding_plot= file.path(output_dir, paste("vaf_ccoding",samplename,".pdf", sep = ""))


###################################################
# open summ.json file and convert it into humam readable format
###################################################

# this function opens SAMPLE_ID_summ.jason and converts it into R table
open_tree <- function(json_summ_file, out_json_to_Rtable) {
  out_res <- NULL
  for (j in 1:length(json_summ_file[["trees"]])) {
    tree_focal <- json_summ_file[["trees"]][j]
    tree_focal_statA <- as.data.frame(t(unlist(sapply(tree_focal, function(x) x[c("clustering_index", "branching_index", "llh", "linearity_index")]))))
    colnames(tree_focal_statA) <- c("clustering_index", "branching_index", "llh", "linearity_index")
    tree_focal_statA$tree_id <- j - 1
    rownames(tree_focal_statA) <- NULL


    tree_focal_statB <- as.data.frame(sapply(tree_focal, function(x) x[3]))[1, -c(3, 6, 9, 12, 15, 18, 21, 24, 27, 30)]
    # tree_focal_statB<-as.data.frame(sapply(tree_focal,function(x)x[3]))[1,!(grepl("cellular_prevalence",colnames(tree_focal_statB)))]
    colnames(tree_focal_statB) <- sub("^[^.]*.", "", colnames(tree_focal_statB))
    stat_both <- cbind(tree_focal_statA, tree_focal_statB)
    out_res <- bind_rows(stat_both, out_res)
    out_res_ordered <- out_res[order(out_res$tree_id), ]
  } # for j loop

  density <- json_summ_file["tree_densities"]
  density_unlist <- data.frame("density" = unlist(density))
  row.names(density_unlist) <- sub("^[^.]*.", "", row.names(density_unlist))

  density_unlist$tree_id <- row.names(density_unlist)
  row.names(density_unlist) <- NULL

  final_table <- merge(out_res_ordered, density_unlist, by.x = "tree_id", by.y = "tree_id") ## add tree densities to all tress table
  write.table(final_table, file = out_json_to_Rtable, col.names = TRUE, row.names = FALSE, sep = "\t", quote = FALSE)
  return(final_table)
}


result_tree <- open_tree(result1, out_json_to_Rtable)


###################################################
# extrcats the best tree
###################################################
# the best tree is the tree with the highest density

best_tree_id <- function(R_table, density) {
  best <- R_table[which.max(R_table$density), ]
  best_tree_focal_name <- best$tree_id
  best_tree_id <- paste(best$tree_id, "json", sep = ".")
  return(best_tree_id)
  return(best_tree_focal_name)
}
best_tree_fileID <- best_tree_id(result_tree, density)


#######################################################################
# extract the stats (SNvs and CNVs assigned to each population) from the best tree
#######################################################################

open_best_tree <- function(trees_out, best_tree_id) {
  unzip(trees_out, files = best_tree_id, exdir = dirname(trees_out), overwrite = TRUE)
  best_tree_path <- paste0(dirname(trees_out), "/", best_tree_id)
  rr <- fromJSON(file = best_tree_path)
  return(rr)
}
rr <- open_best_tree(trees_out, best_tree_fileID)


#######################################################################
# annotate point mutations and CNVs in the best tree
#######################################################################
best_focal <- result1[["trees"]][as.numeric(gsub(".json", "", best_tree_fileID)) + 1]
tree_structure <- as.data.frame(sapply(best_focal, function(x) x["structure"])) ## [6]
tree_roots <- best_focal[[1]]$structure$`0`


merge_both <- function(result1, best_tree_fileID, tree_structure) {
  best_tree <- as.numeric(gsub(".json", "", best_tree_fileID))
  best_focal <- result1[["trees"]][best_tree + 1]
  tree_focal_statB <- as.data.frame(sapply(best_focal, function(x) x["populations"])) ## [3]
  qq <- tree_focal_statB[, grep("cellular_prevalence", colnames(tree_focal_statB))] %>%
    rownames_to_column("sample") %>%
    pivot_longer(-sample,
      names_to = "population",
      values_to = "cellular_prevalence"
    ) %>%
    mutate(population = str_remove(str_remove(population, ".*populations[.]"), "[.]cellular_prevalence")) %>%
    mutate(is_root = ifelse(population %in% tree_roots, TRUE, FALSE)) %>%
    group_by(sample) %>%
    mutate(
      purity = sum(cellular_prevalence[is_root]),
      CCF = cellular_prevalence / purity
    )

  return(qq)
}

both_samples <- merge_both(result1, best_tree_fileID, tree_structure)


write_tsv(both_samples, CCF_table)


ssm <- function(stat_best_tree, ssm_pre, ssm_to_trees, tree_structure, maf_list) {
  out_res_ssm <- NULL
  for (i in 1:length(stat_best_tree$mut_assignments)) {
    focal <- (stat_best_tree$mut_assignments)[i]

    focal_ssms <- data.frame(sapply(focal, function(x) x[1]))

    colnames(focal_ssms) <- sub("^[^.]*.", "", colnames(focal_ssms))
    focal_ssms$phyloWGS_population <- i
    ssm_assign <- merge(ssm_pre, focal_ssms, by.x = "id", by.y = "ssms")[, c("id", "gene", "phyloWGS_population")]
    ssm_assign_spi <- separate(ssm_assign, col = gene, into = c("Chromosome", "Start_Position"), sep = "_", convert = FALSE) %>%
      mutate(Start_Position = as.numeric(Start_Position))
    if (chr_prefixed) {
      ssm_assign_spi$Chromosome <- str_c("chr", as.character(ssm_assign_spi$Chromosome))
    }


    out_res_ssm <- rbind(ssm_assign_spi, out_res_ssm)
  } ## i loop

  ssm_assign_with_maf <- lapply(maf_list, function(x) {
    maf <- read_tsv(x,
      col_types = cols(Chromosome = col_character())
    ) %>%
      # PhyloWGS changes the start position of deletions. This makes the maf start position match that in the PhyloWGS SSM table.
      mutate(Start_Position = ifelse(Variant_Type == "DEL", Start_Position - 1, Start_Position))
    maf <- out_res_ssm %>%
      left_join(maf, by = c("Chromosome", "Start_Position")) %>%
      # Restore the true MAF start postion after the hack above
      mutate(Start_Position = ifelse(Variant_Type == "DEL", Start_Position + 1, Start_Position)) %>%
      select(colnames(maf), everything())
  })
  out_res_ssm <- rbindlist(ssm_assign_with_maf) %>%
    mutate(clonal_status = case_when(
      phyloWGS_population %in% tree_roots & length(tree_roots) > 1 ~ "polyclonal",
      phyloWGS_population %in% tree_roots & length(tree_roots) == 1 ~ "clonal",
      TRUE ~ "subclonal"
    ))

  return(out_res_ssm)
}

ss <- ssm(rr, ssm_pre, ssm_to_trees, tree_structure, mafs)

write_tsv(ss, ssm_to_trees, na = "")


###########################################################
## load mut file to extrcat CNVs start and end positions
###########################################################

cnv <- function(stat_best_tree, cnv_pre, mutation_file, cnv_to_trees) {
  out_res_cnv <-
    bind_rows(lapply(1:length(stat_best_tree$mut_assignments), function(x) {
      data.frame(cnvs = stat_best_tree$mut_assignments[[x]]$cnvs) %>% mutate(phyloWGS_population = x)
    }))



  # return(out_res_cnv)
  out_res_mut <- NULL
  for (cn in 1:length(result_mut$cnvs)) {
    focal_mut_cnv <- (result_mut$cnvs)[cn]

    focal_mut <- data.frame(sapply(focal_mut_cnv, function(x) x[1]))[1, ]
    colnames(focal_mut) <- sub("^[^.]*.", "", colnames(focal_mut))
    focal_mut$cnv_id <- names(focal_mut_cnv)
    out_res_mut <- bind_rows(focal_mut, out_res_mut)
  } ## cn loop

  both_cnvs <- merge(out_res_cnv, out_res_mut, by.x = "cnvs", by.y = "cnv_id") %>%
    select(
      cnvs, phyloWGS_population, physical_cnvs.chrom,
      physical_cnvs.start, physical_cnvs.end,
      physical_cnvs.major_cn, physical_cnvs.minor_cn, physical_cnvs.cell_prev
    )
}


cnv <- cnv(rr, cnv_pre, result_mut, cnv_to_trees)
write.table(cnv, file = cnv_to_trees, col.names = TRUE, row.names = FALSE, sep = "\t", quote = FALSE)

##################### plot the results ####################
###########################################################
#### Slope chart the best tree, cellular prevalence #######

plot_cp <- function(both_samples, cellular_prevalence_plot) {
  pdf(cellular_prevalence_plot, width = 8, height = 8)
  plotA <- ggplot(data = both_samples, aes(x = sample, y = cellular_prevalence, group = population)) +
    geom_line(aes(color = population), size = 2) +
    labs(title = paste("Best Tree", gsub(".json", "", best_tree_fileID), sep = " ")) +
    geom_point(aes(color = population), size = 4) +
    #  Labelling as desired
    xlab("Sample") +
    ylab("Cellular prevalence") +
    theme(
      panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
      panel.background = element_blank(), axis.line = element_line(colour = "black"),
      legend.key = element_rect(fill = NA, colour = NA, size = 0.25), axis.text = element_text(size = 16), axis.title = element_text(size = 18),
      plot.margin = margin(1, 1., 1, 1.5, "cm"), axis.title.y = element_text(margin = margin(t = 70, r = 20, b = 50, l = 10))
    )
  print(plotA)
  dev.off()
}

plot_cp(both_samples, cellular_prevalence_plot)

###########################################################
#### Slope chart the best tree, CCF #######

plot_cp <- function(both_samples, CCF_plot) {
  pdf(CCF_plot, width = 8, height = 8)
  plotB <- ggplot(data = both_samples[both_samples$population != 0, ], aes(x = sample, y = CCF, group = population)) +
    geom_line(aes(color = population), size = 2) +
    labs(title = paste("Best Tree", gsub(".json", "", best_tree_fileID), sep = " ")) +
    geom_point(aes(color = population), size = 4) +
    #  Labelling as desired
    xlab("Sample") +
    ylab("CCF") +
    theme(
      panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
      panel.background = element_blank(), axis.line = element_line(colour = "black"),
      legend.key = element_rect(fill = NA, colour = NA, size = 0.25), axis.text = element_text(size = 16), axis.title = element_text(size = 18),
      plot.margin = margin(1, 1., 1, 1.5, "cm"), axis.title.y = element_text(margin = margin(t = 70, r = 20, b = 50, l = 10))
    )
  print(plotB)
  dev.off()
}

plot_cp(both_samples, CCF_plot)


#############################################
##### Slope chart the best tree (VAF) #######


plot_vaf <- function(ss, VAF_plot) {
  pdf(VAF_plot, width = 8, height = 8)
  plotC <- ss %>%
    select(Hugo_Symbol, Chromosome, Tumor_Sample_Barcode, Start_Position, t_depth, t_alt_count, populations = phyloWGS_population) %>%
    mutate(
      VAF = t_alt_count / t_depth,
      populations = as.factor(populations)
    ) %>%
    filter(!is.na(Tumor_Sample_Barcode)) %>%
    ggplot(aes(
      x = Tumor_Sample_Barcode,
      y = VAF,
      group = interaction(populations, Start_Position),
      color = populations
    )) +
    geom_line(aes(color = populations), size = 0.2, alpha = 0.4) +
    labs(title = paste("Best Tree", gsub(".json", "", best_tree_fileID), sep = " ")) +
    xlab("Sample") +
    ylab("VAF") +
    guides(colour = guide_legend(override.aes = list(alpha = 3))) +
    theme(
      panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
      panel.background = element_blank(), axis.line = element_line(colour = "black"),
      legend.key = element_rect(fill = NA, colour = NA, size = 2), axis.text = element_text(size = 16), axis.title = element_text(size = 18),
      plot.margin = margin(1, 1., 1, 1.5, "cm"), axis.title.y = element_text(margin = margin(t = 70, r = 20, b = 50, l = 10))
    )
  print(plotC)
  dev.off()
}



plot_vaf(ss, VAF_plot)

#############################################################
##### Slope chart the best tree (VAF, coding regions, nonsense, missense and splicing sites) #######

drivers <- read_tsv(driver_genes, col_names = "gene") %>% pull(gene)

plot_vaf_coding <- function(maf, VAF_coding_plot) {
  coding <- ss %>%
    select(Hugo_Symbol, HGVSp_Short, Chromosome, Tumor_Sample_Barcode, Variant_Classification, Start_Position, t_depth, t_alt_count, populations = phyloWGS_population) %>%
    mutate(
      VAF = t_alt_count / t_depth,
      populations = as.factor(populations),
      Tumor_Sample_Barcode = factor(Tumor_Sample_Barcode, levels = sample_order)
    ) %>%
    filter(
      !is.na(Tumor_Sample_Barcode),
      !Variant_Classification %in% c("Silent", "RNA", "IGR", "Intron", "5'Flank", "3'Flank", "5'UTR")
    ) %>%
    mutate(label = ifelse(!is.na(HGVSp_Short), str_c(Hugo_Symbol, "_", HGVSp_Short), str_c(Hugo_Symbol, "_", Variant_Classification)))
  pdf(VAF_coding_plot, width = 8, height = 8)
  plotD <- coding %>%
    ggplot(aes(
      x = Tumor_Sample_Barcode,
      y = VAF,
      group = interaction(populations, Start_Position),
      color = populations
    )) +
    geom_line(aes(color = populations), size = 0.5, alpha = 0.4) +
    geom_text_repel(
      data = filter(group_by(coding, Hugo_Symbol, Start_Position), VAF == max(VAF), Tumor_Sample_Barcode == sample_order[1], Hugo_Symbol %in% drivers),
      aes(
        label = label,
        x = Tumor_Sample_Barcode,
        y = VAF
      ),
      nudge_x = -0.2,
      size = 4
    ) +
    geom_text_repel(
      data = filter(group_by(coding, Hugo_Symbol, Start_Position), VAF == max(VAF), Tumor_Sample_Barcode == sample_order[length(sample_order)], Hugo_Symbol %in% drivers),
      aes(
        label = label,
        x = Tumor_Sample_Barcode,
        y = VAF
      ),
      nudge_x = 0.2,
      size = 4
    ) +
    labs(title = paste("Best Tree", gsub(".json", "", best_tree_fileID), sep = " ")) +
    xlab("Sample") +
    ylab("VAF") +
    guides(colour = guide_legend(override.aes = list(alpha = 3))) +
    theme(
      panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
      panel.background = element_blank(), axis.line = element_line(colour = "black"),
      legend.key = element_rect(fill = NA, colour = NA, size = 2), axis.text = element_text(size = 16), axis.title = element_text(size = 18),
      plot.margin = margin(1, 1., 1, 1.5, "cm"), axis.title.y = element_text(margin = margin(t = 70, r = 20, b = 50, l = 10))
    )

  print(plotD)
  dev.off()
}


plot_vaf_coding(ss, VAF_coding_plot)

#############################################
#####       Draw the best tree        #######
#############################################


tree_structure_long <- tree_structure %>%
  pivot_longer(everything(),
    names_to = "parent",
    values_to = "node"
  ) %>%
  mutate(parent = str_remove_all(parent, ".*[.]")) %>%
  distinct()


positions_x <- function(parents) {
  x <- 1:length(unique(parents))
  names(x) <- unique(parents)
  col_vals <- unname(x[parents])
  return(col_vals)
}

tree_structure_long$x <- positions_x(tree_structure_long$parent)

positions_y <- function(tree_df) {
  y <- c("0" = 0.5)
  for (parent in unique(tree_df$parent)) {
    # parent = "1"
    child_index <- 1
    num_children <- nrow(tree_df[tree_df$parent == parent, ])
    if (num_children == 1) {
      child <- tree_df[tree_df$parent == parent, ]$node
      child_y <- unname(y[parent])
      names(child_y) <- child
      y <- c(y, child_y)
    } else {
      children <- tree_df[tree_df$parent == parent, ]$node
      y_max <- unname(y[parent]) + (0.25 / child_index)
      y_min <- unname(y[parent]) - (0.25 / child_index)
      y_range <- seq(y_min, y_max, length.out = length(children))
      names(y_range) <- children
      y <- c(y, y_range)
    }
    child_index <- child_index + 1
  }
  return(y)
}

tree_structure_long$y <- unname(positions_y(tree_structure_long)[as.character(tree_structure_long$node)])

tree_structure_long <- add_row(tree_structure_long, parent = "0", node = 0, x = 0, y = 0.5)

get_ssms <- function(tree_df, best_focal, best_tree_fileID) {
  data <- best_focal[[str_remove_all(best_tree_fileID, "[.].*")]]$populations
  ssm_vec <- c()
  for (node in tree_df$node) {
    # node = "1"
    num_ssms <- data[[as.character(node)]]$num_ssms
    names(num_ssms) <- as.character(node)
    ssm_vec <- c(ssm_vec, num_ssms)
  }
  return(ssm_vec)
}

tree_structure_long$num_ssms <- get_ssms(tree_structure_long, best_focal, best_tree_fileID)[as.character(tree_structure_long$node)]

tree_structure_long <- tree_structure_long %>%
  mutate(parent = as.numeric(parent)) %>%
  left_join(select(tree_structure_long, node, xstart = x, ystart = y),
    by = c("parent" = "node")
  )



plot_tree <- ggplot(
  tree_structure_long,
  aes(
    x = x,
    y = y,
    label = node
  )
) +
  geom_segment(
    inherit.aes = FALSE,
    aes(
      x = xstart,
      xend = x,
      y = ystart,
      yend = y
    )
  ) +
  geom_point(aes(size = num_ssms),
    fill = "white",
    colour = "black",
    pch = 21
  ) +
  geom_text() +
  scale_size(range = c(5, 20)) +
  ylim(0, 1) +
  theme_void() +
  ggtitle(samplename) +
  theme(legend.position = "none")

ggsave(tree_plot, plot_tree, height = 6, width = 6)

############
##### END ##
############
