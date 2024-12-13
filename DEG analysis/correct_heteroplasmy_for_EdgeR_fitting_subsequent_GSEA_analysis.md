correct heteroplasmy for EdgeR fitting and subsequent GSEA analysis
================

## EdgeR YW2 RM refitting

Corrected heteroplasmy for d6 and d15 heteroplasmy are 35% and 22%.

### EdgeR fitting

``` r
library(edgeR)
```

    ## Loading required package: limma

``` r
rawCounts<-read.table('5022D_rawCounts.txt',header=T,row.names=1)
##remove YW2 
rawCounts <- rawCounts[,-2]

dge <- DGEList(counts=rawCounts)
#filter
keep <- filterByExpr(dge)
```

    ## Warning in filterByExpr.DGEList(dge): All samples appear to belong to the same
    ## group.

``` r
dge <- dge[keep, ]
#normalize
dge <- calcNormFactors(dge)


heteroplasmy_levels_corr <- c(0, 1, 1, 0, 0, 0.6, 0.6)#unit: every 35% heteroplasmy change
#heteroplasmy_levels_corr <- c(0, 0.4, 0.4, 0, 0, 0.2, 0.2) #unit: every 100% heteroplasmy 
design <- model.matrix(~heteroplasmy_levels_corr)

dge <- estimateDisp(dge, design)

#fit the model
fit <- glmFit(dge, design)
lrt <- glmLRT(fit, coef="heteroplasmy_levels_corr")
res <- topTags(lrt, n=Inf)
#res$table

##annotation
# To load the unique biomart_data from the saved CSV file
biomart_data_unique <- read.csv("intermediate_results/biomart_data_unique.csv")
# Convert row names of 'res$table' to a column
res$table$ensembl_gene_id <- rownames(res$table)

# Merge the results with the biomart_data_unique by ensembl_gene_id
merged_res <- merge(res$table, biomart_data_unique, by = "ensembl_gene_id", all.x = TRUE)

# to set the row names of the merged results back to the Ensembl gene IDs:
rownames(merged_res) <- merged_res$ensembl_gene_id
# Drop the ensembl_gene_id column as it's now redundant
merged_res$ensembl_gene_id <- NULL

# #export the intermediate results
# write.csv(merged_res, "intermediate_results/edgeR_YM2_RM_correct_heteroplasmy_34_25.csv", row.names=TRUE)
```

### GSEA for KEGG

``` r
merged_results_imported <- read.csv("intermediate_results/edgeR_YM2_RM_correct_heteroplasmy_34_25.csv", header=TRUE,row.names=1)

merged_results_imported$signed_log_pvalue <- -log10(merged_results_imported$PValue) * sign(merged_results_imported$logFC) #using raw p value

#ranked
merged_results_imported <- merged_results_imported[order(merged_results_imported$signed_log_pvalue, decreasing = TRUE), ]

ranked_genes_values <- merged_results_imported$signed_log_pvalue
names(ranked_genes_values) <- merged_results_imported$hgnc_symbol


library(fgsea)
library(msigdbr)

# Get the MSigDB gene sets for Homo sapiens
msigdb_sets <- msigdbr(species = "Homo sapiens")
# Extract the C2 curated gene sets as an example
pathways <- split(msigdb_sets[msigdb_sets$gs_cat == "C2",]$gene_symbol, msigdb_sets[msigdb_sets$gs_cat == "C2",]$gs_name)

# Run GSEA
gsea_results <- fgsea(pathways, ranked_genes_values, minSize=15, maxSize=500, nproc=1)
```

    ## Warning in preparePathwaysAndStats(pathways, stats, minSize, maxSize,
    ## gseaParam, : There are duplicate gene names, fgsea may produce unexpected
    ## results.

    ##   |                                                                              |                                                                      |   0%  |                                                                              |                                                                      |   1%  |                                                                              |=                                                                     |   1%  |                                                                              |=                                                                     |   2%  |                                                                              |==                                                                    |   2%  |                                                                              |==                                                                    |   3%  |                                                                              |===                                                                   |   4%  |                                                                              |===                                                                   |   5%  |                                                                              |====                                                                  |   5%  |                                                                              |====                                                                  |   6%  |                                                                              |=====                                                                 |   6%  |                                                                              |=====                                                                 |   7%  |                                                                              |=====                                                                 |   8%  |                                                                              |======                                                                |   8%  |                                                                              |======                                                                |   9%  |                                                                              |=======                                                               |  10%  |                                                                              |========                                                              |  11%  |                                                                              |========                                                              |  12%  |                                                                              |=========                                                             |  12%  |                                                                              |=========                                                             |  13%  |                                                                              |==========                                                            |  14%  |                                                                              |==========                                                            |  15%  |                                                                              |===========                                                           |  15%  |                                                                              |===========                                                           |  16%  |                                                                              |============                                                          |  17%  |                                                                              |============                                                          |  18%  |                                                                              |=============                                                         |  18%  |                                                                              |=============                                                         |  19%  |                                                                              |==============                                                        |  19%  |                                                                              |==============                                                        |  20%  |                                                                              |==============                                                        |  21%  |                                                                              |===============                                                       |  21%  |                                                                              |===============                                                       |  22%  |                                                                              |================                                                      |  22%  |                                                                              |================                                                      |  23%  |                                                                              |=================                                                     |  24%  |                                                                              |=================                                                     |  25%  |                                                                              |==================                                                    |  25%  |                                                                              |==================                                                    |  26%  |                                                                              |===================                                                   |  27%  |                                                                              |===================                                                   |  28%  |                                                                              |====================                                                  |  28%  |                                                                              |====================                                                  |  29%  |                                                                              |=====================                                                 |  29%  |                                                                              |=====================                                                 |  30%  |                                                                              |=====================                                                 |  31%  |                                                                              |======================                                                |  31%  |                                                                              |======================                                                |  32%  |                                                                              |=======================                                               |  32%  |                                                                              |=======================                                               |  33%  |                                                                              |========================                                              |  34%  |                                                                              |========================                                              |  35%  |                                                                              |=========================                                             |  35%  |                                                                              |=========================                                             |  36%  |                                                                              |==========================                                            |  37%  |                                                                              |==========================                                            |  38%  |                                                                              |===========================                                           |  38%  |                                                                              |===========================                                           |  39%  |                                                                              |============================                                          |  40%  |                                                                              |=============================                                         |  41%  |                                                                              |=============================                                         |  42%  |                                                                              |==============================                                        |  42%  |                                                                              |==============================                                        |  43%  |                                                                              |==============================                                        |  44%  |                                                                              |===============================                                       |  44%  |                                                                              |===============================                                       |  45%  |                                                                              |================================                                      |  45%  |                                                                              |================================                                      |  46%  |                                                                              |=================================                                     |  47%  |                                                                              |=================================                                     |  48%  |                                                                              |==================================                                    |  48%  |                                                                              |==================================                                    |  49%  |                                                                              |===================================                                   |  49%  |                                                                              |===================================                                   |  50%  |                                                                              |===================================                                   |  51%  |                                                                              |====================================                                  |  51%  |                                                                              |====================================                                  |  52%  |                                                                              |=====================================                                 |  52%  |                                                                              |=====================================                                 |  53%  |                                                                              |======================================                                |  54%  |                                                                              |======================================                                |  55%  |                                                                              |=======================================                               |  55%  |                                                                              |=======================================                               |  56%  |                                                                              |========================================                              |  56%  |                                                                              |========================================                              |  57%  |                                                                              |========================================                              |  58%  |                                                                              |=========================================                             |  58%  |                                                                              |=========================================                             |  59%  |                                                                              |==========================================                            |  60%  |                                                                              |===========================================                           |  61%  |                                                                              |===========================================                           |  62%  |                                                                              |============================================                          |  62%  |                                                                              |============================================                          |  63%  |                                                                              |=============================================                         |  64%  |                                                                              |=============================================                         |  65%  |                                                                              |==============================================                        |  65%  |                                                                              |==============================================                        |  66%  |                                                                              |===============================================                       |  67%  |                                                                              |===============================================                       |  68%  |                                                                              |================================================                      |  68%  |                                                                              |================================================                      |  69%  |                                                                              |=================================================                     |  69%  |                                                                              |=================================================                     |  70%  |                                                                              |=================================================                     |  71%  |                                                                              |==================================================                    |  71%  |                                                                              |==================================================                    |  72%  |                                                                              |===================================================                   |  72%  |                                                                              |===================================================                   |  73%  |                                                                              |====================================================                  |  74%  |                                                                              |====================================================                  |  75%  |                                                                              |=====================================================                 |  75%  |                                                                              |=====================================================                 |  76%  |                                                                              |======================================================                |  77%  |                                                                              |======================================================                |  78%  |                                                                              |=======================================================               |  78%  |                                                                              |=======================================================               |  79%  |                                                                              |========================================================              |  79%  |                                                                              |========================================================              |  80%  |                                                                              |========================================================              |  81%  |                                                                              |=========================================================             |  81%  |                                                                              |=========================================================             |  82%  |                                                                              |==========================================================            |  82%  |                                                                              |==========================================================            |  83%  |                                                                              |===========================================================           |  84%  |                                                                              |===========================================================           |  85%  |                                                                              |============================================================          |  85%  |                                                                              |============================================================          |  86%  |                                                                              |=============================================================         |  87%  |                                                                              |=============================================================         |  88%  |                                                                              |==============================================================        |  88%  |                                                                              |==============================================================        |  89%  |                                                                              |===============================================================       |  90%  |                                                                              |================================================================      |  91%  |                                                                              |================================================================      |  92%  |                                                                              |=================================================================     |  92%  |                                                                              |=================================================================     |  93%  |                                                                              |=================================================================     |  94%  |                                                                              |==================================================================    |  94%  |                                                                              |==================================================================    |  95%  |                                                                              |===================================================================   |  95%  |                                                                              |===================================================================   |  96%  |                                                                              |====================================================================  |  97%  |                                                                              |====================================================================  |  98%  |                                                                              |===================================================================== |  98%  |                                                                              |===================================================================== |  99%  |                                                                              |======================================================================|  99%  |                                                                              |======================================================================| 100%

    ## Warning in fgseaMultilevel(pathways = pathways, stats = stats, minSize =
    ## minSize, : For some pathways, in reality P-values are less than 1e-50. You can
    ## set the `eps` argument to zero for better estimation.

``` r
significant_pathways <- gsea_results[gsea_results$padj < 0.05, ]

#make summary of misgDB sets to merge to get hierarchical info
library(dplyr)
```

    ## 
    ## Attaching package: 'dplyr'

    ## The following objects are masked from 'package:stats':
    ## 
    ##     filter, lag

    ## The following objects are masked from 'package:base':
    ## 
    ##     intersect, setdiff, setequal, union

``` r
msigdb_summary <- msigdb_sets %>%
  group_by(gs_name) %>%
  summarize(description = first(gs_description), gs_subcat = first(gs_subcat))  # keep gs description column

significant_pathways <- merge(significant_pathways, msigdb_summary, by.x="pathway", by.y="gs_name")

#order pathway by NES
significant_pathways <- significant_pathways[order(significant_pathways$NES), ]


#hierarchicals
unique_subcats <- unique(significant_pathways$gs_subcat)
#unique_subcats


library(dplyr)
library(ggplot2)

plot_gsea_for_group <- function(data, title) {
    data <- arrange(data, desc(NES))  # Explicitly order data by NES in descending order using dplyr's arrange
    data$pathway <- factor(data$pathway, levels=rev(data$pathway))  # Set the factor levels explicitly in reverse order
    ggplot(data, aes(pathway, NES)) +
      geom_bar(stat="identity", aes(fill=padj), width=0.7) +
      coord_flip() +
      scale_fill_gradient2(low="blue", high="red", mid="white", midpoint=0.025, limit=c(0,0.05), name="Adj. P-value") +
      labs(title=title, x="Pathway", y="Normalized Enrichment Score (NES)") +
      theme_minimal() +
      theme(axis.text.y=element_text(size=6))
}


# 
# # Save plots to a PDF
# pdf("intermediate_results/edgeR_YW2_rm_correct_heteroplasmy_GSEA_results_by_hierarchy.pdf", width=10, height=8)
# 
# 
# # Since unique_cats only contains "C2", we can focus on the subcategories
# for(subcat in unique_subcats) {
#     subcat_data <- subset(significant_pathways, gs_subcat == subcat)
#     p <- plot_gsea_for_group(subcat_data, paste("C2", subcat, sep=": "))
#     print(p)
# }
# 
# dev.off()  # Close the PDF device


##for KEGG ploting#########################################################
##plot gsea to relocate the pathway names
## use description for plot
plot_gsea_kegg <- function(data, title) {
    
    data <- arrange(data, desc(NES))
    data$description <- factor(data$description, levels=rev(data$description))
    
    ggplot(data, aes(description, NES)) +
      geom_bar(stat="identity", aes(fill=padj), width=0.7) +
      geom_text(aes(label=description, y=ifelse(NES > 0, -0.1, 0.1), 
                    hjust=ifelse(NES > 0, 1, 0)), size=3) +
      coord_flip() +
      scale_fill_gradient2(low="red", mid="purple", high="blue", midpoint=0.025, limit=c(0,0.05), name="Adj. P-value") +
      labs(title=title, x="", y="Normalized Enrichment Score (NES)") +  # Empty x-axis label since we're labeling bars directly
      theme_minimal() +
      theme(axis.text.y=element_blank(),  # Hide y-axis text since we're labeling bars directly
            axis.ticks.y=element_blank())
}



# # Extract KEGG pathways subset
# kegg_data <- subset(significant_pathways, gs_subcat == "CP:KEGG")


# #export the all the significant kegg pathways below:
# library(writexl)
# 
# # Convert list columns to string representations (if they exist)
# if("leadingEdge" %in% colnames(kegg_data)) {
#   kegg_data$leadingEdge <- sapply(kegg_data$leadingEdge, function(x) paste(x, collapse=","))
# }
# # Save as XLSX
# write_xlsx(kegg_data, "intermediate_results/KEGG_significant_pathways_edgeR_YW2_RM_correct_heteroplasmy.xlsx")
# 


##read in the kegg_data:


library(readxl)

kegg_data <- read_xlsx("intermediate_results/KEGG_significant_pathways_edgeR_YW2_RM_correct_heteroplasmy.xlsx")


##KEGG to plot
library(dplyr)

# List of pathways to exclude
exclude_pathways <- c("KEGG_CARDIAC_MUSCLE_CONTRACTION", 
                      "KEGG_PATHOGENIC_ESCHERICHIA_COLI_INFECTION", 
                      "KEGG_HYPERTROPHIC_CARDIOMYOPATHY_HCM", 
                      "KEGG_ARRHYTHMOGENIC_RIGHT_VENTRICULAR_CARDIOMYOPATHY_ARVC", 
                      "KEGG_DILATED_CARDIOMYOPATHY",
                      "KEGG_VIBRIO_CHOLERAE_INFECTION")

# Filter out these pathways
filtered_pathways <- kegg_data %>%
  filter(!pathway %in% exclude_pathways)

# Now you can plot using the 'filtered_pathways' dataset


# Plot KEGG pathways
p_kegg <- plot_gsea_kegg(filtered_pathways, "C2: KEGG")


#print(p_kegg)

ggsave(filename="plots and tables to present/kegg_GSEA_edgeR_YW2_RM_correct_heteroplasmy.png", plot=p_kegg, width=8, height=5) #width 8 is minimal to keep all pathway names
```


