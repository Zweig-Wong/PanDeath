load2(SP = "240716-TCGA_NMF", info = T, path = "Output/NMF/TCGA_PANCAN_nmf.R")
source("Function/Calculation/Dimensionality_Reduction.R", echo=F)
source("Function/Plot/Plot_Heatmap.R", echo=F)
source("Function/Phenotype/Signature_Profile.R", echo=F)

## Figure 1A
TCGA_PANCAN_nmf <- Dimensionality_Reduction_S4(
    TCGA_PANCAN_nmf, 
    Signature_name = "CD_expr_raw", 
    dimension = 2, method = "umap",
    n_neighbors = 15, 
    metric = "euclidean",
    spread = 1, 
    min_dist = 0.01, 
    local_connectivity = 1,
    bandwidth = 1, 
    n_epochs = 500 
)
Cluster_visualization_S4(
    TCGA_PANCAN_nmf, 
    Projection_name = "CD_expr_raw_umap", Projection_method = "UMAP",
    save = T, plotname = "1-Raw_data", cmap = Pancancer_cmap, 
    dir = "Output/Fig1A-Umap/", wh = 1, scale = 1
)
TCGA_PANCAN_nmf <- Dimensionality_Reduction_S4(
    TCGA_PANCAN_nmf, 
    Signature_name = "CD_expr_std_nor", 
    dimension = 2, method = "umap",
    n_neighbors = 100, 
    metric = "euclidean", 
    spread = 0.5, 
    min_dist = 1,
    local_connectivity = 1, 
    bandwidth = 1, 
    n_epochs = 500 
)
Cluster_visualization_S4(
    TCGA_PANCAN_nmf, 
    Projection_name = "CD_expr_std_nor_umap", Projection_method = "UMAP",
    save = T, plotname = "2-Processed_data", cmap = Pancancer_cmap, dir = "Output/Fig1A-Umap/", wh = 1, scale = 1
)

## Figure 1B
TCGA_PANCAN_nmf <- Dimensionality_Reduction_S4(
    TCGA_PANCAN_nmf, 
    Signature_name = "nmf_CD_h", 
    dimension = 2, method = "umap",
    n_neighbors = 100, 
    metric = "euclidean", 
    spread = 0.5, 
    min_dist = 1,
    local_connectivity = 1, 
    bandwidth = 1, 
    n_epochs = 500 
)
Cluster_visualization_S4(
    TCGA_PANCAN_nmf, 
    Projection_name = "nmf_CD_h_umap", Projection_method = "UMAP",
    Group_name = "nmf_CD_predict", Group_items = c("CD_subtype"),
    save = T, plotname = "3-Cluster", cmap = Cluster_cmap, dir = "Output/Fig1A-Umap/", wh = 1, scale = 1
)

## Figure 1C
Proportion_S4(
    TCGA_PANCAN_nmf, 
    Group_name.x = "Clinical", Group_item.x = "Cancer",
    Group_name.y = "nmf_CD_predict", Group_item.y = "CD_subtype", 
    summary = T, test = "Chi", 
    plot = T, plot_style = "Stalk", plot_name = "Cluster_EACH_Cancer", 
    save = T, cmap = Cluster_cmap, wh = 200/75, scale = 2, dir = "Output/Fig1C-Proportion/"
)

## Figure 1D
Expression_mat <- TCGA_PANCAN_nmf@Signature$nmf_CD_h[, c("ID_sample", "CD5", "CD4", "CD3", "CD2", "CD1")] %>% column_to_rownames("ID_sample") %>% as.matrix() %>% t()
Annotation_row_df <- transpose(Score[order(-ID_sample)], keep.names = "CD", make.names = "ID_sample") %>% 
    column_to_rownames("CD") %>% 
    as.matrix() %>% 
    scale() %>% 
    t() %>% 
    as.data.frame() %>% 
    select(c("Pyroptosis", "Necroptosis", "Ferroptosis", "Apoptosis"))

Annotation_column_df <- merge(
    Get_group(TCGA_PANCAN_nmf, "nmf_CD_predict", "CD_subtype"), 
    Get_group(TCGA_PANCAN_nmf, "Clinical", "Cancer"), 
    by = "ID_sample", all.x = T
)[order(CD_subtype, Cancer)] %>% column_to_rownames("ID_sample")
Column_split <- Annotation_column_df[["CD_subtype"]]

Heatmap <- Plot_Heatmap(
    Expression_mat, 
    Annotation_cell_color = colorRamp2(c(0.2, 0.25, 0.3), c("#618cac", "#ffffff", "#f6adb8")),
    Annotation_column_df = Annotation_column_df,
    Annotation_column_color = list(CD_subtype = Cluster_cmap, Cancer = Pancancer_cmap),
    Annotation_row_df = Annotation_row_df,
    Annotation_row_color = list(
        Pyroptosis = colorRamp2(c(-1.5, 0.2, 1), c("#618cac", "#ffffff", "#f6adb8")),
        Necroptosis = colorRamp2(c(-1.5, 0.2, 1), c("#618cac", "#ffffff", "#f6adb8")),
        Ferroptosis = colorRamp2(c(-1.5, 0.2, 1), c("#618cac", "#ffffff", "#f6adb8")),
        Apoptosis = colorRamp2(c(-1.5, 0.2, 1), c("#618cac", "#ffffff", "#f6adb8"))
    ),
    Column_name = F, Row_name = T,
    Column_split = Column_split, 
    plot_name = paste0("Heatmap"), save = T, dir = paste0("Output/Fig1D-Heatmap/")
)

Test_list <- list()
for (Cluster in paste0("Cluster_", c(1:5))) {
    Score_cluster <- Score[Cluster, on = "CD_subtype"][, CD_subtype := NULL] %>% 
        melt(., variable.name = "Cell_death", value.name = "ssGSEA_Score")
    Score_cluster$Cell_death <- factor(Score_cluster$Cell_death, levels = c("Pyroptosis", "Necroptosis", "Ferroptosis", "Apoptosis"))
    ## Plot
        Plot <- ggplot(Score_cluster, mapping = aes_string(x = "Cell_death", y = "ssGSEA_Score", fill = "Cell_death")) + 
            stat_boxplot(mapping = aes_string(x = "Cell_death", y = "ssGSEA_Score", fill = "Cell_death"),
                        width = 0.5, 
                        position = position_dodge(0.5)) +   
            geom_boxplot(aes_string(fill = "Cell_death"),                        
                        position = position_dodge(0.5),                
                        width = 0.5) +
            scale_fill_manual(values = cmap) +
            scale_color_manual(values = cmap) +
            stat_compare_means(aes_string(group = "Cell_death"), ref.group = ".all.", method = "wilcox.test", label = "p.signif", show.legend = F, hide.ns = T) +
            geom_hline(yintercept = 0, linetype = "dashed", color = "gray", size = 1) +
            theme(plot.title = element_text(hjust = 0.5)) +
            theme_pubr() +
            labs(color = "") +
            theme(
                axis.title.x = element_blank(),
                axis.text.x = element_blank(),
                axis.ticks.x = element_blank(),
                axis.title.y = element_blank(),
                axis.text.y = element_blank(),
                axis.ticks.y = element_blank(),
                axis.line = element_blank(),
                legend.position = "none",
                panel.background = element_rect(fill = adjust_transparency(Cluster_cmap[Cluster], alpha = 0.3)),
            ) +
            scale_y_continuous(limits = c(-4, 4), breaks = c(-2, 0, 2))   
        assign(paste0(Cluster, "_CD_Score"), Plot)
        save2(data = paste0(Cluster, "_CD_Score"), dir = "Output/Fig1D-CD_score/", mode = "w", format = "img", wh = 1, scale = 1)
    ## Test
        Test_list[[Cluster]] <- CompareGroups(Score_cluster, "Cell_death", "ssGSEA_Score")
}