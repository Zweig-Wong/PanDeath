source("Function/Plot/Plot_Heatmap.R", echo=F)
source("Function/Phenotype/Signature_Profile.R", echo=F)

## Figure 4A
DAMPR_marker <- Get_marker(DAMPR, subset = list(
    Receptor = "[Category == 'Receptor']"
))
DAMPR_expr <- RNA[, c("ID_sample", DAMPR_marker$Gene), with = FALSE]

DAMPR_expr_mat <- DAMPR_expr %>% column_to_rownames("ID_sample") %>% as.matrix() %>% scale() %>% t()
Annotation_row_df <- DAMPR_marker[, .(Gene, Term)] %>% 
    column_to_rownames("Gene")
Annotation_column_df <- TCGA_PANCAN_nmf@Group$nmf_CD_predict[, .(ID_sample, CD_subtype)] %>% 
    column_to_rownames("ID_sample")    
Heatmap <- Plot_Heatmap(
    DAMPR_expr_mat, 
    Annotation_cell_color = colorRamp2(
        seq(-1.5, 1.5, length.out = 10), 
        c("#374f97", "#5377ae", "#85a7cb", "#bacfe2", "#e3eaef", 
            "#f0e5da", "#e4bea2", "#c88467", "#a1493d", "#811925")
    ),
    Annotation_row_df = Annotation_row_df,
    Annotation_column_df = Annotation_column_df,
    Column_name = F, Row_name = T,
    Row_split = Annotation_row_df[["Term"]], Row_cluster = T,
    Column_split = Annotation_column_df[["CD_subtype"]], Column_cluster = F,
    plot_name = "DAMPR_heatmap", save = T, dir = "Output/Fig4A-DAMPR/", wh = 1.1, scale = 2.5
)

## Figure 4B
IM_marker <- Get_marker(Immune_Modulator)
IM_expr <- TCGA_PANCAN_RNA@Signature$RNA_raw[, c("ID_sample", IM_marker$Gene), with = FALSE]
IM_expr <- merge(IM_expr, TCGA_PANCAN_nmf@Group$nmf_CD_predict, by = "ID_sample")[, -c("ID_sample")]
IM_expr <- IM_expr[, lapply(.SD, median, na.rm = TRUE), by = "CD_subtype", .SDcols = -c("CD_subtype")][order(CD_subtype)]
IM_expr_mat <- IM_expr %>% column_to_rownames("CD_subtype") %>% as.matrix() %>% scale() %>% t()
Annotation_row_df <- IM_marker[, .(Gene, Term)] %>% 
    column_to_rownames("Gene")
Annotation_row_df$Term <- factor(Annotation_row_df$Term, levels = unique(Annotation_row_df$Term))
Heatmap <- Plot_Heatmap(
    IM_expr_mat, 
    Annotation_cell_color = colorRamp2(
        seq(-1.5, 1.5, length.out = 10), 
        c("#374f97", "#5377ae", "#85a7cb", "#bacfe2", "#e3eaef", "#f0e5da", "#e4bea2", "#c88467", "#a1493d", "#811925")
    ),
    Annotation_row_df = Annotation_row_df,
    Column_name = T, Row_name = T,
    Row_split = Annotation_row_df[["Term"]], Row_cluster = F,
    plot_name = "Immune Modulator", save = T, dir = "Output/Fig4B-Modulators/", wh = 1/2, scale = 2
)

## Figure 4C
Immune_cell <- c(
    "Neutrophils",
    "T_Cells_CD8", 
    "T_Cells_CD4", 
    "NK_Cells",
    "B_Cells", 
    "Plasma_Cells", 
    "Dendritic_Cells",
    "Macrophages"
)
Data <- merge(TCGA_PANCAN_nmf@Group$nmf_CD_predict, TCGA_PANCAN_nmf@Signature$Immune_cell[, c("ID_sample", Immune_cell), with = FALSE], all.x = T, by = "ID_sample")
Mean <- Data[, lapply(.SD, mean, na.rm = TRUE), .SDcols = Immune_cell, by = "CD_subtype"] 
Immune_cell_value_mat <- Mean[, c("CD_subtype", ..Immune_cell)] %>% 
    column_to_rownames("CD_subtype") %>% as.matrix() %>% scale() %>% t()
Immune_cell_value_mat <- Immune_cell_value_mat[, paste0("Cluster_", 1:5)]

CompareGroups <- function(dt, Group_item, Signature_item) {
    results <- rbindlist(lapply(unique(dt[[Group_item]]), function(g) {data.table(
        Group = g, p.value = t.test(
        dt[get(Group_item) == g, get(Signature_item)],
        dt[get(Group_item) != g, get(Signature_item)],
        alternative = "two.sided", exact = FALSE
    )$p.value)}))
    results[, p.adj := p.adjust(p.value, method = "BH")]
    results[, signif := fcase(p.adj < 0.0001, "****", p.adj < 0.001, "***", p.adj < 0.01, "**", p.adj < 0.05, "*", default = "ns")]
    results <- results[order(Group)]
    return(results)
}
Test_list <- list()
for (Cell in Immune_cell) {
    Test_list[[Cell]] <- CompareGroups(Data, "CD_subtype", Cell)
}
Test <- rbindlist(Test_list, idcol = "Signature")[, p_adj_scale := fcase(
    p.adj >= 0.05, 0.25,
    p.adj < 0.05 & p.adj >= 0.01, 0.5,
    p.adj < 0.01 & p.adj >= 0.001, 0.75,
    p.adj < 0.001, 1
)]

Immune_cell_p_mat <- dcast(Test, Group ~ Signature, value.var = "p_adj_scale")[
    paste0("Cluster_", 1:5), on = "Group"][, c("Group", ..Immune_cell)] %>% 
    column_to_rownames("Group") %>% 
    as.matrix() %>% 
    t()
Immune_cell_p_mat <- Immune_cell_p_mat[, paste0("Cluster_", 1:5)]

col_fun <- colorRamp2(
        seq(-1, 1, length.out = 10), 
        c("#374f97", "#5377ae", "#85a7cb", "#bacfe2", "#e3eaef", 
            "#f0e5da", "#e4bea2", "#c88467", "#a1493d", "#811925")
    )
cell_fun <- function(j, i, x, y, width, height, fill){
    grid.circle(x = x, y = y, r = abs(Immune_cell_p_mat[i, j]/3) * min(unit.c(width, height)), 
        gp = gpar(fill = col_fun(Immune_cell_value_mat[i, j]), col = NA))
}
Immune_cell_Summary <- Heatmap(
    Immune_cell_value_mat, 
    name = "Immune_cell", 
    col = col_fun, na_col = "white",
    rect_gp = gpar(type = "none"), 
    cell_fun = cell_fun,
    cluster_columns = F,
    cluster_rows = F,
    row_names_side = 'left',
    border = T,
    show_heatmap_legend = T
)

## Figure 4D
CTL <- TCGA_PANCAN_Immune$Immune_thearapy[, c("ID_sample", "TIDE_CTL")]
TCGA_PANCAN_nmf <- Add_signature(TCGA_PANCAN_nmf, Signature_dt = CTL)
Signature_Profile_S4(
    TCGA_PANCAN_nmf, 
    Signature_name = "CTL", 
    Group_name = "nmf_CD_predict", 
    plot = T, cmap = Cluster_cmap, save = T, plot_style = "violin_split", plot_name = "CTL", dir = "Output/Fig4D/", wh = 1
)
CAF <- TCGA_PANCAN_Immune$Immune_thearapy[, c("ID_sample", "TIDE_CAF")]
TCGA_PANCAN_nmf <- Add_signature(TCGA_PANCAN_nmf, Signature_dt = CAF)
Signature_Profile_S4(
    TCGA_PANCAN_nmf, 
    Signature_name = "CAF", 
    Group_name = "nmf_CD_predict", 
    plot = T, cmap = Cluster_cmap, save = T, plot_style = "violin_split", plot_name = "CAF", dir = "Output/Fig4D/", wh = 1
)
for (pathway in c("TGFb", "INFg")) {
    Score <- Calculate_Score(Marker = Immune_Pathway@Marker[pathway, on = "Term"], Expr = RNA, method = "ssGSEA", save = F)
    sourcedata_list[[pathway]] <- Score
    TCGA_PANCAN_nmf <- Add_signature(TCGA_PANCAN_nmf, Signature_dt = Score, Signature_name = pathway)
    Signature_Profile_S4(
        TCGA_PANCAN_nmf, 
        Signature_name = pathway, 
        Group_name = "nmf_CD_predict", 
        plot = T, cmap = Cluster_cmap, save = T, plot_style = "violin_split", plot_name = pathway, dir = "Output/Fig4D/", wh = 1
    )
}

## Figure 4E
Cell_cmap <- c(
    "Epithelial" = "#fe4f79",
    "Fibroblast" = "#9f776f",
    "Pericyte" = "#e3ba86",
    "Endothelial" = "#f9c8c9",

    "CD8T" = "#ef8fa2",
    "CD4T" = "#4fb3c9",
    "NK" = "#b36fab",

    "Plasma" = "#00a248",
    "Bcell" = "#a3d768",

    "DC" = "#146fab",
    "Mast" = "#6a8ec9",
    "Mono_Mphi" = "#1e8b9a"
)
Major_umap <- DimPlot(object, reduction = 'umap', group.by = "CellType", label = TRUE, pt.size = 0.75) +
    scale_color_manual(values = Cell_cmap) 

## Figure 4G
Signature_Landscape(
    object = SC_PANCAN, 
    Signature_name = c("PCD_Score"),
    Group_name = "PCD_Subtype", Group_item = "Subtype", 
    test_method = "wilcox.test", 
    dir = "Output/Fig4G_PCD_subtype/", plot_name = "Landscape",
    wh = 0.81,
    scale = 1
)

## Figure 4H
for (item in c("Proliferation", "Stemness", "Cellular_Senescence", "EMT")) {
    Signature_Profile_S4(
        SC_PANCAN, 
        Signature_name = "Tumor_pathway_sample", Signature_item = item,
        Group_name = "PCD_Subtype", Group_item = "Subtype", 
        plot = T, cmap = Cluster_cmap, save = T, plot_style = "violin_split", plot_name = item, dir = "Output/Fig4H_Tumor/", wh = 1
    )
}

## Figure 4I
SC_PANCAN@Signature$DC_marker_sample <- SC_PANCAN@Signature$DC_marker_sample[, c("ID_sample",
    "CLEC9A", "XCR1",  # cDC1
    "CD1C", "SIRPA", # cDC2
    "CCR7", "CLEC4C",  # pDC
    "IRF7", "LAMP3", # migDC
    "TLR2", "NLRP3", "CGAS", # DAMPRs
    "HLA_A", "HLA_B", "HLA_C" # HLA
)]
Signature_Landscape(
    object = SC_PANCAN, 
    Signature_name = c("DC_marker_sample"), 
    Group_name = "PCD_Subtype", Group_item = "Subtype", 
    test_method = "wilcox.test", 
    dir = "Output/Fig4I_DC/", plot_name = "Marker",
    wh = 0.35,
    scale = 2
)
for (item in c("DAMPR", "AP_KEGG")) {
    Signature_Profile_S4(
        SC_PANCAN, 
        Signature_name = "DC_pathway_sample", Signature_item = item,
        Group_name = "PCD_Subtype", Group_item = "Subtype", 
        plot = T, cmap = Cluster_cmap, save = T, plot_style = "violin_split", plot_name = item, dir = "Output/Fig4I_DC/", wh = 1
    )
}

## Figure 4J
SC_PANCAN@Signature$CD8T_marker_sample <- SC_PANCAN@Signature$CD8T_marker_sample[, c("ID_sample",
    "GZMA", "GZMB", "GZMH", "GZMK", "PRF1", "GNLY", # Cytotoxicity
    "ICOS", "TNFRSF9", # Co-stimulator
    "PDCD1", "CTLA4", # Exhausted
    "IFNG", "TGFB1", # Cytokine
    "IL7R", "TCF7" # Naive
)]
Signature_Landscape(
    object = SC_PANCAN, 
    Signature_name = c("CD8T_marker_sample"), 
    Group_name = "PCD_Subtype", Group_item = "Subtype", 
    test_method = "wilcox.test", 
    dir = "Output/Fig4J_CD8T/", plot_name = "Marker2",
    wh = 0.35,
    scale = 2
)
for (item in c("Cytotoxicity", "IFNg")) {
    Signature_Profile_S4(
        SC_PANCAN, 
        Signature_name = "CD8T_pathway_sample", Signature_item = item,
        Group_name = "PCD_Subtype", Group_item = "Subtype", 
        plot = T, cmap = Cluster_cmap, save = T, plot_style = "violin_split", plot_name = item, dir = "Output/Fig4J_CD8T/", wh = 1
    )
}

## Figure 4K
SC_PANCAN@Signature$CD4T_marker_sample <- SC_PANCAN@Signature$CD4T_marker_sample[, c("ID_sample",
    "TNF", "IL2", # Th1
    "CXCL13", "IL21", #Tfh
    "FOXP3", "IL2RA", # Treg
    "ICOS", "TNFRSF9", # Co-stimulator
    "PDCD1", "CTLA4", # Exhausted
    "IFNG", "TGFB1", # Cytokine
    "IL7R", "TCF7" # Naive
)]
Signature_Landscape(
    object = SC_PANCAN, 
    Signature_name = c("CD4T_marker_sample"), 
    Group_name = "PCD_Subtype", Group_item = "Subtype", 
    test_method = "wilcox.test", 
    dir = "Output/Fig4K_CD4T/", plot_name = "Marker2",
    wh = 0.35,
    scale = 2
)
for (item in c("Th1", "IFNg")) {
    Signature_Profile_S4(
        SC_PANCAN, 
        Signature_name = "CD4T_pathway_sample", Signature_item = item,
        Group_name = "PCD_Subtype", Group_item = "Subtype", 
        plot = T, cmap = Cluster_cmap, save = T, plot_style = "violin_split", plot_name = item, dir = "Output/Fig4K_CD4T/", wh = 1
    )
}

## Figure 4L
SC_PANCAN@Signature$CAF_marker_sample <- SC_PANCAN@Signature$CAF_marker_sample[, c("ID_sample",
    "FAP", "DCN",
    "COL1A1", "COL1A2", 
    "TGFB1", "TGFB2", "TGFB3", "LTBP1", "LTBP2",
    "TGFBR1", "TGFBR2", "TGFBR3", "SMAD2", "SMAD3"
)]
Signature_Landscape(
    object = SC_PANCAN, 
    Signature_name = c("CAF_marker_sample"), 
    Group_name = "PCD_Subtype", Group_item = "Subtype", 
    test_method = "wilcox.test", 
    dir = "Output/Fig4L_CAF/", plot_name = "Marker2",
    wh = 0.35,
    scale = 2
)
for (item in c("CAF", "TGFb")) {
    Signature_Profile_S4(
        SC_PANCAN, 
        Signature_name = "CAF_pathway_sample", Signature_item = item,
        Group_name = "PCD_Subtype", Group_item = "Subtype", 
        plot = T, cmap = Cluster_cmap, save = T, plot_style = "violin_split", plot_name = item, dir = "Output/Fig4L_CAF/", wh = 1
    )
}

## Figure 4M
Monocle_time <- plot_cells(
    cds_CD8T, color_cells_by = "pseudotime", 
    label_groups_by_cluster = FALSE, label_cell_groups = FALSE,
    label_branch_points = FALSE, label_roots = TRUE, label_leaves = TRUE, graph_label_size = 1.5, 
    show_trajectory_graph = TRUE, trajectory_graph_color = "black", trajectory_graph_segment_size = 0.75,
    cell_size = 1.25, alpha = 0.5
) + theme(
    axis.text.x = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks.x = element_blank(),
    axis.ticks.y = element_blank(),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    axis.line.x = element_blank(),
    axis.line.y = element_blank(),
    legend.position = "none"
)
Monocle_time$layers[[1]]$aes_params$colour <- "white"
metadata <- as.data.table(colData(cds_CD8T), keep.rownames = "Cell")
metadata[, Subtype2 := fcase(
    Subtype == "PCD-P", "PCD-P",
    Subtype == "PCD-A", "PCD-A",
    default = "Others"
)]
colData(cds_CD8T)$Subtype2 <- metadata$Subtype2
Monocle_subtype <- plot_cells(
    cds_CD8T, color_cells_by = "Subtype2",
    label_groups_by_cluster = FALSE, label_cell_groups = FALSE,
    label_branch_points = FALSE, label_roots = TRUE, label_leaves = TRUE, graph_label_size = 1.5, 
    show_trajectory_graph = TRUE, trajectory_graph_color = "black", trajectory_graph_segment_size = 0.75,
    cell_size = 1.25, alpha = 0.5
) + scale_color_manual(values = c("PCD-P" = "#fdcc83", "PCD-A" = "#6292b0", "Others" = "transparent")) + theme(
    axis.text.x = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks.x = element_blank(),
    axis.ticks.y = element_blank(),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    axis.line.x = element_blank(),
    axis.line.y = element_blank(),
    legend.position = "none"
)
Monocle_subtype$layers[[1]]$aes_params$colour <- "white"

## Figure 4N
Signature_Profile_S4(
    SC_PANCAN, 
    Signature_name = "CD8T_pathway_sample", Signature_item = "Naive",
    Group_name = "PCD_Subtype", Group_item = "Subtype", 
    plot = T, cmap = Cluster_cmap, save = T, plot_style = "violin_split", plot_name = item, dir = dir, wh = 1
)

## Figure 4O
Pseudotime <- merge(
    as.data.table(colData(cds_CD8T), keep.rownames = "Cell")[, .(Cell, ID_sample)],
    as.data.table(pseudotime(cds_CD8T), 
    keep.rownames = "Cell"
)[, .(Cell = V1, Pseudotime = V2)], by = "Cell", all.x = T)[, Cell := NULL][, lapply(.SD, median, na.rm = TRUE), by = "ID_sample"]
Pseudotime[sapply(Pseudotime, is.infinite)] <- NA
Pseudotime$Pseudotime_scaled <- scale(Pseudotime$Pseudotime)
Pseudotime <- merge(Pseudotime, PCD_Subtype, by = "ID_sample", all.x = T)[order(Subtype)]
setcolorder(Pseudotime, c("ID_sample", "Subtype"))
save2(data = Pseudotime, dir = dir, mode = "w", format = "tsv")
Signature_Profile(
    Signature_dt = Pseudotime, Signature_name = "Pseudotime", Signature_item = "Pseudotime_scaled",
    Group_dt = PCD_Subtype, Group_name = "PCD_Subtype", Group_item = "Subtype", 
    plot = T, cmap = Cluster_cmap, save = T, plot_style = "violin_split", plot_name = "Subtype_Pseudotime", dir = dir, wh = 1
)

## Figure 4P
Monocle_time <- plot_cells(
    cds_CD4T, color_cells_by = "pseudotime", 
    label_groups_by_cluster = FALSE, label_cell_groups = FALSE,
    label_branch_points = FALSE, label_roots = TRUE, label_leaves = TRUE, graph_label_size = 1.5, 
    show_trajectory_graph = TRUE, trajectory_graph_color = "black", trajectory_graph_segment_size = 0.75,
    cell_size = 1.25, alpha = 0.5
) + theme(
    axis.text.x = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks.x = element_blank(),
    axis.ticks.y = element_blank(),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    axis.line.x = element_blank(),
    axis.line.y = element_blank(),
    legend.position = "none"
)
Monocle_time$layers[[1]]$aes_params$colour <- "white"
metadata <- as.data.table(colData(cds_CD4T), keep.rownames = "Cell")
metadata[, Subtype2 := fcase(
    Subtype == "PCD-P", "PCD-P",
    Subtype == "PCD-A", "PCD-A",
    default = "Others"
)]
colData(cds_CD4T)$Subtype2 <- metadata$Subtype2
Monocle_subtype <- plot_cells(
    cds_CD4T, color_cells_by = "Subtype2",
    label_groups_by_cluster = FALSE, label_cell_groups = FALSE,
    label_branch_points = FALSE, label_roots = TRUE, label_leaves = TRUE, graph_label_size = 1.5, 
    show_trajectory_graph = TRUE, trajectory_graph_color = "black", trajectory_graph_segment_size = 0.75,
    cell_size = 1.25, alpha = 0.5
) + scale_color_manual(values = c("PCD-P" = "#fdcc83", "PCD-A" = "#6292b0", "Others" = "transparent")) + theme(
    axis.text.x = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks.x = element_blank(),
    axis.ticks.y = element_blank(),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    axis.line.x = element_blank(),
    axis.line.y = element_blank(),
    legend.position = "none"
)
Monocle_subtype$layers[[1]]$aes_params$colour <- "white"

## Figure 4Q
Signature_Profile_S4(
    SC_PANCAN, 
    Signature_name = "CD4T_pathway_sample", Signature_item = "Naive",
    Group_name = "PCD_Subtype", Group_item = "Subtype", 
    plot = T, cmap = Cluster_cmap, save = T, plot_style = "violin_split", plot_name = item, dir = dir, wh = 1
)

## Figure 4R
Pseudotime <- merge(
    as.data.table(colData(cds_CD4T), keep.rownames = "Cell")[, .(Cell, ID_sample)],
    as.data.table(pseudotime(cds_CD4T), 
    keep.rownames = "Cell"
)[, .(Cell = V1, Pseudotime = V2)], by = "Cell", all.x = T)[, Cell := NULL][, lapply(.SD, median, na.rm = TRUE), by = "ID_sample"]
Pseudotime[sapply(Pseudotime, is.infinite)] <- NA
Pseudotime$Pseudotime_scaled <- scale(Pseudotime$Pseudotime)
Pseudotime <- merge(Pseudotime, PCD_Subtype, by = "ID_sample", all.x = T)[order(Subtype)]
setcolorder(Pseudotime, c("ID_sample", "Subtype"))
save2(data = Pseudotime, dir = dir, mode = "w", format = "tsv")
Signature_Profile(
    Signature_dt = Pseudotime, Signature_name = "Pseudotime", Signature_item = "Pseudotime_scaled",
    Group_dt = PCD_Subtype, Group_name = "PCD_Subtype", Group_item = "Subtype", 
    plot = T, cmap = Cluster_cmap, save = T, plot_style = "violin_split", plot_name = "Subtype_Pseudotime", dir = dir, wh = 1
)