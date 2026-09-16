## Figure 5A
ICI_PANCAN_nmf <- Dimensionality_Reduction_S4(
    ICI_PANCAN_nmf, 
    Signature_name = "CD_expr_raw", 
    dimension = 2, 
    method = "umap",
    n_neighbors = 15, 
    metric = "euclidean", 
    spread = 0.5, 
    min_dist = 0.5, 
    local_connectivity = 1, 
    bandwidth = 1, 
    n_epochs = 500 
)
Cluster_visualization_S4(
    ICI_PANCAN_nmf, 
    Projection_name = "CD_expr_raw_umap", Projection_method = "UMAP",
    plot_Group = T, Group_name = "Clinical", Group_items = c("Cancer"), 
    save = T, plotname = "Umap_raw", cmap = Pancancer_cmap, dir = "Output/Fig5A-NMF/Umap/", wh = 1, scale = 1
)
ICI_PANCAN_nmf <- Dimensionality_Reduction_S4(ICI_PANCAN_nmf, 
    Signature_name = "CD_expr_std_nor", 
    dimension = 2, 
    method = "umap",
    n_neighbors = 100, 
    metric = "manhattan", 
    spread = 1, 
    min_dist = 0.125, 
    local_connectivity = 10, 
    bandwidth = 10, 
    n_epochs = 500 
)
Cluster_visualization_S4(
    ICI_PANCAN_nmf, 
    Projection_name = "CD_expr_std_nor_umap", Projection_method = "UMAP",
    plot_Group = T, Group_name = "Clinical", Group_items = c("Cancer"),
    save = T, plotname = "Umap_processed", cmap = Pancancer_cmap, dir = "Output/Fig5A-NMF/Umap/", wh = 1, scale = 1
)
ICI_PANCAN_nmf <- Dimensionality_Reduction_S4(
    ICI_PANCAN_nmf, 
    Signature_name = "nmf_CD_h", 
    dimension = 2, 
    method = "umap",
    n_neighbors = 100, 
    metric = "manhattan", 
    spread = 0.75, 
    min_dist = 0.25, 
    local_connectivity = 10, 
    bandwidth = 10, 
    n_epochs = 500 
)
Cluster_visualization_S4(
    ICI_PANCAN_nmf, 
    Projection_name = "nmf_CD_h_umap", Projection_method = "UMAP",
    plot_Group = T, Group_name = "nmf_CD_predict", Group_items = c("CD_subtype"),
    save = T, plotname = "Umap_cluster", cmap = Cluster_cmap, dir = "Output/Fig5A-NMF/Umap/", wh = 1, scale = 1
)

## Figure 5C
Proportion_S4(
    ICI_PANCAN_nmf, 
    Group_name.x = "Clinical", Group_item.x = "Cancer",
    Group_name.y = "nmf_CD_predict", Group_item.y = "CD_subtype", 
    summary = T, 
    plot = T, plot_style = "Stalk", plot_name = "Cluster_EACH_Cancer", 
    save = T, cmap = Cluster_cmap, wh = 1, scale = 1, dir = "Output/Fig5C-Proportion/"
)

## Figure 5D
Expression_mat <- ICI_PANCAN_nmf@Signature$nmf_CD_h[, c("ID_sample", "CD5", "CD4", "CD3", "CD2", "CD1")] %>% 
    column_to_rownames("ID_sample") %>% 
    as.matrix() %>% 
    t()
Annotation_row_df <- transpose(Score[order(-ID_sample)], keep.names = "CD", make.names = "ID_sample") %>% 
    column_to_rownames("CD") %>% 
    as.matrix() %>% 
    scale() %>% 
    t() %>% 
    as.data.frame() %>% 
    select(c("Pyroptosis", "Necroptosis", "Ferroptosis", "Apoptosis"))

Annotation_column_df <- merge(
    Get_group(ICI_PANCAN_nmf, "nmf_CD_predict", "CD_subtype"), 
    Get_group(ICI_PANCAN_nmf, "Clinical", "Cancer"), 
    by = "ID_sample", all.x = T
)[order(CD_subtype, Cancer)] %>% column_to_rownames("ID_sample")
Column_split <- Annotation_column_df[["CD_subtype"]]

Heatmap <- Plot_Heatmap(
    Expression_mat, 
    Annotation_cell_color = colorRamp2(c(0.175, 0.25, 0.325), c("#618cac", "#ffffff", "#f6adb8")),
    Annotation_column_df = Annotation_column_df,
    Annotation_column_color = list(CD_subtype = Cluster_cmap, Cancer = Pancancer_cmap),
    Annotation_row_df = Annotation_row_df,
    Annotation_row_color = list(
        Pyroptosis = colorRamp2(c(-1.5, 0, 1.5), c("#618cac", "#ffffff", "#f6adb8")),
        Necroptosis = colorRamp2(c(-1.5, 0, 1.5), c("#618cac", "#ffffff", "#f6adb8")),
        Ferroptosis = colorRamp2(c(-1.5, 0, 1.5), c("#618cac", "#ffffff", "#f6adb8")),
        Apoptosis = colorRamp2(c(-1.5, 0, 1.5), c("#618cac", "#ffffff", "#f6adb8"))
    ),
    Column_name = F, Row_name = T,
    Column_split = Column_split, wh = 2.25, scale = 1.5,
    plot_name = paste0("Heatmap"), save = T, dir = paste0("Output/Fig5D-Heatmap/")
)
for (Cluster in paste0("Cluster_", c(1:5))) {
    Score_cluster <- Score[Cluster, on = "CD_subtype"][, CD_subtype := NULL] %>% 
        melt(., variable.name = "Cell_death", value.name = "ssGSEA_Score")
    Score_cluster$Cell_death <- factor(Score_cluster$Cell_death, levels = c("Pyroptosis", "Necroptosis", "Ferroptosis", "Apoptosis"))
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
            axis.title.y = element_blank(),
            axis.line = element_blank(),
            legend.position = "none",
            panel.background = element_rect(fill = adjust_transparency(Cluster_cmap[Cluster], alpha = 0.3)),
        ) +
        scale_y_continuous(limits = c(-4, 4), breaks = c(-2, 0, 2))   
}

## Figure 5E
load2(SP = "240716-TCGA_NMF", path = "Output/NMF/nmf_CD_w.R")
TCGA_W <- nmf_CD_w
setnames(TCGA_W, paste0(names(TCGA_W), "_TCGA"))
setnames(TCGA_W, "Gene_TCGA", "Gene")
TCGA_W_mat <- TCGA_W %>% column_to_rownames("Gene") %>% as.matrix()
load2(SP = "260303-ICI_NMF", path = "Output/2-NMF/nmf_CD_w.R")
ICI_W <- nmf_CD_w
setnames(ICI_W, paste0(names(ICI_W), "_ICI"))
setnames(ICI_W, "Gene_ICI", "Gene")
ICI_W_mat <- ICI_W[TCGA_W$Gene, on = "Gene"] %>% column_to_rownames("Gene") %>% as.matrix()
Correlation <- corr.test(x = ICI_W_mat, y = TCGA_W_mat, method = "spearman", adjust = "holm", minlength = 100000000)$r
col_fun <- colorRamp2(c(-1, 0, 1), c("#618cac", "#ffffff", "#f6adb8"))
cell_fun <- function(j, i, x, y, width, height, fill){
    if(i >= j) {
        grid.circle(x = x, y = y, r = abs((Correlation[i, j])/2) * min(unit.c(width, height)), 
            gp = gpar(fill = col_fun(Correlation[i, j]), col = NA))
        grid.text(sprintf("%.2f", Correlation[i, j]), x, y, gp = gpar(fontsize = 15))
    }
}
Correlation_Heatmap <- Heatmap(
    Correlation, 
    name = "Correlation", 
    col = col_fun, na_col = "white",
    rect_gp = gpar(type = "none"), 
    cell_fun = cell_fun,
    show_column_names = F,
    show_row_names = F,
    cluster_columns = F,
    cluster_rows = F,
    row_names_side = 'left',
    border = T,
    show_heatmap_legend = F
)

## Figure 5F
Signature_Landscape(
    object = ICI_PANCAN_nmf, 
    Signature_names = c("DAMPR", "Cytokine", "CTL_CAF"),
    Group_name = "nmf_CD_predict", Group_item = "CD_subtype", 
    dir = "Output/Fig5F-TME/", plot_name = "Immune",
    wh = 1.25,
    scale = 1
)

## Figure 5G
Proportion_S4(
    ICI_PANCAN_nmf, 
    Group_name.x = "nmf_CD_predict", Group_item.x = "CD_subtype",
    Group_name.y = "Clinical", Group_item.y = "Response_binary_merge",
    summary = F, test = "Chi", plot = T, plot_style = "Sanky", cmap = Response_binary_cmap, 
    save = T, wh = 1, scale = 1, plot_name = "Response_binary_merge", dir = "Output/Fig5G-Response/"
)

## Figure 5H
Survival_S4(
    ICI_PANCAN_nmf, 
    Event = "OS", mode = "KM", 
    Group_name = "nmf_CD_predict",
    summary = F, plot = T, cmap = Cluster_cmap,
    merge_plot = F, plotname = "OS", save = T, dir = "Output/Fig5H-Outcome/"
)
RMST_list <- list()
for (Cluster in paste0("Cluster_", 1:5)) {
    Group_dt <- ICI_PANCAN_nmf@Group$nmf_CD_predict[, .(
        ID_patient, 
        CD_subtype = fifelse(CD_subtype == Cluster, Cluster, "Others")
    )]
    Group_dt$CD_subtype <- factor(Group_dt$CD_subtype, levels = c("Others", Cluster))
    Group_dt <- merge(ICI_PANCAN_nmf@Group$Clinical, Group_dt, by = "ID_patient")
    data <- na.omit(merge(ICI_PANCAN_nmf@Clinical, Group_dt, by = "ID_patient", all = T)[, .(ID_patient, OS, OS_time, CD_subtype)])
    fit <- rmst2(
        time = data$OS_time,
        status = data$OS,
        arm = ifelse(data$CD_subtype == Cluster, 1, 0),
        tau = 4*365
    )
    RMST_list[[Cluster]] <- data.table(  
        RMST = fit$unadjusted.result[1, 1],
        CI_Low = fit$unadjusted.result[1, 2],
        CI_High = fit$unadjusted.result[1, 3],
        P = fit$unadjusted.result[1, 4],
        tau = fit$tau
    )
}
RMST_table <- rbindlist(RMST_list, idcol = "Cluster")

## Figure 4I
KM_list <- list()
for (Cluster in paste0("Cluster_", 1:5)) {
    Group_dt <- ICI_PANCAN_nmf@Group$nmf_CD_predict[, .(
        ID_patient, ID_sample, 
        CD_subtype = fifelse(CD_subtype == Cluster, Cluster, "Others")
    )]
    ICI_PANCAN_nmf <- Add_group(
        ICI_PANCAN_nmf, 
        from_Group_dt = T, Group_name = Cluster, Group_dt = Group_dt,
        Factor_level_list = list(Cluster = c(Cluster, "Others"))
    )
    Cluster_Survival <- Survival_S4(
        ICI_PANCAN_nmf, 
        Event = "OS", mode = "KM", 
        Group_name = Cluster, Group_item = "CD_subtype", 
        summary = T, plot = T
    )
    KM_result <- data.table(
        median_Cluster = Cluster_Survival[Group == Cluster]$`50_quantile`,
        median_Others = Cluster_Survival[Group == "Others"]$`50_quantile`,
        p_value = Cluster_Survival[1, ]$log_rank_pvalue
    )
    KM_result[, median_Diff := median_Cluster - median_Others]
    KM_result[, p_signif := fcase(
        p_value < 0.05, "*", 
        p_value < 0.01, "**", 
        p_value < 0.001, "***", 
        default = ""
    )]
    KM_list[[Cluster]] <- KM_result
}
KM_data <- rbindlist(KM_list, idcol = "Group")  
Plot_median_survival <- ggplot(KM_data) + 
    geom_linerange(mapping = aes(x = Group, ymin = median_Others, ymax = median_Cluster, color = Group), show.legend = F, linewidth = 2) +
    geom_text(mapping = aes(x = Group, y = median_Cluster + 1, label = round(p_value, digits = 3), color = "black"), size = 5, show.legend = F) +
    scale_color_manual(values = Cluster_cmap) +
    scale_y_continuous(limits = c(100, 1000)) +
    xlab("") +
    ylab("Median Survival (Days)") +
    theme_bw() +
    theme(
        axis.line = element_blank(),
        legend.position = "none",
        panel.border = element_rect(colour = "black", fill = NA, size = 1),
    ) +
    coord_flip()

Cox_list <- list()
for (Cluster in paste0("Cluster_", 1:5)) {
    Group_dt <- ICI_PANCAN_nmf@Group$nmf_CD_predict[, .(
        ID_patient, 
        CD_subtype = fifelse(CD_subtype == Cluster, Cluster, "Others")
    )]
    Group_dt$CD_subtype <- factor(Group_dt$CD_subtype, levels = c("Others", Cluster))
    Group_dt <- merge(ICI_PANCAN_nmf@Group$Clinical, Group_dt, by = "ID_patient")
    Cox_list[[Cluster]] <- Survival(
        Survival_dt = ICI_PANCAN_nmf@Clinical[, .(ID_patient, OS, OS_time)], Event = "OS", mode = "multicox", 
        Group_name = "Cluster", Group_dt = Group_dt, Group_items = c("Cancer", "CD_subtype"),
        summary = T, plot = F, save = F
    )["CD_subtype", on = "Type"]
}
Cox_table <- rbindlist(Cox_list)
Cox_table$HR_log2 <- log2(Cox_table$HR)
Cox_table$CI_upper_log2 <- log2(Cox_table$CI_upper)
Cox_table$CI_lower_log2 <- log2(Cox_table$CI_lower)
Plot_Cox <- ggplot(Cox_table[var != "Others", ]) + 
    geom_hline(yintercept = 0, linetype = "dashed", linewidth = 1) +
    geom_linerange(mapping = aes(x = var, ymin = CI_lower_log2, ymax = CI_upper_log2, color = var), show.legend = F, linewidth = 1) +
    geom_point(mapping = aes(x = var, y = HR_log2, color = var), size = 5, show.legend = F) +
    geom_text(mapping = aes(x = var, y = HR_log2, label = round(p_value, digits = 3), color = "black"), size = 5, show.legend = F) +
    scale_color_manual(values = Cluster_cmap) +
    scale_y_continuous(limits = c(-1.22, 1.22)) +
    xlab("") +
    ylab("Log2(Harzard ratio)") +
    theme_bw() +
    theme(
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.line = element_blank(),
        legend.position = "none",
        panel.border = element_rect(colour = "black", fill = NA, size = 1),
    ) +
    coord_flip()

Chisq_outcome <- list()
for (Cluster in paste0("Cluster_", 1:5)) {
    Group_dt <- ICI_PANCAN_nmf@Group$nmf_CD_predict[, .(
        ID_sample, ID_patient,
        CD_subtype = fifelse(CD_subtype == Cluster, Cluster, "Others")
    )]
    ICI_PANCAN_nmf <- Add_group(
        ICI_PANCAN_nmf, 
        from_Group_dt = T, Group_name = Cluster, Group_dt = Group_dt,
        Factor_level_list = list(Cluster = c(Cluster, "Others"))
    )
    Proportion_S4(
        ICI_PANCAN_nmf, 
        Group_name.x = Cluster, Group_item.x = "CD_subtype",
        Group_name.y = "Clinical", Group_item.y = "Response_binary_merge",
        summary = F, plot = T, plot_style = "Sanky", cmap = Response_binary_cmap, 
        save = T, wh = 2, scale = 0.5, plot_name = Cluster, dir = "Output/Fig5I-Median/"
    )
    table <- dcast(na.omit(merge(
        ICI_PANCAN_nmf@Group[[Cluster]][, .(ID_patient, CD_subtype)],
        ICI_PANCAN_nmf@Group$Clinical[, .(ID_patient, Response_binary_merge)], 
        by = "ID_patient"
    ))[, .N, by = c("CD_subtype", "Response_binary_merge")], CD_subtype~Response_binary_merge, value.var = "N")[order(CD_subtype)] %>% 
    column_to_rownames("CD_subtype") %>% as.matrix()
    result <- chisq.test(table)
    Chisq_outcome[[Cluster]] <- data.table(Statistics = result$statistic, p.value = result$p.value)
}
Chisq_outcome <- rbindlist(Chisq_outcome, idcol = "Cluster")

## Figure 5J
DAMPR <- ICI_PANCAN_nmf@Signature$DAMPR[Group_dt$ID_sample, on = "ID_sample"][, .(ID_sample, Receptor)]
Signature_Profile(
    Signature_dt = DAMPR, 
    Group_dt = Group_dt, Group_item = "Collect_time", paired = T,
    plot = T, save = T, cmap = c("Pre" = "#057d96", "Post" = "#9e0000"), 
    plot_style = "boxplot_merge", plot_name = "DAMPR", wh = 1, dir = "Output/Fig5J-PrePos/"
)

## Figure 5K
Score <- Score[Group_dt$ID_sample, on = "ID_sample"][, .(ID_sample, Pyroptosis, Necroptosis, Ferroptosis, Apoptosis)]
Signature_Profile(
    Signature_dt = Score, 
    Group_dt = Group_dt, Group_item = "Collect_time", paired = T,
    plot = T, save = T, cmap = c("Pre" = "#057d96", "Post" = "#9e0000"), 
    plot_style = "boxplot_merge", plot_name = "PCD_Score", wh = 1, dir = "Output/Fig5K-PrePost/"
)

## Figure 5L
Pre_Post <- ggplot(data = data, aes_string(x = "Collect_time", stratum = "CD_subtype", alluvium = "subject", y = "N", fill = "CD_subtype")) +
    geom_flow() +
    geom_stratum(width = 0.25, alpha = 0.5) +
    geom_text(aes(label = paste0(sprintf("%0.1f%%", round(data$Frequency, digits = 4)*100))), stat = "stratum", size = 3) +
    scale_fill_manual(values = Cluster_cmap) +
    scale_x_discrete(expand = c(0.1, 0.1)) +
    theme_pubr()

## Figure 5M
Differential(
    ICI_PANCAN_nmf@Signature$CD_expr_std[Group_dt$ID_sample, on = "ID_sample"][, -c("ID_patient")], 
    Group_dt = Group_dt, Group_item = "Collect_time", Group_level = c("Pre", "Post"),
    Genes_in_interest = 20, method = "limma", paired = T,
    logFC = 0.2, p = 0.05,
    summary = F, plot = T, save = T, dir = "Output/Fig5F-PrePost/"
) 
