pacman::p_load(limma, enrichplot, clusterProfiler, ggridges)
source("Function/Phenotype/Signature_Profile.R", echo=F)
source("Function/Phenotype/Survival.R", echo=F)

## Figure 3A
load2(SP = "240716-TCGA_NMF", path = "Output/NMF/TCGA_PANCAN_nmf.R")
Hallmarks <- read.gmt("../../../System/DataCenter/Marker/Hallmark/Raw/H_Hallmarks.v2024.1.Hs.symbols.gmt")

Pathway_GSEA <- function(
    Expr, 
    Group_name = "", Group_dt = NULL, Group_item = NULL, Group1 = "High", Group2 = "Low", 
    GeneSet, GeneSet_name,
    plot = F, plot_style = "ridgeplot", plot_name = NULL, plot_GenSets, cmap = NULL, wh = 1, scale = 1, save = F, dir = "GSEA/"
) {
    dir.create(dir, recursive = T, showWarnings = F)
    Expr <- transpose(Expr, keep.names = "Gene", make.names = "ID_sample")
    ## Differential analysis
        ## Set contrast
            Group_dt <- Group_dt[colnames(Expr)[-1], on = "ID_sample"]
            group_list <- factor(Group_dt[[Group_item]], levels = c(Group1, Group2))
            design <- model.matrix(~0+group_list)
            colnames(design) <- levels(group_list)
            rownames(design) <- colnames(Expr)[-1]
            cont.matrix <- makeContrasts(contrasts = paste0(Group1, "-", Group2), levels = design)
        ## Limma
            fit <- lmFit(Expr, design)
            fit <- contrasts.fit(fit, cont.matrix)
            fit <- eBayes(fit, trend = TRUE)
        ## Result
            DEG <- na.omit(topTable(fit, coef = 1, n = Inf, adjust.method = "BH", sort.by = "logFC"))
            DEG <- DEG[order(DEG$logFC, decreasing = T), ]
            write.table(DEG, file = paste0(dir, "DEG.tsv"), row.names = F, quote = F, sep = "\t")
    ## GSEA
        ranks <- DEG$logFC
        names(ranks) <- DEG$Gene
        GSEA_result <- GSEA(geneList = ranks, TERM2GENE = GeneSet, pvalueCutoff = 1, verbose = FALSE)
        assign(paste0("GSEA_", GeneSet_name), GSEA_result, envir = .GlobalEnv)
        save2(data = paste0("GSEA_", GeneSet_name), dir = dir, mode = "d")
        GSEA_result_table <- as.data.table(as.data.frame(GSEA_result))
        GSEA_result_table <- GSEA_result_table[order(NES, decreasing = T), ]
        assign(paste0("GSEA_", GeneSet_name, "_table"), GSEA_result_table, envir = .GlobalEnv)
        save2(data = paste0("GSEA_", GeneSet_name, "_table"), dir = dir, mode = "dw", format = "tsv")
    ## Plot
        GSEA_result@result$NES[GSEA_result@result$NES > 3] <- 3
        GSEA_result@result$NES[GSEA_result@result$NES < -3] <- -3
        GSEA_result@result$Description <- factor(GSEA_result@result$Description, levels = c(plot_GenSets, setdiff(unique(GSEA_result@result$Description), plot_GenSets)))
        plot <- ridgeplot(GSEA_result, showCategory = plot_GenSets, orderBy = "Description", decreasing = T, fill = "NES") +
            scale_fill_gradientn(colors = c("#618cac", "#ffffff", "#f6adb8"), values = c(0, 0.5, 1), limit = c(-3, 3)) +
            scale_x_continuous(limits = c(-1.5, 1.5), breaks = seq(-1.5, 1.5, 0.5)) +
            geom_density_ridges(quantile_lines = TRUE, alpha = 0.75, quantiles = 2) +
            theme(
                axis.text.x = element_blank(),
                axis.ticks.x = element_blank(),
                axis.text.y = element_blank(),
                axis.ticks.y = element_blank(),
                legend.position = "none"
            ) 
} 

plot_GenSets_list <- list(
    "Immune" = c(
        "HALLMARK_INFLAMM ATORY_RESPONSE",
        "HALLMARK_INTERFERON_GAMMA_RESPONSE",
        "HALLMARK_INTERFERON_ALPHA_RESPONSE",
        "HALLMARK_TNFA_SIGNALING_VIA_NFKB",
        "HALLMARK_TGF_BETA_SIGNALING",
        "HALLMARK_IL2_STAT5_SIGNALING",
        "HALLMARK_IL6_JAK_STAT3_SIGNALING"
    ),
    "CellCycle" = c(
        "HALLMARK_E2F_TARGETS",
        "HALLMARK_G2M_CHECKPOINT",
        "HALLMARK_MITOTIC_SPINDLE"
    ),
    "Signaling" = c(
        "HALLMARK_MYC_TARGETS_V1",
        "HALLMARK_MTORC1_SIGNALING",
        "HALLMARK_P53_PATHWAY",
        "HALLMARK_KRAS_SIGNALING_UP",
        "HALLMARK_PI3K_AKT_MTOR_SIGNALING"
    ),
    "Metabolism" = c(
        "HALLMARK_DNA_REPAIR",
        "HALLMARK_UNFOLDED_PROTEIN_RESPONSE",
        "HALLMARK_OXIDATIVE_PHOSPHORYLATION",
        "HALLMARK_REACTIVE_OXYGEN_SPECIES_PATHWAY",
        "HALLMARK_FATTY_ACID_METABOLISM"
    ),
    "EMT" = c(
        "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION",
        "HALLMARK_APICAL_JUNCTION",
        "HALLMARK_ANGIOGENESIS"
    )
)

for (Cluster in paste0("Cluster_", 1:5)) {
    Group_dt <- TCGA_PANCAN_nmf@Group$nmf_CD_predict
    Group_dt[CD_subtype != Cluster, ]$CD_subtype <- "Cluster_0"
    for (plot_GenSets_name in names(plot_GenSets_list)) {
        plot_GenSets <- plot_GenSets_list[[plot_GenSets_name]]
        Pathway_GSEA(
            Expr = TCGA_PANCAN_RNA@Signature$RNA_std, 
            Group_name = "nmf_CD_predict", Group_dt = Group_dt, Group_item = "CD_subtype", Group1 = Cluster, Group2 = "Cluster_0",
            GeneSet = Hallmarks, GeneSet_name = "Hallmarks",
            plot = T, plot_style = "ridgeplot", plot_name = plot_GenSets_name,
            plot_GenSets = plot_GenSets,
            wh = 15/(4*length(plot_GenSets)), scale = 2*length(plot_GenSets)/15, save = T,
            dir = paste0("Output/Fig3A-Hallmark/", Cluster, "/")
        )
    }
}

## Figure 3B
Proliferation <- data.table(
    ID_sample = TCGA_PANCAN_nmf_Profiles@Signature$Proliferation_CellCycle$ID_sample,
    Proliferation = rowMeans(TCGA_PANCAN_nmf_Profiles@Signature$Proliferation_CellCycle[, -c("ID_sample")])
)
TCGA_PANCAN_nmf_Profiles <- Add_signature(TCGA_PANCAN_nmf_Profiles, Signature_dt = Proliferation)
Signature_Profile_S4(
    TCGA_PANCAN_nmf_Profiles, 
    Signature_name = "Proliferation", 
    Group_name = "nmf_CD_predict", 
    plot = T, cmap = Cluster_cmap, save = T, plot_style = "violin_split", plot_name = "Proliferation", dir = "Output/Fig3B-Malignancy/", wh = 1
)

Senescence <- data.table(
    ID_sample = TCGA_PANCAN_nmf_Profiles@Signature$Cellular_Senescence$ID_sample,
    Senescence = rowMeans(TCGA_PANCAN_nmf_Profiles@Signature$Cellular_Senescence[, -c("ID_sample")])
)
TCGA_PANCAN_nmf_Profiles <- Add_signature(TCGA_PANCAN_nmf_Profiles, Signature_dt = Senescence)

Signature_Profile_S4(
    TCGA_PANCAN_nmf_Profiles, 
    Signature_name = "Senescence", 
    Group_name = "nmf_CD_predict", 
    plot = T, cmap = Cluster_cmap, save = T, plot_style = "violin_split", plot_name = "Senescence", dir = "Output/Fig3B-Malignancy/", wh = 1
)

Signature_Profile_S4(
    TCGA_PANCAN_nmf_Profiles, 
    Signature_name = "Stemness_mirandaCancerStemnessIntratumoral2019", 
    Group_name = "nmf_CD_predict", 
    plot = T, cmap = Cluster_cmap, save = T, plot_style = "violin_split", plot_name = "Stemness", dir = "Output/Fig3B-Malignancy/", wh = 1
)

Signature_Profile_S4(
    TCGA_PANCAN_nmf_Profiles, 
    Signature_name = "EMT_jainPancancerQuantitationEpithelialmesenchymal2022", 
    Group_name = "nmf_CD_predict", 
    plot = T, cmap = Cluster_cmap, save = T, plot_style = "violin_split", plot_name = "EMT", dir = "Output/Fig3B-Malignancy/", wh = 1
)

## Figure 3C
Group_dt <- merge(TCGA_PANCAN_nmf@Group$Clinical, TCGA_PANCAN_nmf@Group$nmf_CD_predict, by = "ID_sample")
Signature_dt <- TCGA_PANCAN_nmf@Signature$nmf_CD_h

Cox_list <- list()
for (CD in c(paste0("CD", 1:5))) {
    Cox_list[[CD]] <- Survival(
        Survival_dt = TCGA_PANCAN_nmf@Survival, Event = "OS", Right_sensor = 1825, mode = "multicox", 
        Group_name = "CD_subtype", Group_dt = Group_dt, Group_items = c("Cancer"),
        Signature_name = CD, Signature_dt = Signature_dt, Signature_items = c(CD),
        summary = T, plot = F, save = F
    )[CD, on = "Type"]
}
Cox_table <- rbindlist(Cox_list)
Cox_table$Type <- "nmf_CD_h"
Cox_table$HR_log2 <- log2(Cox_table$HR)
Cox_table$CI_upper_log2 <- log2(Cox_table$CI_upper)
Cox_table$CI_lower_log2 <- log2(Cox_table$CI_lower)

Plot <- ggplot(Cox_table) + 
    geom_hline(yintercept = 0, linetype = "dashed", linewidth = 1) +
    geom_linerange(mapping = aes(x = var, ymin = CI_lower_log2, ymax = CI_upper_log2, color = var), show.legend = F, linewidth = 1) +
    geom_point(mapping = aes(x = var, y = HR_log2, color = var), size = 5, show.legend = F) +
    scale_color_manual(values = CD_cmap) +
    scale_y_continuous(limits = c(-3.5, 3.5)) +
    xlab("") +
    ylab("Log2(Harzard ratio)") +
    theme_bw() +
    theme(
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.line = element_blank(),
        legend.position = "none",
        panel.border = element_rect(colour = "black", fill = NA, size = 1),
    ) +
    coord_flip()

## Figure 3D
for (cancer in c("HNSC", "LGG", "LUAD")) {
    TCGA_PANCAN_nmf_cancer <- Filter_sample(TCGA_PANCAN_nmf, info = T, filter = list(
        Filter_cancer = paste0("[disease_primary_disease == '", cancer, "']")
    ))
    Cluster_Survival <- Survival_S4(
        TCGA_PANCAN_nmf_cancer, 
        Event = "OS", Right_sensor = 1825, mode = "KM", 
        Group_name = "nmf_CD_predict",
        summary = F, plot = T, cmap = Cluster_cmap, 
        merge_plot = F, plotname = paste0(cancer, "_OS"), save = T, dir = dir
    )
    RMST_list <- list()
    for (Cluster in paste0("Cluster_", 1:5)) {
        Group_dt <- TCGA_PANCAN_nmf_cancer@Group$nmf_CD_predict[, .(
            ID_sample, 
            CD_subtype = fifelse(CD_subtype == Cluster, Cluster, "Others")
        )]
        Group_dt$CD_subtype <- factor(Group_dt$CD_subtype, levels = c("Others", Cluster))
        Group_dt <- merge(TCGA_PANCAN_nmf_cancer@Group$Clinical, Group_dt, by = "ID_sample")
        data <- na.omit(merge(TCGA_PANCAN_nmf_cancer@Survival, Group_dt, by = "ID_sample", all = T)[, .(ID_sample, OS, OS_time, CD_subtype)])
        fit <- rmst2(
            time = data$OS_time,
            status = data$OS,
            arm = ifelse(data$CD_subtype == Cluster, 1, 0),
            tau = 1825
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
}
