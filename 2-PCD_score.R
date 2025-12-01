# 2. PCD Scores
dir.create("Output/2-PCD_score")
pacman::p_load(data.table, ggplot2, ggpubr, GSVA)
cmap <- setNames(
    c("#fdcc83", "#f6afb0", "#71b3a8", "#6292b0"), 
    c("Pyroptosis", "Necroptosis", "Ferroptosis", "Apoptosis")
)
PCD_Gene <- fread("Utility/PCD_Gene.tsv")
PCD_Expr <- fread("Data/CPTAC/PCD_Expr.tsv")
PCD_Subtype <- fread("Output/1-Infer_PCD_subtype/Prediction.tsv")[, .(ID_sample, PCD_subtype = Subtype)]
# PCD_Subtype <- fread("Data/CPTAC/PCD_Subtype_NMF.tsv") # NMF-derived subtype

## Data pre-processing
PCD_Expr <- PCD_Expr[
    , lapply(.SD, function(x) (x - min(x))/(max(x) - min(x))), .SDcols = -c("ID_sample")][
        , ID_sample := PCD_Expr$ID_sample
]
PCD_Expr <- transpose(PCD_Expr, keep.names = "Gene", make.names = "ID_sample") %>% 
    column_to_rownames("Gene") %>% 
    as.matrix() 
PCD_Gene_list <- list()
for (name in unique(PCD_Gene[["PCD"]])) {
    PCD_Gene_list[[name]] <- intersect(PCD_Gene[name, on = "PCD"]$Gene, rownames(PCD_Expr))
}

## Calculate ssGSEA score
PCD_Score <- scale(t(gsva(ssgseaParam(PCD_Expr, PCD_Gene_list)))) %>% 
    as.data.frame() %>% 
    rownames_to_column("ID_sample") %>% 
    as.data.table()        
PCD_Score <- merge(PCD_Subtype, PCD_Score, all.x = T, by = "ID_sample")
fwrite(PCD_Score, "Output/2-PCD_score/PCD_Score.tsv", sep = "\t")

## Plot
for (Subtype in c("CD-P", "CD-PN", "CD-NF", "CD-FA", "CD-A")) {
    PCD_Score_Subtype <- PCD_Score[Subtype, on = "PCD_subtype"][, PCD_subtype := NULL] %>% 
        melt(., variable.name = "PCD", value.name = "ssGSEA_Score")
    PCD_Score_Subtype$PCD <- factor(PCD_Score_Subtype$PCD, levels = c("Pyroptosis", "Necroptosis", "Ferroptosis", "Apoptosis"))
    Plot <- ggplot(PCD_Score_Subtype, mapping = aes(x = PCD, y = ssGSEA_Score, fill = PCD)) + 
        stat_boxplot(mapping = aes(x = PCD, y = ssGSEA_Score, fill = PCD), width = 0.5, position = position_dodge(0.5)) +   
        geom_boxplot(aes(fill = PCD), position = position_dodge(0.5), width = 0.5) +
        scale_fill_manual(values = cmap) +
        scale_color_manual(values = cmap) +
        stat_compare_means(aes(group = PCD), ref.group = ".all.", method = "wilcox.test", label = "p.signif", show.legend = F, hide.ns = T) +
        geom_hline(yintercept = 0, linetype = "dashed", color = "gray", size = 1) +
        theme(plot.title = element_text(hjust = 0.5)) +
        theme_pubr() +
        labs(color = "") +
        scale_y_continuous(limits = c(-4, 4), breaks = c(-2, 0, 2))   
    png(paste0("Output/2-PCD_score/", Subtype, "_Score.png"), width = 1800, height = 1800, res = 300)
    print(Plot)
    dev.off()
}
