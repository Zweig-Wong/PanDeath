pacman::p_load(ComplexHeatmap, circlize, ggpubr, maftools)

## Figure 2A
maf <- read.maf("TCGA_PANCAN/Raw/mc3.v0.2.8.PUBLIC.maf.gz")
maf_summary <- as.data.table(getSampleSummary(maf))[, total := NULL]
maf_summary[, ID_sample := str_sub(Tumor_Sample_Barcode, 1, 15)]
oncostrip(maf = maf, genes = CD_gene, writeMatrix = TRUE)
CD_gene_mutation <- fread("Output/Fig2A-Oncoplot/onco_matrix.txt")
setnames(CD_gene_mutation, "V1", "Gene")
setnames(CD_gene_mutation, as.character(maf_summary$Tumor_Sample_Barcode), maf_summary$ID_sample)

Sample <- intersect(TCGA_PANCAN_nmf@Clinical$ID_sample, colnames(CD_gene_mutation))
CD_gene_mutation <- CD_gene_mutation[, c("Gene", ..Sample)]
CD_gene_mutation_mat <- CD_gene_mutation %>% column_to_rownames("Gene") %>% as.matrix()
CD_gene_mutation_mat[CD_gene_mutation_mat == "0"] <- ""

Mutaion_annotation_list <- list()
for (col in colnames(CD_gene_mutation)[-1]) {
    Mutaion_annotation_list[[col]] <- unique(CD_gene_mutation[[col]])
}
Mutaion_annotation <- unique(unlist(Mutaion_annotation_list))
Mutaion_annotation <- as.factor(setdiff(Mutaion_annotation, c("", "0")))
cmap <- setNames(ColorCard$Article$NatMed_1[1:length(Mutaion_annotation)], Mutaion_annotation)
alter_fun <- list(
    background = function(...) NULL,
    Multi_Hit = alter_graphic("rect", fill = cmap[["Multi_Hit"]]),
    Frame_Shift_Del = alter_graphic("rect", fill = cmap[["Frame_Shift_Del"]]),
    Frame_Shift_Ins = alter_graphic("rect", fill = cmap[["Frame_Shift_Ins"]]),
    In_Frame_Del = alter_graphic("rect", fill = cmap[["In_Frame_Del"]]),
    In_Frame_Ins = alter_graphic("rect", fill = cmap[["In_Frame_Ins"]]),
    Missense_Mutation = alter_graphic("rect", fill = cmap[["Missense_Mutation"]]),
    Nonsense_Mutation = alter_graphic("rect", fill = cmap[["Nonsense_Mutation"]]),
    Nonstop_Mutation = alter_graphic("rect", fill = cmap[["Nonstop_Mutation"]]),
    Splice_Site = alter_graphic("rect", fill = cmap[["Splice_Site"]]),
    Translation_Start_Site = alter_graphic("rect", fill = cmap[["Translation_Start_Site"]])
)

Annotation_column_df <- Reduce(function(x, y) merge(x, y, all.x = TRUE, by = "ID_sample"), list(
    Get_group(TCGA_PANCAN_nmf, "nmf_CD_predict", "CD_subtype"),
    TCGA_PANCAN_CNV$Ploidy[, .(ID_sample, WGD, ITH)],
    TCGA_PANCAN_CNV$LOH[, .(ID_sample, LOH = Fraction_of_segs_with_LOH)],
    TCGA_PANCAN_CNV$HRD[, .(ID_sample, HRD)]
))[order(CD_subtype)][Sample, on = "ID_sample"]
Annotation_column_df <- Annotation_column_df %>% column_to_rownames("ID_sample")

Annotation_column_top <- HeatmapAnnotation(
    cbar = anno_oncoprint_barplot(show_fraction = TRUE),
    df = Annotation_column_df[c("CD_subtype")],
    col = list(
        CD_subtype = Cluster_cmap
    ),        
    na_col = "#ffffff"
)
Annotation_column_bottom <- HeatmapAnnotation(
    cbar = anno_oncoprint_barplot(show_fraction = TRUE),
    df = Annotation_column_df[c("WGD", "LOH", "HRD", "ITH")],
    col = list(
        WGD = colorRamp2(c(0, 2), c("white", "#c38a8b")),
        LOH = colorRamp2(c(0, 1), c("white", "#935ea3")),
        HRD = colorRamp2(c(0, 100), c("white", "#45677c")),
        ITH = colorRamp2(c(0, 1), c("white", "#508078"))
    ),
    na_col = "#ffffff"
)
CD_gene_mutation_mat <- CD_gene_mutation_mat[, rownames(Annotation_column_df)]
Column_split <- Annotation_column_df[["CD_subtype"]]

Annotation_row_df <- CD_Merged[CD_gene_mutation$Gene[1:25], on = "Gene"][-c(2, 22, 23)] %>% 
    column_to_rownames("Gene")

Annotation_row <- rowAnnotation(
    cbar = anno_oncoprint_barplot(show_fraction = TRUE),
    df = Annotation_row_df,
    col = list(
        cell_death = c(
            "Pyroptosis" = "#fdcc83",
            "Necroptosis" = "#f6afb0",
            "Ferroptosis" = "#71b3a8",
            "Apoptosis" = "#6292b0"
        )
    ),
    na_col = "#ffffff"
)
CD_gene_mutation_mat <- CD_gene_mutation_mat[rownames(Annotation_row_df), ]
Row_split <- Annotation_row_df$CD_subtype

CD_gene_mutation_oncoprint <- oncoPrint(
    CD_gene_mutation_mat, remove_empty_columns = TRUE, 
    alter_fun = alter_fun, col = cmap,
    pct_side = "right", row_names_side = "left",
    top_annotation = Annotation_column_top, bottom_annotation = Annotation_column_bottom,
    right_annotation = Annotation_row,
    column_split = Column_split
)

TopGene <- CD_gene_mutation$Gene[1:25]
Sample <- intersect(TCGA_PANCAN_nmf@Clinical$ID_sample, colnames(CD_gene_mutation))
CD_gene_mutation <- CD_gene_mutation[, c("Gene", ..Sample)]
CD_gene_mutation <- transpose(CD_gene_mutation, keep.names = "ID_sample", make.names = "Gene")
CD_gene_mutation <- merge(CD_gene_mutation, TCGA_PANCAN_nmf@Group$nmf_CD_predict, all.y = T, by = "ID_sample")

CD_gene_mutation_01 <- CD_gene_mutation[
    , lapply(.SD, function(x) fifelse(x == 0 | x == "", 0, 1)), .SDcols = -c("ID_sample", "CD_subtype")][ 
        , CD_subtype := CD_gene_mutation$CD_subtype]

CD_gene_mutation_all <- transpose(
    CD_gene_mutation_01[
        , c(lapply(.SD, sum, na.rm = TRUE), N = .N), .SDcols = -c("CD_subtype")][
            , lapply(.SD, function(x) x/N)][
                , N:= NULL][
                    , CD_subtype := "ALL"]
    , keep.names = "Gene", make.names = "CD_subtype"
)

CD_gene_mutation_subtype <- transpose(
    CD_gene_mutation_01[
        , c(lapply(.SD, sum, na.rm = TRUE), .N), by = "CD_subtype"][
            , lapply(.SD, function(x) x/N), by = "CD_subtype"][
                , N:= NULL][
                    order(CD_subtype)]
    , keep.names = "Gene", make.names = "CD_subtype"

)

CD_gene_mutation_diff <- merge(CD_gene_mutation_subtype, CD_gene_mutation_all, by = "Gene")
CD_gene_mutation_diff <- CD_gene_mutation_diff[, lapply(.SD, function(x) x - ALL), .SDcols = -c("Gene")][, ALL := NULL][, Gene := CD_gene_mutation_diff$Gene]
setcolorder(CD_gene_mutation_diff, "Gene")

CD_gene_mutation_mat <- CD_gene_mutation_diff %>% column_to_rownames("Gene") %>% as.matrix()
CD_gene_mutation_mat <- CD_gene_mutation_mat*100

Annotation_row_df <- CD_Merged[
c(
    "TP53", "HUWE1", "BIRC6", "CREBBP", "CTNNB1", "IGF2R", "SPTAN1",
    "NLRP3", "ERBB3", "NFE2L2", "TLR4", "ERBB2", "CASP8", "BAP1", "BRCA1", 
    "NTRK3", "ITGB4", "ROCK1", "DPYD", "HGF", "NLRP9", "NLRP2", "NLRP1", "DHX9", "MADD"
), on = "Gene"][-c(2, 14, 15)] %>% column_to_rownames("Gene")

Heatmap_mutation <- Plot_Heatmap(
    CD_gene_mutation_mat, Annotation_cell_value = T,
    Annotation_cell_color = colorRamp2(seq(2, -2, length.out = 5), ColorCard$Contrast$PG[2:6]),
    Row_name = T, 
    Annotation_row_df = Annotation_row_df, 
    Annotation_row_color = list(cell_death = c("Pyroptosis" = "#fdcc83", "Necroptosis" = "#f6afb0", "Ferroptosis" = "#71b3a8", "Apoptosis" = "#6292b0")), Column_name = T, Row_cluster = F,
    plot_name = "CD_gene_mutation_summary", save = T, dir = "Output/Fig2A-Oncoplot/", wh = 0.31, scale = 2.5
)    


## Figure 2B
Gene_CNV <- fread("TCGA_PANCAN/Raw/ISAR_GISTIC.all_thresholded.by_genes.txt.gz")
setnames(Gene_CNV, substr(colnames(Gene_CNV), 1, 15))
column <- c("Gene Symbol", "Locus ID", "Cytoband", intersect(TCGA_PANCAN_nmf@Clinical$ID_sample, colnames(Gene_CNV)))
Gene_CNV <- Gene_CNV[, ..column]

Sample <- intersect(TCGA_PANCAN_nmf@Clinical$ID_sample, colnames(Gene_CNV))
CNV <- Gene_CNV[, c("Gene Symbol", ..Sample)]
CNV <- transpose(CNV[CD_gene, on = "Gene Symbol"], keep.names = "ID_sample", make.names = "Gene Symbol")
CNV <- merge(CNV, TCGA_PANCAN_nmf@Group$nmf_CD_predict, all.y = T, by = "ID_sample")
CNV_amp <- CNV[
    , lapply(.SD, function(x) fifelse(x > 1, 1, 0)), .SDcols = -c("ID_sample", "CD_subtype")][ 
        , CD_subtype := CNV$CD_subtype]
CNV_amp_all <- transpose(
    CNV_amp[
        , c(lapply(.SD, sum, na.rm = TRUE), N = .N), .SDcols = -c("CD_subtype")][
            , lapply(.SD, function(x) x/N)][
                , N:= NULL][
                    , CD_subtype := "ALL"]
    , keep.names = "Gene", make.names = "CD_subtype"
)
CNV_amp_subtype <- transpose(
    CNV_amp[
        , c(lapply(.SD, sum, na.rm = TRUE), .N), by = "CD_subtype"][
            , lapply(.SD, function(x) x/N), by = "CD_subtype"][
                , N:= NULL][
                    order(CD_subtype)]
    , keep.names = "Gene", make.names = "CD_subtype"
)

CNV_del <- CNV[
    , lapply(.SD, function(x) fifelse(x < -1, 1, 0)), .SDcols = -c("ID_sample", "CD_subtype")][ 
        , CD_subtype := CNV$CD_subtype]
CNV_del_all <- transpose(
    CNV_del[
        , c(lapply(.SD, sum, na.rm = TRUE), N = .N), .SDcols = -c("CD_subtype")][
            , lapply(.SD, function(x) x/N)][
                , N:= NULL][
                    , CD_subtype := "ALL"]
    , keep.names = "Gene", make.names = "CD_subtype"
)
CNV_del_subtype <- transpose(
    CNV_del[
        , c(lapply(.SD, sum, na.rm = TRUE), .N), by = "CD_subtype"][
            , lapply(.SD, function(x) x/N), by = "CD_subtype"][
                , N:= NULL][
                    order(CD_subtype)]
    , keep.names = "Gene", make.names = "CD_subtype"
)

Gene_list <- c("GSDMD", "FADD", "TFRC", "TGFB2", "CASP3", "BCL2", "TP53", "GPX4")
PiePlot_list <- list()
Data_list <- list()
for (Gene in Gene_list) {
    data <- transpose(rbindlist(list(
        AMP = CNV_amp_subtype[Gene, on = "Gene"], 
        DEL = CNV_del_subtype[Gene, on = "Gene"]
    ), idcol = "CNV")[, Gene := NULL], keep.names = "CD_subtype", make.names = "CNV")[, Blank := 1 - AMP - DEL]
    
    Data_list[[Gene]] <- list()
    PiePlot_list[[Gene]] <- list()
    for (Cluster in paste0("Cluster_", 1:5)) {
        data_cluster <- melt(data[Cluster, on = "CD_subtype"], id.vars = "CD_subtype")
        data_cluster$variable <- factor(data_cluster$variable, levels = c("Blank", "DEL", "AMP"))
        Data_list[[Gene]][[Cluster]] <- data_cluster
        PiePlot_list[[Gene]][[Cluster]] <- ggplot(data_cluster, aes(x = CD_subtype, y = value, fill = variable)) +
            geom_bar(stat = "identity", width = 1) +
            coord_polar(theta = "y") +
            scale_fill_manual(values = c("AMP" = "#f6aeaf", "DEL" = "#6192af", "Blank" = "#f5f6f6")) +
            labs(title = Gene) +
            geom_text(aes(label = scales::percent(value, accuracy = 0.01)), position = position_stack(vjust = 0.5)) +
            theme_pubr() +
            theme(
                axis.title.x = element_blank(),
                axis.text.x = element_blank(),
                axis.ticks.x = element_blank(),
                axis.title.y = element_blank(),
                axis.text.y = element_blank(),
                axis.ticks.y = element_blank(),
                axis.line = element_blank(),
                legend.position = "none",
                plot.title = element_text(hjust = 0.5)
            )
    }
}

plot_list <- lapply(paste0("Cluster_", 1:5), function(Cluster) {
    row_plots <- lapply(Gene_list, function(Gene) {
        PiePlot_list[[Gene]][[Cluster]]
    })
    ggarrange(plotlist = row_plots, ncol = length(row_plots), nrow = 1)
})
CNV_Pie <- ggarrange(plotlist = plot_list, ncol = 1, nrow = length(plot_list))


## Figure 2C
load("TCGA_PANCAN/Processed/TCGA_PANCAN_SNV.R")
mutation_name <- c("total_perMB_log", "total_SNV_perMB_log", "Silent_Mutation_perMB_log", "Nonsilent_Mutation_perMB_log")
Mutation_load_log <- merge(
    TCGA_PANCAN_SNV[, c("ID_sample", ..mutation_name)],
    TCGA_PANCAN_nmf@Group$nmf_CD_predict,
    all.y = TRUE, by = "ID_sample"
)

SNV_high_list <- list()
for (mutation in mutation_name) {
    Mutation_load_log[, paste0(mutation, "_high") := ifelse(get(mutation) > 1, "YES", "NO")]
    SNV_prop <- Proportion(
        Group_dt.x = Mutation_load_log, Group_item.x = "CD_subtype",
        Group_dt.y = Mutation_load_log, Group_item.y = paste0(mutation, "_high"), 
        summary = T, save = F
    )
    SNV_prop <- SNV_prop["YES", on = paste0(mutation, "_high")]
    SNV_prop[[paste0(mutation, "_high")]] <- NULL
    SNV_high_list[[mutation]] <- SNV_prop
}
SNV_high <- rbindlist(SNV_high_list, idcol = "SNV")
SNV_high$SNV <- factor(SNV_high$SNV, levels = mutation_name)

SNV_high_plot <- ggplot(data = SNV_high, aes(x = CD_subtype, y = Proportion, fill = CD_subtype)) +
    facet_wrap("SNV", nrow = 1) +
    geom_bar(stat = "identity", width = 0.8) +
    geom_text(mapping = aes(x = CD_subtype, y = Proportion + 0.002, label = percent(Proportion, 0.1)), size = 5, show.legend = F) +
    scale_fill_manual(values = Cluster_cmap) +
    scale_y_continuous(expand = c(0, 0)) +
    theme_pubr() +
    theme(
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.title.y = element_blank(),
        legend.position = "none",
        strip.background = element_rect(color = "white", fill = "white", linetype = "blank"),
        panel.spacing = unit(3, "lines") 
    )

## Figure 2D
Signature_names <- c("Ploidy", "Aneuploidy", "Segments", "LOH", "HRD", "MSI")
TCGA_PANCAN_nmf <- Add_signature(TCGA_PANCAN_nmf, Signature_dt_list = TCGA_PANCAN_CNV[Signature_names])
Signature_Landscape(
    object = TCGA_PANCAN_nmf, 
    Signature_names = Signature_names,
    Group_name = "nmf_CD_predict", Group_item = "CD_subtype", 
    test_method = "t.test", 
    dir = "Output/Fig2D-Genomic Stress/Genomic/", plot_name = "Genomic",
    wh = 0.75,
    scale = 2
)

## Figure 2E
load("TCGA_PANCAN/Processed/TCGA_PANCAN_SV.R")
TCGA_PANCAN_SV <- TCGA_PANCAN_SV[, -c("ID_patient")]
save2(data = TCGA_PANCAN_SV, dir = "Output/Fig2E-Genomic Stress/SV/", mode = "dw", format = "tsv")
load2(SP = "260311-Fig2-Stress", path = "Output/Fig2E-Genomic Stress/SV/TCGA_PANCAN_SV.R")
TCGA_PANCAN_nmf <- Add_signature(TCGA_PANCAN_nmf, Signature_dt = TCGA_PANCAN_SV)
Signature_Profile_S4(
    TCGA_PANCAN_nmf, 
    Signature_name = "TCGA_PANCAN_SV", 
    Group_name = "nmf_CD_predict", 
    plot = T, cmap = Cluster_cmap, save = T, plot_style = "violin_split", plot_name = "SV", wh = 1, dir = "Output/Fig2E-Genomic Stress/SV/"
)
Signature_Profile_S4(
    TCGA_PANCAN_nmf, 
    Signature_name = "TCGA_PANCAN_SV", 
    Group_name = "nmf_CD_predict", 
    plot = T, cmap = Cluster_cmap, save = T, plot_style = "violin_split", plot_name = "SV_cartesian", wh = 1, dir = "Output/Fig2E-Genomic Stress/SV/"
)


## Figure 2F
load2(SP = "Marker", info = F, path = "DDR")
DDR_marker <- Get_marker(DDR, subset = list(Pathway = "[Category == 'Pathway']"))[, Term := "DDR"]
save2(data = DDR_marker, dir = "Output/Fig2F-Genomic Stress/DDR/", mode = "dw", format = "tsv")
DDR_score <- Calculate_Score(
    Marker = DDR_marker, 
    Expr = TCGA_PANCAN_RNA@Signature$RNA_std, 
    method = "ssGSEA", save = T, dir = "Output/Fig2F-Genomic Stress/DDR/"
)
load2(SP = "260311-Fig2-Stress", path = "Output/Fig2F-Genomic Stress/DDR/Score.R")
DDR_score <- Score
TCGA_PANCAN_nmf <- Add_signature(TCGA_PANCAN_nmf, Signature_dt = DDR_score)
Signature_Profile_S4(
    TCGA_PANCAN_nmf, 
    Signature_name = "DDR_score", 
    Group_name = "nmf_CD_predict", 
    plot = T, cmap = Cluster_cmap, save = T, plot_style = "violin_split", plot_name = "DDR", wh = 1, dir = "Output/Fig2F-Genomic Stress/DDR/"
)

## Figure 2G
load2(SP = "Marker", info = F, path = "Oxidative_Stress")
Oxidative_Stress_marker <- Get_marker(Oxidative_Stress, subset = list(Oxidative = "[Term == 'Oxidative']"))
save2(data = Oxidative_Stress_marker, dir = "Output/Fig2G-Metabolic Stress/Oxidative_Stress/", mode = "dw", format = "tsv")
Oxidative_Stress_score <- Calculate_Score(
    Marker = Oxidative_Stress_marker, 
    Expr = TCGA_PANCAN_RNA@Signature$RNA_std, 
    method = "ssGSEA", save = T, dir = "Output/Fig2G-Metabolic Stress/Oxidative_Stress/"
)
load2(SP = "260311-Fig2-Stress", path = "Output/Fig2G-Metabolic Stress/Oxidative_Stress/Score.R")
Oxidative_Stress_score <- Score
TCGA_PANCAN_nmf <- Add_signature(TCGA_PANCAN_nmf, Signature_dt = Oxidative_Stress_score)
Signature_Profile_S4(
    TCGA_PANCAN_nmf, 
    Signature_name = "Oxidative_Stress_score", 
    Group_name = "nmf_CD_predict", 
    plot = T, cmap = Cluster_cmap, save = T, plot_style = "violin_split", plot_name = "Oxidative_Stress", wh = 1, dir = "Output/Fig2G-Metabolic Stress/Oxidative_Stress/"
)

## Figure 2H
load2(SP = "Marker", info = F, path = "Hypoxic_Stress")
Hypoxic_Stress_marker <- Get_marker(Hypoxic_Stress, subset = list(Hypoxia = "[Term == 'Hypoxia']"))
save2(data = Hypoxic_Stress_marker, dir = "Output/Fig2H-Metabolic Stress/Hypoxic_Stress/", mode = "dw", format = "tsv")
Hypoxic_Stress_score <- Calculate_Score(
    Marker = Hypoxic_Stress_marker, 
    Expr = TCGA_PANCAN_RNA@Signature$RNA_std, 
    method = "ssGSEA", save = T, dir = "Output/Fig2H-Metabolic Stress/Hypoxic_Stress/"
)
load2(SP = "260311-Fig2-Stress", path = "Output/Fig2H-Metabolic Stress/Hypoxic_Stress/Score.R")
Hypoxic_Stress_score <- Score
TCGA_PANCAN_nmf <- Add_signature(TCGA_PANCAN_nmf, Signature_dt = Hypoxic_Stress_score)
Signature_Profile_S4(
    TCGA_PANCAN_nmf, 
    Signature_name = "Hypoxic_Stress_score", 
    Group_name = "nmf_CD_predict", 
    plot = T, cmap = Cluster_cmap, save = T, plot_style = "violin_split", plot_name = "Hypoxic_Stress", wh = 1, dir = "Output/Fig2H-Metabolic Stress/Hypoxic_Stress/"
)

## Figure 2I
load2(SP = "Marker", info = F, path = "Autophagy")
Autophagy_marker <- Get_marker(Autophagy, subset = list(Autophagy = "[Term == 'Autophagy' & Source == 'REACTOME', ]"))
save2(data = Autophagy_marker, dir = "Output/Fig2I-Metabolic Stress/Autophagy/", mode = "dw", format = "tsv")
Autophagy_score <- Calculate_Score(
    Marker = Autophagy_marker, 
    Expr = TCGA_PANCAN_RNA@Signature$RNA_std, 
    method = "ssGSEA", save = T, dir = "Output/Fig2I-Metabolic Stress/Autophagy/"
)
load2(SP = "260311-Fig2-Stress", path = "Output/Fig2I-Metabolic Stress/Autophagy/Score.R")
Autophagy_score <- Score
TCGA_PANCAN_nmf <- Add_signature(TCGA_PANCAN_nmf, Signature_dt = Autophagy_score)
Signature_Profile_S4(
    TCGA_PANCAN_nmf, 
    Signature_name = "Autophagy_score", 
    Group_name = "nmf_CD_predict", 
    plot = T, cmap = Cluster_cmap, save = T, plot_style = "violin_split", plot_name = "Autophagy", wh = 1, dir = "Output/Fig2I-Metabolic Stress/Autophagy/"
)

## Figure 2J
load2(SP = "Marker", info = F, path = "UPR")
UPR_marker <- Get_marker(UPR)[, Term := "UPR"]
save2(data = UPR_marker, dir = "Output/Fig2J-Metabolic Stress/UPR/", mode = "dw", format = "tsv")
UPR_score <- Calculate_Score(
    Marker = UPR_marker, 
    Expr = TCGA_PANCAN_RNA@Signature$RNA_std, 
    method = "ssGSEA", save = T, dir = "Output/Fig2J-Metabolic Stress/UPR/"
)
load2(SP = "260311-Fig2-Stress", path = "Output/Fig2J-Metabolic Stress/UPR/Score.R")
UPR_score <- Score
TCGA_PANCAN_nmf <- Add_signature(TCGA_PANCAN_nmf, Signature_dt = UPR_score)
Signature_Profile_S4(
    TCGA_PANCAN_nmf, 
    Signature_name = "UPR_score", 
    Group_name = "nmf_CD_predict", 
    plot = T, cmap = Cluster_cmap, save = T, plot_style = "violin_split", plot_name = "UPR", wh = 1, dir = "Output/Fig2J-Metabolic Stress/UPR/"
)