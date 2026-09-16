## Figure 6A-C
for (Cohort in c("TCGA", "CPTAC", "ICI")) {
    Prediction <- get(paste0("Prediction_", Cohort))
    Prediction$ID <- as.character(1:nrow(Prediction))
    Prediction$True_Label <- factor(as.character(Prediction$True_Label))
    Prediction$Predicted_Label <- factor(as.character(Prediction$Predicted_Label))
    Prediction[, Max_Probability := do.call(pmax, .SD), .SDcols = patterns("^Probability_Cluster_")]
    Prediction <- Prediction[order(True_Label, -Max_Probability)]
    Prediction$ID <- factor(Prediction$ID, levels = Prediction$ID)
    sourcedata <- Prediction[, .(
        `Sample ID` = ID_sample,
        `Actual Subtype` = fcase(
            True_Label == "Cluster_1", "PCD-P",
            True_Label == "Cluster_2", "PCD-PN",
            True_Label == "Cluster_3", "PCD-NF",
            True_Label == "Cluster_4", "PCD-FA",
            True_Label == "Cluster_5", "PCD-A",
            default = NA
        ),
        `Predicted Subtype` = fcase(
            Predicted_Label == "Cluster_1", "PCD-P",
            Predicted_Label == "Cluster_2", "PCD-PN",
            Predicted_Label == "Cluster_3", "PCD-NF",
            Predicted_Label == "Cluster_4", "PCD-FA",
            Predicted_Label == "Cluster_5", "PCD-A",
            default = NA
        ),
        `Prediction Probability of PCD-P` = Probability_Cluster_1,
        `Prediction Probability of PCD-PN` = Probability_Cluster_2,
        `Prediction Probability of PCD-NF` = Probability_Cluster_3,
        `Prediction Probability of PCD-FA` = Probability_Cluster_4,
        `Prediction Probability of PCD-A` = Probability_Cluster_5
    )]
    save2(data = sourcedata, dir = paste0("Output/", Cohort, "/Concordance/"), mode = "dw", format = "tsv")
    Probability <- ggplot(Prediction, aes(x = ID, y = Max_Probability, fill = True_Label)) +
        geom_bar(stat = "identity") +
        scale_fill_manual(values = Cluster_cmap) +
        theme_pubr() +
        theme(
            axis.text.y = element_blank(),
            axis.ticks.y = element_blank(),
            axis.title.y = element_blank(),
            legend.position = "none"
        ) +
        coord_flip()
    save2(data = Probability, dir = paste0("Output/", Cohort, "/Concordance/"), mode = "w", format = "img", wh = 1/8, scale = 1)
    data <- Prediction[, .(ID, True_Label, Predicted_Label)] 
    data$Freq <- 1
    Concordance <- ggplot(data = data, aes_string(axis1 = "True_Label", axis2 = "Predicted_Label", y = "Freq")) +
        geom_alluvium(aes_string(fill = "True_Label", order = "ID"), width = 0.1) +
        geom_stratum(aes_string(fill = "True_Label"), width = 0.1, color = "black") +
        # geom_label(stat = "stratum", aes(label = after_stat(stratum))) +
        scale_fill_manual(values = Cluster_cmap) +
        scale_x_discrete(expand = c(0.1, 0.1)) +
        theme_pubr() +
        theme(
            # axis.text.x = element_blank(),
            axis.text.y = element_blank(),
            axis.ticks.x = element_blank(),
            axis.ticks.y = element_blank(),
            axis.title.x = element_blank(),
            axis.title.y = element_blank(),
            axis.line.x = element_blank(),
            axis.line.y = element_blank(),
            legend.position = "none"
        )
    save2(data = Concordance, dir = paste0("Output/", Cohort, "/Concordance/"), mode = "w", format = "img", wh = 1/2.5, scale = 1)
}

## Figure 6D
Summary <- Signature_Profile_S4(
    TCGA_PANCAN_nmf, 
    Signature_name = "Signature", 
    Group_name = "nmf_CD_predict", 
    plot = F, summary = T, save = T, dir = "Output/Fig6D-Summary/", plot_name = "Signature"
)
Summary_dcast <- dcast(Summary, CD_subtype ~ Signature_name, value.var = "average")
Summary_dcast <- Summary_dcast[, .( CD_subtype, Proliferation, Stemness, Senescence, JAK_STAT, WNT, MAPK, PI3K, NFkB, TGFb, EMT, p53, EGFR)] %>% column_to_rownames("CD_subtype")

Summary_dcast <- rbind(
    max = rep(1, ncol(Summary_dcast)),
    min = rep(-1, ncol(Summary_dcast)),
    Summary_dcast
)

for (Cluster in paste0("Cluster_", 1:5)) {
    data <- Summary_dcast[c("max", "min", Cluster), ]
    wh <- 1
    scale <- 0.75
    path_png <- paste0("Output/Fig6D-Summary/", Cluster, "_radar.png")
    path_pdf <- paste0("Output/Fig6D-Summary/", Cluster, "_radar.pdf")
    pdf(path_pdf, width = 6*wh*scale, height = 6*scale)
        radarchart(
            data, axistype = 1,
            pcol = Cluster_cmap[[Cluster]], pfcol = scales::alpha(Cluster_cmap[[Cluster]], 0.5), plwd = 2, plty = 1,
            cglcol = "grey", cglty = 1, cglwd = 0.8,
            axislabcol = "grey", 
            vlcex = 0.7, vlabels = colnames(data),
            caxislabels = c(-1, -0.5, 0, 0.5, 1)
        )
    dev.off()
}

