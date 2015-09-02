##----------------------------------
## required libraries
##----------------------------------
library(topGO)
library(org.Hs.eg.db)
library(Rgraphviz)

library(dplyr)
library(tidyr)
library(data.table)
library(foreach)
library(stringr)
library(ggplot2)
library(RColorBrewer)

##----------------------------------
## functions
##----------------------------------
find.dm.events <- function(DT) {
    if("gene_symbol" %in% colnames(DT)){
        DT %>%
        filter(abs(delta_psi_mean) >= 0.05,
               Control_n / max(Control_n, na.rm = TRUE)  >= 0.75,
               DM1_n / max(DM1_n, na.rm = TRUE) >= 0.75,
               DM1_n_sig / DM1_n >= 0.25) %>%
                   select(gene_symbol, event_name, isoforms, Control_psi_mean, Control_psi_sd, Control_n,
                          DM1_psi_mean, DM1_psi_sd, DM1_n, delta_psi_mean, DM1_n_sig) %>%
                              arrange(desc(abs(delta_psi_mean)))
    } else {
        DT %>%
        filter(abs(delta_psi_mean) >= 0.05,
               Control_n / max(Control_n, na.rm = TRUE)  >= 0.75,
               DM1_n / max(DM1_n, na.rm = TRUE) >= 0.75,
               DM1_n_sig / DM1_n >= 0.25) %>%
                   select(event_name, isoforms, Control_psi_mean, Control_psi_sd, Control_n,
                          DM1_psi_mean, DM1_psi_sd, DM1_n, delta_psi_mean, DM1_n_sig) %>%
                              arrange(desc(abs(delta_psi_mean)))
    }
}

load("~/Projects/DMseq/bin/sum.InformativeMisoCounts.R")

##----------------------------------
## comparison data
##----------------------------------
event_type <- "nonUTRevents.multi"

allControls_res <- tbl_dt(fread(paste("~/Projects/DMseq/results/allControls/allControls", event_type, "results.txt", sep = "_")))
allControls_res <- allControls_res %>% mutate(delta_psi = Quad_psi_mean - Tibialis_psi_mean)

quad_vs_tibialis <- allControls_res %>%
    select(gene_symbol, event_name, Quad_psi_mean, Quad_n, Tibialis_psi_mean, Tibialis_n, Quad_vs_Tibialis_n_sig, isoforms, delta_psi) %>%
        filter(abs(delta_psi) >= 0.05,
               Quad_n / max(Quad_n, na.rm = TRUE) >= 0.75,
               Tibialis_n / max(Tibialis_n, na.rm = TRUE) >= 0.75,
               Quad_vs_Tibialis_n_sig / Quad_n >= 0.25) %>%
                   arrange(desc(Quad_vs_Tibialis_n_sig))

## compare quad and tibialis psi values from healthy indivudals
allControls_res %>%
    filter(Quad_n / max(Quad_n, na.rm = TRUE) >= 0.75,
           Tibialis_n / max(Tibialis_n, na.rm = TRUE) >= 0.75) %>%
        group_by(isoforms) %>%
            summarize(delta_psi = Quad_psi_mean - Tibialis_psi_mean, mean_psi = (Quad_psi_mean + Tibialis_psi_mean)/2) %>%
                ggplot(aes(x = mean_psi, y = delta_psi,
                           color = isoforms %in% quad_vs_tibialis$isoforms,
                           alpha = isoforms %in% quad_vs_tibialis$isoforms)) +
                    geom_point() +
                        scale_alpha_discrete(range = c(0.25, 0.9)) +
                            guides(col = guide_legend(title = "Differentially Regulated"),
                                   alpha = FALSE) +
                                labs(x = "Mean PSI",
                                     y = "Delta PSI \n Healthy Quadricep vs Healthy Tibialis") +
                                    theme_bw(20)


allControls_res %>%
    filter(Quad_n / max(Quad_n, na.rm = TRUE) >= 0.75,
           Tibialis_n / max(Tibialis_n, na.rm = TRUE) >= 0.75) %>%
        group_by(isoforms) %>%
            summarize(delta_psi = Quad_psi_mean - Tibialis_psi_mean, mean_psi = (Quad_psi_mean + Tibialis_psi_mean)/2) %>%
                ggplot(aes(x = delta_psi, color = isoforms %in% quad_vs_tibialis$isoforms)) +
                    stat_ecdf() +
                        labs(x = "CDF",
                             y = "Delta PSI \n Healthy Quadricep vs Healthy Tibialis") +
                            theme_bw(20)
##
## Load DM vs control data for each tissue
##
tibialis_pdata <- tbl_dt(fread("~/Projects/DMseq/data/DM_tibialis_pdata.txt"))
tibialis_res <- tbl_dt(fread(paste("~/Projects/DMseq/results/tibialis/tibialis", event_type, "results.txt", sep = "_")))

quadricep_pdata <- tbl_dt(fread("~/Projects/DMseq/data/DM_quadricep_pdata.txt"))
quadricep_res <- tbl_dt(fread(paste("~/Projects/DMseq/results/quadricep/quadricep", event_type, "results.txt", sep = "_")))

dm_tibialis <- find.dm.events(tibialis_res)
dm_quadricep <- find.dm.events(quadricep_res)

## plot all dm vs control delta psi
setnames(tibialis_res, c(paste("tibialis_", names(tibialis_res)[1:17], sep = ""), "isoforms"))
setnames(quadricep_res, c(paste("quad_", names(quadricep_res)[1:17], sep = ""), "isoforms"))
common_events <- inner_join(quadricep_res, tibialis_res, by = "isoforms")

common_events %>% 
    filter(tibialis_Control_n / max(tibialis_Control_n, na.rm = TRUE) >= 0.75,
           tibialis_DM1_n / max(tibialis_DM1_n, na.rm = TRUE) >= 0.75,
           quad_Control_n / max(quad_Control_n, na.rm = TRUE) >= 0.75,
           quad_DM1_n / max(quad_DM1_n, na.rm = TRUE) >= 0.75) %>%
        ggplot(aes(x=tibialis_delta_psi_mean, y=quad_delta_psi_mean,
                   alpha = isoforms %in% union(dm_tibialis$isoforms, dm_quadricep$isoforms),
                   color = isoforms %in% union(dm_tibialis$isoforms, dm_quadricep$isoforms))) +
            geom_point()+
                scale_alpha_discrete(range = c(0.25, 0.9)) +
                    guides(col = guide_legend(title = "Dysregulated in DM1"),
                           alpha = FALSE) +
                        labs(x = "Tibialis Delta PSI \n DM1 vs. Control",
                             y = "Quadricep Delta PSI \n DM1 vs. Control") +
                            theme_bw(20)

##----------------------------------
## load psi value data
##----------------------------------
## con_data <- tbl_dt(fread(paste("~/Projects/DMseq/data/allSamples", event_type, "consolidatedSummaries.txt", sep = "_")))
## pdata <- tbl_dt(fread("~/Projects/DMseq/data/DM_sample_pdata.txt"))
## setnames(pdata, c("sample", "diagnosis", "tissue", "group", "read_length"))

## con_data <- left_join(con_data, select(pdata, sample, tissue, diagnosis), by = "sample")
## con_data <- con_data %>% filter(tissue %in% c("Tibialis", "Quad"))
## con_data[, informativeCounts := sum.InformativeMisoCounts(con_data$counts)]
## f_con_data <- con_data %>% filter(informativeCounts >= 20)

##----------------------------------
## compare delta psi values
##----------------------------------
dm_dysreg_events_union <- union(dm_tibialis$isoforms, dm_quadricep$isoforms)
dm_dysreg_events_intersect <- intersect(dm_tibialis$isoforms, dm_quadricep$isoforms)

event_set1 <- intersect(dm_dysreg_events_intersect, quad_vs_tibialis$isoforms)
event_set2 <- intersect(dm_dysreg_events_union, quad_vs_tibialis$isoforms)
event_set3 <- union(dm_dysreg_events_union, quad_vs_tibialis$isoforms)

event_set <- dm_dysreg_events_union
deltapsi_data <- Reduce(function(...) merge(..., by = "isoforms", all = TRUE),
                        list(select(allControls_res, gene_symbol, event_name, isoforms, delta_psi) %>% filter(isoforms %in% event_set),
                             select(tibialis_res, isoforms, delta_psi_mean) %>% filter(isoforms %in% event_set),
                             select(quadricep_res, isoforms, delta_psi_mean) %>% filter(isoforms %in% event_set)))
setnames(deltapsi_data, c("isoforms", "gene_symbol", "event_name", "Quad_vs_Tibialis_deltapsi", "Tibialis", "Quad"))

ggplot(deltapsi_data, aes(x = Tibialis, y = Quad, colour = cut(deltapsi_data$Quad_vs_Tibialis_deltapsi, seq(-1,1, 0.2)))) +
    geom_point() +
        geom_abline(yintercept = 0, slope = 1, linetype = "longdash") +
            labs(y = "Quad Delta PSI \n DM versus Control",
                 x = "Tibialis Delta PSI \n DM versus Control") +
                     scale_color_discrete(name = "Quad_vs_Tibialis_deltapsi")

fit <- lm(Tibialis ~ Quad, data = deltapsi_data)
cor.test(fit$residuals, deltapsi_data$Quad_vs_Tibialis_deltapsi)

ggplot() +
    geom_point(aes(y = deltapsi_data$Quad_vs_Tibialis_deltapsi, x = fit$residuals)) +
        #geom_abline(yintercept = 0, slope = 1, linetype = "longdash") +
            labs(x = "Residuals \n lm(DM_Quad_deltapsi ~ DM_Tibialis_deltapsi",
                 y = "Delta PSI \n Quad versus Tibials")


deltapsi_data_long <- gather(deltapsi_data, key = "Tissue", value = "DM_vs_Control_deltapsi", 5:6)

ggplot(deltapsi_data_long, aes(x = Quad_vs_Tibialis_deltapsi, y = DM_vs_Control_deltapsi, colour = Tissue)) +
    geom_point() +
        geom_hline(yintercept = 0, linetype = "longdash") +
            geom_vline(xintercept = 0, linetype = "longdash") +
                labs(y = "Delta PSI \n DM versus Control",
                     x = "Delta PSI \n Control Quad versus Control Tibialis")

ggplot(deltapsi_data_long, aes(y = DM_vs_Control_deltapsi, x = cut(deltapsi_data_long$Quad_vs_Tibialis_deltapsi, seq(-1,1, 0.1)))) +
    geom_boxplot(aes(fill = Tissue)) +
        labs(y = "Delta PSI \n DM versus Control",
             x = "Delta PSI \n Control Quad versus Control Tibialis")


ggplot(deltapsi_data_long, aes(x = DM_vs_Control_deltapsi, fill = cut(deltapsi_data_long$Quad_vs_Tibialis_deltapsi, seq(-1,1, 0.1)))) +
    geom_bar(position = "stack") +
        scale_fill_discrete(name = "Quad_vs_Tibialis_deltapsi")
    
##----------------------------------
## ontology analysis
##----------------------------------
all_muscle_genes <- filter(allControls_res, Quad_n / max(Quad_n, na.rm = TRUE) >= 0.75,
                           Tibialis_n / max(Tibialis_n, na.rm = TRUE) >= 0.75)$gene_symbol %>%
                               unique()

## For genes differentially regulated in healthy quad compared to healthy tibialis
sig_genes <- unique(quad_vs_tibialis$gene_symbol)
## different in all comparisons
sig_genes <- Reduce(intersect, list(dm_tibialis$gene_symbol, dm_quadricep$gene_symbol, quad_vs_tibialis$gene_symbol))
## different in DM quad OR tibialis AND between controls
sig_genes <- intersect(union(dm_tibialis$gene_symbol, dm_quadricep$gene_symbol), quad_vs_tibialis$gene_symbol)

ontology_class <- "BP"
ontology_class <- "MF"

myGO2genes <- AnnotationDbi::select(org.Hs.eg.db, keys = all_muscle_genes, columns=c("ENSEMBL", "GO"), keytype="SYMBOL")
myGO2genes <- myGO2genes %>% filter(!is.na(ENSEMBL)) %>% tbl_df
myGO2genesList <- tapply(filter(myGO2genes, ONTOLOGY == ontology_class)$ENSEMBL, filter(myGO2genes, ONTOLOGY == ontology_class)$GO, FUN = c)

ensemblIDs <- myGO2genes$ENSEMBL[match(all_muscle_genes, myGO2genes$SYMBOL)]
geneList <- factor(as.integer(all_muscle_genes %in% sig_genes))
names(geneList) <- ensemblIDs

GOdata <- new("topGOdata", description = "GO analysis of genes with differential splicing",
              ontology = ontology_class, allGenes = geneList, nodeSize = 5,
              annot = annFUN.GO2genes, GO2genes = myGO2genesList)

resultFisher <- runTest(GOdata, algorithm = "classic", statistic = "fisher")
resultFisher.elim <- runTest(GOdata, algorithm = "elim", statistic = "fisher")

allRes <- GenTable(GOdata, classicFisher = resultFisher,
                   elimFisher = resultFisher.elim,
                   orderBy = "elimFisher", ranksOf = "classicFisher",
                   topNodes = 20)

kable(allRes)

## plot GO graph
showSigOfNodes(GOdata, score(resultFisher), firstSigNodes = 5, useInfo = "all")

##----------------------------------
## compare to gene set that correlates
## with DM phenotypes
##----------------------------------
tibialis_cor <- tbl_dt(fread(paste("~/Projects/DMseq/results/tibialis/tibialis", event_type, "cor.txt", sep = "_")))

strength_cor_events <- tibialis_cor %>%
    filter(pval_HG.QMT < 0.05 & cor_HG.QMT > 0.5 | pval_ADF.QMT < 0.05 & cor_ADF.QMT > 0.5 | pval_Actual.Strength.6pt.Scale < 0.05 & cor_Actual.Strength.6pt.Scale > 0.5,
           isoforms %in% dm_tibialis$isoforms) %>% select(event_name, gene_symbol, cor_HG.QMT, cor_ADF.QMT, cor_Actual.Strength.6pt.Scale, isoforms)

strength_cor_events %>% filter(gene_symbol %in% sig_genes)

strength_cor_events %>% filter(isoforms %in% intersect(dm_tibialis$isoforms, quad_vs_tibialis$isoforms)) %>% arrange(desc(cor_ADF.QMT))
