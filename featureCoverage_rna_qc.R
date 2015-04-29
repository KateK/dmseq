##----------------------------------
## required libraries
##----------------------------------
library(ggplot2)
library(dplyr)
library(reshape2)


##----------------------------------
## 
##----------------------------------
feature_cov <- tbl_df(read.table("~/Projects/DMseq/data/allSamples_consolidatedAlignmentStats.txt", header = T, sep = "\t"))

plotsubset <- feature_cov %>% select(SAMPLE, PCT_RIBOSOMAL_BASES, PCT_CODING_BASES, PCT_UTR_BASES, PCT_INTRONIC_BASES, PCT_INTERGENIC_BASES)

feature_pcts <- melt(plotsubset, id.vars = "SAMPLE", measure.vars = c("PCT_RIBOSOMAL_BASES", "PCT_CODING_BASES", "PCT_UTR_BASES", "PCT_INTRONIC_BASES", "PCT_INTERGENIC_BASES"), variable.name = "REGION", value.name = "PERCENT")

pdata_file <- "~/Projects/DMseq/data/DM_sample_pdata.txt"
pdata <- read.table(pdata_file, header = TRUE, row.names = NULL, sep = "\t") 
names(pdata)[1] <- "SAMPLE"
feature_pcts_pdata <- left_join(feature_pcts, pdata, by = "SAMPLE")

## stacked barplot
postscript(file = "~/Projects/DMseq/data/allSamples_alignmenStats.eps")

ggplot(feature_pcts, aes(x = SAMPLE, y = PERCENT*100, fill = REGION)) +
    geom_bar(stat = "identity", ) +
        labs(y = "PERCENT", x = "") +
            scale_y_continuous(expand = c(0,0)) +
                scale_x_discrete(expand = c(0,0)) +
                    theme_bw(15) +
                        theme(axis.text.x=element_text(angle=90, hjust=1, vjust=0.5, size=6),
                              plot.margin = unit(c(1,0,1,0), "cm"),
                              legend.key.size = unit(0.25, "cm"),
                              legend.key.width = unit(0.25, "cm"),
                              legend.text = element_text(size = 7, hjust = 1, vjust = 1))

dev.off()


## feature_pcts %>% group_by(SAMPLE) %>% summarize(sum(PERCENT))

ggplot(feature_pcts_pdata, aes(x = Diagnosis, y = PERCENT*100, fill = Diagnosis)) +
    geom_boxplot() +
        labs(y = "PERCENT", x = "") +
            guides(fill = FALSE) +
                facet_wrap(~REGION, ncol = 5) +
                    theme_grey(15)
