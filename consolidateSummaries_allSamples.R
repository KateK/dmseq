library(data.table)
library(dplyr)

##----------------------------------
## load event data
##----------------------------------
groups <- c("tibialis", "quadricep", "heart", "DMseq")
event_types <- c("nonUTRevents.multi", "ALE", "AFE", "polyA")

for (j in 1:length(event_types)) {
    event_type <- event_types[j]
    output_file <- paste("~/MISO/single_end_mode/summaries/allSamples_", event_type, "_consolidatedSummaries.txt", sep = "")
    allSamples <- data.table()
    for (i in 1: length(groups)) {
        group <- groups[i]
        con_data <- fread(paste("~/MISO/single_end_mode/summaries/", group, "/", group, "_", event_type, "_consolidatedSummaries.txt", sep = ""))
        allSamples <- rbind(allSamples, con_data)
    }
    write.table(allSamples, output_file, row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")
}
