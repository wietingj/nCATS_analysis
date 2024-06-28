methy_tabix <- file.path(tempdir(), "methy_data.bgz")

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("NanoMethViz")
library(NanoMethViz)

hg38_exontable <- get_exons_hg38() 

methy_calls <- c("/path/to/ID01_methylation_calls.tsv",
                 "/path/to/ID02_methylation_calls.tsv",
                 "/path/to/ID03_methylation_calls.tsv",
                 "/path/to/ID04_methylation_calls.tsv",
                 "/path/to/ID05_methylation_calls.tsv",
                 "/path/to/ID06_methylation_calls.tsv",
                 "/path/to/ID07_methylation_calls.tsv",
                 "/path/to/ID08_methylation_calls.tsv")
# ... to be continued ID09-ID40 ...

create_tabix_file(methy_calls, methy_tabix)

#group annotations 

group <- c("HFA-male",
           "NC-female",
           "HFA-male",
           "HFA-male",
           "HFA-female",
           "HFA-male",
           "HFA-female",
           "NC-male",
           "NC-female",
           "HFA-female",
           "NC-female",
           "NC-male",
           "NC-female",
           "NC-male",
           "HFA-male",
           "NC-male",
           "HFA-male",
           "NC-female",
           "HFA-male",
           "HFA-female",
           "NC-male",
           "HFA-female",
           "HFA-female",
           "NC-female",
           "HFA-female",
           "NC-male",
           "NC-female",
           "NC-male",
           "NC-male",
           "HFA-male",
           "HFA-female",
           "HFA-male",
           "NC-male",
           "NC-male",
           "HFA-female",
           "HFA-female",
           "HFA-male",
           "NC-female",
           "NC-female",
           "NC-female")

sample_anno <- data.frame(sample, group, stringsAsFactors = FALSE) 

nmeth_results <- NanoMethResult(methy_tabix, sample_anno, hg38_exontable)

BSSeq_nmeth_results <- methy_to_bsseq(nmeth_results, out_folder = tempdir())

BSSeq_nmeth_results <- methy_to_bsseq(nmeth_results)

#Basic plot MAGEL2
# guides 23639316 23651466
plot_region(
  nmeth_results,
  "chr15",
  23639316,
  23651466,
  anno_regions = NULL,
  binary_threshold = NULL,
  avg_method = "mean",
  spaghetti = T,
  heatmap = F,
  span = NULL,
  window_prop = 0,
  palette = ggplot2::scale_colour_brewer(palette = "Set1"),
  line_size = 1
)

#subset via BSseq function subsetByOverlaps (can either be BSseq or GRanges)
BSSeq_MAGEL2subset <- subsetByOverlaps(BSSeq_nmeth_results, GRanges(seqnames = "chr15",
                                                                  ranges = IRanges(start = 23639316, end = 23651466)))

#DML-Testing MAGEL2
BiocManager::install("DSS")
library("DSS")

# this is hfa only
malegroup <- c("ID01_methylation_calls",
               "ID03_methylation_calls",
               "ID04_methylation_calls",
               "ID06_methylation_calls",
               "ID15_methylation_calls",
               "ID17_methylation_calls",
               "ID19_methylation_calls",
               "ID30_methylation_calls",
               "ID32_methylation_calls",
               "ID37_methylation_calls")

femalegroup <- c("ID05_methylation_calls",
                 "ID07_methylation_calls",
                 "ID10_methylation_calls",
                 "ID20_methylation_calls",
                 "ID22_methylation_calls",
                 "ID23_methylation_calls",
                 "ID25_methylation_calls",
                 "ID31_methylation_calls",
                 "ID35_methylation_calls",
                 "ID36_methylation_calls")

#DML Testing
dml_MAGEL2_hfaonly_sex <- DMLtest(BSSeq_MAGEL2subset, group1 = malegroup, 
                    group2 = femalegroup,
                    smoothing=T, 
                    ncores=1)
dml_MAGEL2_sex_hfaonly <- callDML(dml_MAGEL2_hfaonly_sex, delta = 0.05, p.threshold = 0.05)
dmr_MAGEL2_sex_hfaonly <- callDMR(dml_MAGEL2_hfaonly_sex)

#coverage
#whole MAGEL2 
region <- GRanges(seqnames = c("chr15"), 
          ranges = IRanges(start = 23639316, end = 23651466))

#overall
coverage_MAGEL2 <- getCoverage(BSSeq_nmeth_results, regions = region, what = "perRegionAverage")
mean(coverage_MAGEL2)
median(coverage_MAGEL2)
max(coverage_MAGEL2)
min(coverage_MAGEL2)
sd(coverage_MAGEL2)

#bygroup
dataframe_coverage_MAGEL2 <- data.frame(coverage_MAGEL2)
dataframe_coverage_MAGEL2_tra <- transpose(dataframe_coverage_MAGEL2)
rownames(dataframe_coverage_MAGEL2_tra) <- sample
colnames(dataframe_coverage_MAGEL2_tra) <- "mean_coverage"
dataframe_coverage_MAGEL2_tra$Group <- group
aggregate(dataframe_coverage_MAGEL2_tra$mean_coverage, list(dataframe_coverage_MAGEL2_tra$Group), FUN=mean)
aggregate(dataframe_coverage_MAGEL2_tra$mean_coverage, list(dataframe_coverage_MAGEL2_tra$Group), FUN=median)
aggregate(dataframe_coverage_MAGEL2_tra$mean_coverage, list(dataframe_coverage_MAGEL2_tra$Group), FUN=sd)
aggregate(dataframe_coverage_MAGEL2_tra$mean_coverage, list(dataframe_coverage_MAGEL2_tra$Group), FUN=min)
aggregate(dataframe_coverage_MAGEL2_tra$mean_coverage, list(dataframe_coverage_MAGEL2_tra$Group), FUN=max)
t.test(mean_coverage ~ Group, data = dataframe_coverage_MAGEL2_tra_frgm1, var.equal = FALSE, alternative = "two.sided")

#"nested"approach
#fragment1
region_fragm1 <- GRanges(seqnames = c("chr15"), 
                         ranges = IRanges(start = 23639316, end = 23647107))
#fragment2
region_fragm2 <- GRanges(seqnames = c("chr15"), 
                         ranges = IRanges(start = 23647429, end = 23651466))

#fragment1
coverage_MAGEL2_frgm1 <- getCoverage(BSSeq_nmeth_results, regions = region_fragm1, what = "perRegionAverage")
mean(coverage_MAGEL2_frgm1)
median(coverage_MAGEL2_frgm1)
max(coverage_MAGEL2_frgm1)
min(coverage_MAGEL2_frgm1)
sd(coverage_MAGEL2_frgm1)

#fragment2
coverage_MAGEL2_frgm2 <- getCoverage(BSSeq_nmeth_results, regions = region_fragm2, what = "perRegionAverage")
mean(coverage_MAGEL2_frgm2)
median(coverage_MAGEL2_frgm2)
max(coverage_MAGEL2_frgm2)
min(coverage_MAGEL2_frgm2)
sd(coverage_MAGEL2_frgm2)

#bygroup frgm1
dataframe_coverage_MAGEL2_frgm1 <- data.frame(coverage_MAGEL2_frgm1)
dataframe_coverage_MAGEL2_tra_frgm1 <- transpose(dataframe_coverage_MAGEL2_frgm1)
rownames(dataframe_coverage_MAGEL2_tra_frgm1) <- sample
colnames(dataframe_coverage_MAGEL2_tra_frgm1) <- "mean_coverage"
dataframe_coverage_MAGEL2_tra_frgm1$Group <- group
aggregate(dataframe_coverage_MAGEL2_tra_frgm1$mean_coverage, list(dataframe_coverage_MAGEL2_tra_frgm1$Group), FUN=mean)
aggregate(dataframe_coverage_MAGEL2_tra_frgm1$mean_coverage, list(dataframe_coverage_MAGEL2_tra_frgm1$Group), FUN=median)
aggregate(dataframe_coverage_MAGEL2_tra_frgm1$mean_coverage, list(dataframe_coverage_MAGEL2_tra_frgm1$Group), FUN=sd)
aggregate(dataframe_coverage_MAGEL2_tra_frgm1$mean_coverage, list(dataframe_coverage_MAGEL2_tra_frgm1$Group), FUN=min)
aggregate(dataframe_coverage_MAGEL2_tra_frgm1$mean_coverage, list(dataframe_coverage_MAGEL2_tra_frgm1$Group), FUN=max)
t.test(mean_coverage ~ Group, data = dataframe_coverage_MAGEL2_tra_frgm1, var.equal = FALSE, alternative = "two.sided")

#bygroup frgm2
dataframe_coverage_MAGEL2_frgm2 <- data.frame(coverage_MAGEL2_frgm2)
dataframe_coverage_MAGEL2_tra_frgm2 <- transpose(dataframe_coverage_MAGEL2_frgm2)
rownames(dataframe_coverage_MAGEL2_tra_frgm2) <- sample
colnames(dataframe_coverage_MAGEL2_tra_frgm2) <- "mean_coverage"
dataframe_coverage_MAGEL2_tra_frgm2$Group <- group

aggregate(dataframe_coverage_MAGEL2_tra_frgm2$mean_coverage, list(dataframe_coverage_MAGEL2_tra_frgm2$Group), FUN=mean)
aggregate(dataframe_coverage_MAGEL2_tra_frgm2$mean_coverage, list(dataframe_coverage_MAGEL2_tra_frgm2$Group), FUN=median)
aggregate(dataframe_coverage_MAGEL2_tra_frgm2$mean_coverage, list(dataframe_coverage_MAGEL2_tra_frgm2$Group), FUN=sd)
aggregate(dataframe_coverage_MAGEL2_tra_frgm2$mean_coverage, list(dataframe_coverage_MAGEL2_tra_frgm2$Group), FUN=min)
aggregate(dataframe_coverage_MAGEL2_tra_frgm2$mean_coverage, list(dataframe_coverage_MAGEL2_tra_frgm2$Group), FUN=max)

t.test(mean_coverage ~ Group, data = dataframe_coverage_MAGEL2_tra_frgm2, var.equal = FALSE, alternative = "two.sided")

#plots
dataframe_coverage_MAGEL2 <- data.frame(coverage_MAGEL2)
dataframe_coverage_MAGEL2_tra <- transpose(dataframe_coverage_MAGEL2)
rownames(dataframe_coverage_MAGEL2_tra) <- sample
colnames(dataframe_coverage_MAGEL2_tra) <- "mean_coverage"

#group annotation
dataframe_coverage_MAGEL2_tra$Group <- group
dataframe_MAGEL2_coverage_bybase_frgm1_tra$Group <- group

#histogram by group
MAGEL2_coverage_histogram <- ggplot(dataframe_coverage_MAGEL2_tra, aes(x=mean_coverage, fill=group)) +
                           geom_histogram(binwidth = 2, colour="black", position = 'identity') +
                           scale_fill_manual(values=c("grey30", "grey80")) +
                           geom_vline(aes(xintercept=mean(mean_coverage, na.rm=T)),   # Ignore NA values for mean
                           color="black", linetype="dashed", linewidth=1)+ 
                           labs(y="count", 
                           x="Coverage", 
                           title="Coverage per sample",
                           fill="Group")

#perBaseCoverage
coverage_MAGEL2_bybase <- getCoverage(BSSeq_MAGEL2subset, regions = region, what = "perBase")
dataframe_MAGEL2_coverage_bybase <- as.data.frame(coverage_MAGEL2_bybase)

#CpGannotation from table
CpGs_MAGEL2 <- dml_MAGEL2$pos 

dataframe_MAGEL2_coverage_bybase$mean <- apply(dataframe_MAGEL2_coverage_bybase[,1:40], 1, mean)
dataframe_MAGEL2_coverage_bybase$CpG_position <- CpGs_MAGEL2

#plot coverage across CpGpos 

MAGEL2_coverage_byCpG <- ggplot(dataframe_MAGEL2_coverage_bybase, aes(x=CpG_position, y=mean)) + 
                         geom_point()+ 
                         xlim(c(23639316, 23651466)) + 
                         ylim(c(0, 60)) + 
                         geom_vline(aes(xintercept= 23647429),
                                     color="red", linetype="dashed", linewidth=0.5)+
                         geom_vline(aes(xintercept= 23647448),
                                     color="red", linetype="dashed", linewidth=0.5)+
                         geom_vline(aes(xintercept= 23651466),
                                   color="red", linetype="dashed", linewidth=0.5)+
                         geom_vline(aes(xintercept= 23651447),
                                   color="red", linetype="dashed", linewidth=0.5)+
                         geom_vline(aes(xintercept= 23647107),
                                     color="blue", linetype="dashed", linewidth=0.5)+
                         geom_vline(aes(xintercept= 23647088),
                                     color="blue", linetype="dashed", linewidth=0.5) + 
                         geom_vline(aes(xintercept= 23639316),
                                   color="blue", linetype="dashed", linewidth=0.5)+
                         geom_vline(aes(xintercept= 23639336),
                                   color="blue", linetype="dashed", linewidth=0.5) +
                         ggtitle("Coverage across CpG positions") +
                         ylab("Coverage") +
                         xlab("Genomic position (chr15:)")

#combined figure for publication
ggarrange(MAGEL2_coverage_histogram + labs(tag = "A"),
          MAGEL2_coverage_byCpG + labs(tag = "B"),
          nrow = 2, ncol = 1, common.legend = F)

#dataframe with mean methy per sample by region
rawmethCI_MAGEL2_regions <- getMeth(BSSeq_MAGEL2subset, regions = region, type = "raw", what = "perBase", confint = F)
rawmethCI_MAGEL2_df_regions <- as.data.frame(rawmethCI_MAGEL2_regions)
#
CpGs_MAGEL2 <- dml_MAGEL2$pos 
#
rawmethCI_MAGEL2_df_regions_tra <- transpose(rawmethCI_MAGEL2_df_regions)
colnames(rawmethCI_MAGEL2_df_regions_tra) <- CpGs_MAGEL2
rownames(rawmethCI_MAGEL2_df_regions_tra) <- sample_anno$sample

write_xlsx(rawmethCI_MAGEL2_df_regions_tra, "/path/to/HFA_MAGEL2_methy_overall_byCpG_perID.xlsx")

# shorten xlxs with mean columns of DMR region

rawmeth_MAGEL2_DMR_mean <- read_xlsx("path/to/HFA_MAGEL2_DMRonly_forboxplot.xlsx", col_names = T)
rawmeth_MAGEL2_DMR_mean$Group <- group

rawmeth_MAGEL2_DMR_mean_long <- melt(rawmeth_MAGEL2_DMR_mean)

boxplot <- ggplot(rawmeth_MAGEL2_DMR_mean_long, aes(x= Group, 
                                     y= value, 
                                     fill= Group)) +
  geom_violin(width=1.0) +
  geom_boxplot(width=0.3, color="black", alpha=0.5) +
  labs(title = "Differentially Methylated Region (Chr15:23647640 - 23647939)",
       x = "Group",
       y = "Methylation Rate") +
  theme(plot.title = element_text(hjust = 0.5, size = 12, face = "bold")) 

boxplot + scale_fill_manual(values=c("#E41A1C","#377EB8","#4DAF4A","#984EA3"))