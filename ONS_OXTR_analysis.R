# Methylation Analysis ONS OXTR

# dependencies
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("NanoMethViz")
library("NanoMethViz")

BiocManager::install("DSS")
library("DSS")

# data

methy_tabix <- file.path(tempdir(), "methy_data.bgz")

methy_calls <- c("/pathway/to/file/ID01_methylation_calls.tsv",
                 "/pathway/to/file/ID02_methylation_calls.tsv",
                 "/pathway/to/file/ID03_methylation_calls.tsv",
                 "/pathway/to/file/ID04_methylation_calls.tsv",
                 "/pathway/to/file/ID05_methylation_calls.tsv",
                 "/pathway/to/file/(...)")

create_tabix_file(methy_calls, methy_tabix)

hg38_exontable <- get_exons_hg38() 

#table w/ annotation of sample IDs, group 

sample <- c(
  "ID01_methylation_calls",
  "ID02_methylation_calls",
  "ID03_methylation_calls",
  "ID04_methylation_calls",
  "ID05_methylation_calls",
  "(...)")

group <- c("HFA",
           "Control",
           "HFA",
           "HFA",
           "HFA",
           "(...)")

sample_anno <- data.frame(sample, group, stringsAsFactors = FALSE) 

#NanoMethResults file
nmeth_results <- NanoMethResult(methy_tabix, sample_anno, hg38_exontable)

#plot methylation overview by group / spaghetti plot
plot_region(
  nmeth_results,
  "chr3",
  8750381,
  8770434,
  anno_regions = NULL,
  spaghetti = T)

# BSSeq type object for further analysis using NanoMethViz methy_to_bsseq interface

BSSeq_nmeth_results <- methy_to_bsseq(nmeth_results, out_folder = tempdir())

# coverage
# annotations
# OXTR 
region <- GRanges(seqnames = c("chr3"), 
                  ranges = IRanges(start = 8750381, end = 8770434))

#"nested"approach: chr3:  8.770.937 – 8.763.067 and chr3:  8.761.480 – 8.749.344
#fragment1
region_fragm1 <- GRanges(seqnames = c("chr3"), 
                         ranges = IRanges(start = 8749344, end = 8761480))
#fragment2
region_fragm2 <- GRanges(seqnames = c("chr3"), 
                         ranges = IRanges(start = 8763067, end = 8770937))

# OXTR 
coverage_OXTR <- getCoverage(BSSeq_nmeth_results, regions = region, what = "perRegionAverage")
mean(coverage_OXTR)
max(coverage_OXTR)
min(coverage_OXTR)
sd(coverage_OXTR)

#fragment1
coverage_OXTR_frgm1 <- getCoverage(BSSeq_nmeth_results, regions = region_fragm1, what = "perRegionAverage")
mean(coverage_OXTR_frgm1)
max(coverage_OXTR_frgm1)
min(coverage_OXTR_frgm1)
sd(coverage_OXTR_frgm1)

#fragment2
coverage_OXTR_frgm2 <- getCoverage(BSSeq_nmeth_results, regions = region_fragm2, what = "perRegionAverage")
mean(coverage_OXTR_frgm2)
max(coverage_OXTR_frgm2)
min(coverage_OXTR_frgm2)
sd(coverage_OXTR_frgm2)

# bygroup
# OXTR
dataframe_coverage_OXTR <- data.frame(coverage_OXTR)
dataframe_coverage_OXTR_tra <- transpose(dataframe_coverage_OXTR)
rownames(dataframe_coverage_OXTR_tra) <- sample
colnames(dataframe_coverage_OXTR_tra) <- "mean_coverage"
dataframe_coverage_OXTR_tra$Group <- group

aggregate(dataframe_coverage_OXTR_tra$mean_coverage, list(dataframe_coverage_OXTR_tra$Group), FUN=mean)
aggregate(dataframe_coverage_OXTR_tra$mean_coverage, list(dataframe_coverage_OXTR_tra$Group), FUN=sd)
aggregate(dataframe_coverage_OXTR_tra$mean_coverage, list(dataframe_coverage_OXTR_tra$Group), FUN=min)
aggregate(dataframe_coverage_OXTR_tra$mean_coverage, list(dataframe_coverage_OXTR_tra$Group), FUN=max)

t.test(mean_coverage ~ Group, data = dataframe_coverage_OXTR_tra, var.equal = FALSE, alternative = "two.sided")

# frgm1
dataframe_coverage_OXTR_frgm1 <- data.frame(coverage_OXTR_frgm1)
dataframe_coverage_OXTR_tra_frgm1 <- transpose(dataframe_coverage_OXTR_frgm1)
rownames(dataframe_coverage_OXTR_tra_frgm1) <- sample
colnames(dataframe_coverage_OXTR_tra_frgm1) <- "mean_coverage"
dataframe_coverage_OXTR_tra_frgm1$Group <- group

aggregate(dataframe_coverage_OXTR_tra_frgm1$mean_coverage, list(dataframe_coverage_OXTR_tra_frgm1$Group), FUN=mean)
aggregate(dataframe_coverage_OXTR_tra_frgm1$mean_coverage, list(dataframe_coverage_OXTR_tra_frgm1$Group), FUN=sd)
aggregate(dataframe_coverage_OXTR_tra_frgm1$mean_coverage, list(dataframe_coverage_OXTR_tra_frgm1$Group), FUN=min)
aggregate(dataframe_coverage_OXTR_tra_frgm1$mean_coverage, list(dataframe_coverage_OXTR_tra_frgm1$Group), FUN=max)

t.test(mean_coverage ~ Group, data = dataframe_coverage_OXTR_tra_frgm1, var.equal = FALSE, alternative = "two.sided")

# frgm2
dataframe_coverage_OXTR_frgm2 <- data.frame(coverage_OXTR_frgm2)
dataframe_coverage_OXTR_tra_frgm2 <- transpose(dataframe_coverage_OXTR_frgm2)
rownames(dataframe_coverage_OXTR_tra_frgm2) <- sample
colnames(dataframe_coverage_OXTR_tra_frgm2) <- "mean_coverage"
dataframe_coverage_OXTR_tra_frgm2$Group <- group

aggregate(dataframe_coverage_OXTR_tra_frgm2$mean_coverage, list(dataframe_coverage_OXTR_tra_frgm2$Group), FUN=mean)
aggregate(dataframe_coverage_OXTR_tra_frgm2$mean_coverage, list(dataframe_coverage_OXTR_tra_frgm2$Group), FUN=sd)
aggregate(dataframe_coverage_OXTR_tra_frgm2$mean_coverage, list(dataframe_coverage_OXTR_tra_frgm2$Group), FUN=min)
aggregate(dataframe_coverage_OXTR_tra_frgm2$mean_coverage, list(dataframe_coverage_OXTR_tra_frgm2$Group), FUN=max)

t.test(mean_coverage ~ Group, data = dataframe_coverage_OXTR_tra_frgm2, var.equal = FALSE, alternative = "two.sided")

#plots

dataframe_coverage_OXTR <- data.frame(coverage_OXTR)
dataframe_coverage_OXTR_tra <- transpose(dataframe_coverage_OXTR)
rownames(dataframe_coverage_OXTR_tra) <- sample
colnames(dataframe_coverage_OXTR_tra) <- "mean_coverage"
dataframe_coverage_OXTR_tra$Group <- group

# Figure 1A
OXTR_coverage_histogram <- ggplot(dataframe_coverage_OXTR_tra, aes(x=mean_coverage, fill=group)) +
  geom_histogram(binwidth = 2, colour="black", position = 'identity') +
  scale_fill_manual(values=c("grey30", "grey80")) +
  geom_vline(aes(xintercept=mean(mean_coverage, na.rm=T)),   # Ignore NA values for mean
             color="red", linetype="dashed", linewidth=1)+ 
  labs(y="count", 
       x="mean coverage", 
       title="Distribution - Coverage per sample",
       fill="Group")

# Figure 2B DNA conc

# annotate dna_conc <- 
c(71.07,
146.4,
187.3,
84.47,
(...))
dataframe_coverage_OXTR_tra$dna_conc <- dna_conc

OXTR_coverage_corrdnaconc <- ggplot(dataframe_coverage_OXTR_tra, aes(x=mean_coverage, y=dna_conc)) + 
  geom_point()+
  geom_smooth(method=lm)+
  stat_cor(method = "pearson", label.x = 3, label.y = 30)+ 
  labs(y="DNA concentration (ng/µl)", 
       x="mean coverage", 
       title="Correlation - Coverage and DNA concentration")

# Figure 1 C & D
coverage_OXTR_bybase_frgm1 <- getCoverage(BSSeq_OXTRsubset, regions = region_fragm1, what = "perBase")
dataframe_OXTR_coverage_bybase_frgm1 <- as.data.frame(coverage_OXTR_bybase_frgm1)

coverage_OXTR_bybase_frgm2 <- getCoverage(BSSeq_OXTRsubset, regions = region_fragm2, what = "perBase")
dataframe_OXTR_coverage_bybase_frgm2 <- as.data.frame(coverage_OXTR_bybase_frgm2)

OXTR_coverage_byCpG_frgm1 <- ggplot(dataframe_OXTR_coverage_bybase_frgm1, aes(x=CpG_position, y=mean)) + 
  geom_point()+ 
  geom_smooth(method="loess", se=F) + 
  xlim(c(8749344, 8761480)) + 
  ylim(c(15, 25)) + 
  labs(subtitle="Fragment chr3: 8.761.480 – 8.749.344", 
       y="Coverage", 
       x="Genomic position (chr3:)", 
       title="Coverage across CpG positions")

OXTR_coverage_byCpG_frgm2 <- ggplot(dataframe_OXTR_coverage_bybase_frgm2, aes(x=CpG_position, y=mean)) + 
  geom_point()+ 
  geom_smooth(method="loess", se=F) + 
  xlim(c(8763067, 8770937)) + 
  ylim(c(25, 40)) + 
  labs(subtitle="Fragment chr3:  8.770.937 – 8.763.067", 
       y="Coverage", 
       x="Genomic position (chr3:)", 
       title="Coverage across CpG positions")+
  geom_vline(aes(xintercept= 8769033),
             color="red", linetype="dashed", linewidth=0.5)+
  geom_vline(aes(xintercept= 8769438),
             color="red", linetype="dashed", linewidth=0.5)+
  geom_vline(aes(xintercept= 8768329),
             color="blue", linetype="dashed", linewidth=0.5)+
  geom_vline(aes(xintercept= 8767266),
             color="blue", linetype="dashed", linewidth=0.5)

# DM analysis
# subset of target sequence

BSSeq_OXTRsubset <- subsetByOverlaps(BSSeq_nmeth_results, GRanges(seqnames = "chr3",
                                                                  ranges = IRanges(start = 8750381, end = 8770434)))

#group annotations
HFAgroup <- c("ID01_methylation_calls",
              "ID03_methylation_calls",
              "ID04_methylation_calls",
              "ID05_methylation_calls",
              "(...)")

Controlgroup <- c("ID02_methylation_calls",
                  "ID08_methylation_calls",
                  "ID09_methylation_calls",
                  "ID11_methylation_calls",
                  "ID12_methylation_calls",
                  "(...)")
                  
malegroup <- c("ID01_methylation_calls",
               "ID03_methylation_calls",
               "ID04_methylation_calls",
               "ID06_methylation_calls",
               "(...)")
             
femalegroup <- c("ID02_methylation_calls",
                 "ID05_methylation_calls",
                 "ID07_methylation_calls",
                 "ID09_methylation_calls",
                 "(...)")

# DML test
# by group

dml_OXTR <- DMLtest(BSSeq_OXTRsubset, group1 = HFAgroup, 
                    group2 = Controlgroup,
                    smoothing=F, 
                    ncores=1)

# by sex

dml_OXTR_bysex <- DMLtest(BSSeq_OXTRsubset, group1 = malegroup, 
                    group2 = femalegroup,
                    smoothing=F, 
                    ncores=1)

# DML calling
#by group

dmls_nmeth_OXTR_results <- callDML(dml_OXTR, p.threshold = 0.05, delta = 0.1)

# by sex

dmls_nmeth_OXTR_bysex_results <- callDML(dml_OXTR_bysex, p.threshold = 0.05, delta = 0.1)

# DMR calling
# by group

dmrs_nmeth_OXTR_results = callDMR(dml_OXTR, p.threshold = 0.05, delta=0.1)

# by sex

dmrs_nmeth_OXTR_results = callDMR(dml_OXTR_bysex, p.threshold = 0.05, delta=0.1)

# Figure 3
OXTR_MT2_plotregions <- c(GRanges(seqnames = "chr3",
                                  ranges = IRanges(start = 8769033, end = 8769438)),
                          GRanges(seqnames = "chr3",
                                  ranges = IRanges(start = 8767266, end = 8768329)),
                          GRanges(seqnames = "chr3",
                                  ranges = IRanges(start = 8769088, end = 8769088)),
                          GRanges(seqnames = "chr3",
                                  ranges = IRanges(start = 8769111, end = 8769111)),
                          GRanges(seqnames = "chr3",
                                  ranges = IRanges(start = 8769121, end = 8769121)))

# dataframe with mean methylation per sample by region
rawmethCI_OXTR_regions <- getMeth(BSSeq_OXTRsubset, regions = OXTR_MT2_plotregions, type = "raw", what = "perRegion", confint = F)
rawmethCI_OXTR_df_regions <- as.data.frame(rawmethCI_OXTR_regions)

# invert table + group annotation
rawmethCI_OXTR_df_regions_tra <- transpose(rawmethCI_OXTR_df_regions)
rawmethCI_OXTR_df_regions_tra$Group <- group
colnames(rawmethCI_OXTR_df_regions_tra) <- c("MT2", "Exon3", "CpG901", "CpG924", "CpG934", "Group")
rownames(rawmethCI_OXTR_df_regions_tra) <- sample

# boxplot/violinplot
# MT2
MT2 <- ggplot(rawmethCI_OXTR_df_regions_tra, aes(x= Group, y= MT2, fill= Group)) +
  geom_violin(width=1.0) +
  geom_boxplot(width=0.3, color="black", alpha=0.5) +
  labs(title = "MT2 (chr3: 8769033-8769438)",
       x = "Group",
       y = "Methylation Rate") +
  theme(plot.title = element_text(hjust = 0.5, size = 12, face = "bold"))
# Exon3
Exon3 <- ggplot(rawmethCI_OXTR_df_regions_tra, aes(x= Group, y= Exon3, fill= Group)) +
  geom_violin(width=1.0) +
  geom_boxplot(width=0.3, color="black", alpha=0.5) +
  labs(title = "Exon 3 (chr3: 8767266-8768329)",
       x = "Group",
       y = "Methylation Rate") +
  theme(plot.title = element_text(hjust = 0.5, size = 12, face = "bold"))
# CpG901
CpG901 <- ggplot(rawmethCI_OXTR_df_regions_tra, aes(x= Group, y= CpG901, fill= Group)) +
  geom_violin(width=1.0) +
  geom_boxplot(width=0.3, color="black", alpha=0.5) +
  labs(title = "CpG -901 (chr3: 8769088)",
       x = "Group",
       y = "Methylation Rate") +
  theme(plot.title = element_text(hjust = 0.5, size = 12, face = "bold"))
# CpG924
CpG924 <- ggplot(rawmethCI_OXTR_df_regions_tra, aes(x= Group, y= CpG924, fill= Group)) +
  geom_violin(width=1.0) +
  geom_boxplot(width=0.3, color="black", alpha=0.5) +
  labs(title = "CpG -924 (chr3: 8769111)",
       x = "Group",
       y = "Methylation Rate") +
  theme(plot.title = element_text(hjust = 0.5, size = 12, face = "bold"))
# CpG934
CpG934 <- ggplot(rawmethCI_OXTR_df_regions_tra, aes(x= Group, y= CpG934, fill= Group)) +
  geom_violin(width=1.0) +
  geom_boxplot(width=0.3, color="black", alpha=0.5) +
  labs(title = "CpG -934 (chr3: 8769121)",
       x = "Group",
       y = "Methylation Rate") +
  theme(plot.title = element_text(hjust = 0.5, size = 12, face = "bold"))

# correlation plots
BSSeq_OXTR_CpGIsland_subset <- subsetByOverlaps(BSSeq_nmeth_results,
                                                GRanges(seqnames = "chr3",
                                                        ranges = IRanges(start = 8767276, end = 8769594)))

BSSeq_OXTR_MT2_subset <- subsetByOverlaps(BSSeq_nmeth_results,
                                          GRanges(seqnames = "chr3",
                                                  ranges = IRanges(start = 8769033, end = 8769438)))

lmrmatrix_OXTR_CpGIsland <- bsseq_to_log_methy_ratio(BSSeq_OXTR_CpGIsland_subset, prior_count = 2)
lmrmatrix_OXTR_CpGIsland_df <- as.data.frame(lmrmatrix_OXTR_CpGIsland)
lmrmatrix_OXTR_CpGIsland_df_tra <- transpose(lmrmatrix_OXTR_CpGIsland_df)
rownames(lmrmatrix_OXTR_CpGIsland_df_tra) <- sample
lmrmatrix_OXTR_CpGIsland_df_tra$Group <- group
lmrmatrix_OXTR_CpGIsland_df_tra_HFA <- subset(lmrmatrix_OXTR_CpGIsland_df_tra, Group == "HFA")
lmrmatrix_OXTR_CpGIsland_df_tra_Control <- subset(lmrmatrix_OXTR_CpGIsland_df_tra, Group == "Control")
lmrmatrix_OXTR_CpGIsland_df_tra_HFA$Group <- NULL
lmrmatrix_OXTR_CpGIsland_df_tra_Control$Group <- NULL
colnames(lmrmatrix_OXTR_CpGIsland_df_tra_HFA) <- OXTR_Island_CpGs
colnames(lmrmatrix_OXTR_CpGIsland_df_tra_Control) <- OXTR_Island_CpGs
corrOXTR_HFA <- cor(lmrmatrix_OXTR_CpGIsland_df_tra_HFA)
corrOXTR_Control <- cor(lmrmatrix_OXTR_CpGIsland_df_tra_Control)

#MT2 subset
lmrmatrix_OXTR_MT2 <- bsseq_to_log_methy_ratio(BSSeq_OXTR_MT2_subset, prior_count = 2)
lmrmatrix_OXTR_MT2_df <- as.data.frame(lmrmatrix_OXTR_MT2)
lmrmatrix_OXTR_MT2_df_tra <- transpose(lmrmatrix_OXTR_MT2_df)
rownames(lmrmatrix_OXTR_MT2_df_tra) <- sample
lmrmatrix_OXTR_MT2_df_tra$Group <- group
lmrmatrix_OXTR_MT2_df_tra_HFA <- subset(lmrmatrix_OXTR_MT2_df_tra, Group == "HFA")
lmrmatrix_OXTR_MT2_df_tra_Control <- subset(lmrmatrix_OXTR_MT2_df_tra, Group == "Control")
lmrmatrix_OXTR_MT2_df_tra_HFA$Group <- NULL
lmrmatrix_OXTR_MT2_df_tra_Control$Group <- NULL
OXTR_MT2_CpGs <- rownames(lmrmatrix_OXTR_MT2_df)
colnames(lmrmatrix_OXTR_MT2_df_tra_HFA) <- OXTR_Island_CpGs
colnames(lmrmatrix_OXTR_MT2_df_tra_Control) <- OXTR_Island_CpGs
corrOXTR_HFA <- cor(lmrmatrix_OXTR_MT2_df_tra_HFA)
corrOXTR_Control <- cor(lmrmatrix_OXTR_MT2_df_tra_Control)
colnames(corrOXTR_HFA) <- OXTR_MT2_CpGs
colnames(corrOXTR_Control) <- OXTR_MT2_CpGs
rownames(corrOXTR_HFA) <- OXTR_MT2_CpGs
rownames(corrOXTR_Control) <- OXTR_MT2_CpGs

par(mfrow=c(1,2))
corrplot(corrOXTR_HFA, method = "color", type = "upper", tl.pos = "l", tl.col = "black", tl.cex = 0.5, title= "HFA", mar=c(0,0,1,0))
corrplot(corrOXTR_Control, method = "color", type = "upper", tl.pos = "l", tl.col = "black", tl.cex = 0.5, title = "Controls", mar=c(0,0,1,0))
