library(tidypopgen)
library(admixtools)
library(ggplot2)

setwd("/home/patel/embo_popgen_2025/!_MyWork/Day5")
dir("./data")


modern_gt <- tidypopgen::gen_tibble(
  "./data/modern_samples.bed",valid_alleles = c("A","T","C","G"),
  missing_alleles = c("X"))

modern_gt %>% print(n=21)
unique(modern_gt$population)

modern_gt %>% group_by(population) %>% tally() %>% print(n=21)

loci_report <- modern_gt %>%
  group_by(population) %>%
  qc_report_loci()

png("./plots/loci_report.png")
autoplot(loci_report, type = "missing")
dev.off()

modern_gt <- modern_gt %>% 
  select_loci_if(loci_missingness(genotypes)<0.04)

ancient_gt <- tidypopgen::gen_tibble(
  "./data/ancient_samples.vcf",
  valid_alleles = c("A","T","C","G","X"), 
  quiet = TRUE)

ancient_gt <- gt_pseudohaploid(ancient_gt)

ancient_gt$id

ancient_gt$id[ancient_gt$id == "GB20"] <- "Mota"

ancient_gt$population <- ancient_gt$id

merged_dry <- rbind_dry_run(modern_gt, ancient_gt, flip_strand = TRUE)

merged_gt <- rbind(modern_gt, ancient_gt, 
                     flip_strand = TRUE,
                     backingfile = "./data/merged_samples")

merged_gt <- merged_gt %>% group_by(population)

f2_dir <- "./data/f2_tidypopgen"

# compute f2
f2_tidypopgen <- gt_extract_f2(merged_gt, 
                               outdir = "./data/f2_tidypopgen", 
                               overwrite = TRUE)

f2_blocks = f2_from_precomp("./data/f2_tidypopgen")

lbk_modern_panel <- c("Basque", "Bedouin2", "Druze", "Cypriot", "Tuscan",
  "Sardinian", "French", "Spanish", "Onge", "Han", "Mayan", "Mixe", "Surui")

lbk_f3out <- f3(data = f2_blocks, 
# outgroup
                pop1 = "Mbuti",  
                # group of interest
                pop2 = "LBK",
                # comparison set 
                pop3 = lbk_modern_panel)
lbk_f3out

lbk_f3out %>% arrange(desc(est))


ggplot(lbk_f3out, aes(pop3, est)) +
  geom_point() +
  geom_errorbar(aes(ymin = est - 2 * se, ymax = est + 2 * se)) +
  labs(y = "Shared drift with LBK", x = "populations") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))
ggsave("./plots/f3_LBK.png")

lbk_f3out$pop3<-factor(lbk_f3out$pop3, levels = lbk_f3out$pop3[order(lbk_f3out$est)])

ggplot(lbk_f3out, aes(pop3, est)) +
  geom_point() +
  geom_errorbar(aes(ymin = est - 2 * se, ymax = est + 2 * se)) +
  labs(y = "Shared drift with LBK", x = "populations") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))
ggsave("./plots/f3_LBK.png")

# MA1 comparison
lbk_f3out <- f3(data = f2_blocks, 
# outgroup
                pop1 = "Mbuti",  
                # group of interest
                pop2 = "MA1",
                # comparison set 
                pop3 = lbk_modern_panel)
lbk_f3out
lbk_f3out$pop3<-factor(lbk_f3out$pop3, levels = lbk_f3out$pop3[order(lbk_f3out$est)])

ggplot(lbk_f3out, aes(pop3, est)) +
  geom_point() +
  geom_errorbar(aes(ymin = est - 2 * se, ymax = est + 2 * se)) +
  labs(y = "Shared drift with LBK", x = "populations") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))
ggsave("./plots/f3_MA1.png")

aa_f3admix <- f3(data = f2_blocks,
                # admixed pop
                 pop1 = "AA",
                #  contributing pop 1
                 pop2 = "Yoruba",
                #  contributing pop 2
                 pop3 = "French")
aa_f3admix

eurasian_sources <- c("French","Spanish","Sardinian","LBK")

somali_f3admix <- f3(data = f2_blocks,
                     pop1 = "Somali",
                     pop2 = eurasian_sources, 
                     pop3 = "Mota")
somali_f3admix

dinka_f3admix <- f3(data = f2_blocks,
                     pop1 = "Dinka",
                     pop2 = eurasian_sources, 
                     pop3 = "Mota")
dinka_f3admix

dinka_f3admix <- f3(data = f2_blocks,
                     pop1 = "Somali",
                     pop2 = "Dinka", 
                     pop3 = "Mota")
dinka_f3admix

neand_french_f4 <- f4(data = f2_blocks,
                     pop1 = "pan_troglodytes",
                     pop2 = "AltaiNea",
                     pop3 = "Mbuti",
                     pop4 = "French")
neand_french_f4

neand_han_f4 <- f4(data = f2_blocks,
                     pop1 = "pan_troglodytes",
                     pop2 = "AltaiNea",
                     pop3 = "Mbuti",
                     pop4 = "Han")
neand_han_f4

# yamnaya pop relationship to eur pops 
eurasian_sources <- c("French", "Basque", "Spanish")
somali_f3admix <- f4(data = f2_blocks,
                     pop1 = "pan_troglodytes",
                     pop2 = "Yamnaya",
                     pop3 = "LBK",
                     pop4 = eurasian_sources)
somali_f3admix

pops <- c("Han", "pan_troglodytes", "AA","Yoruba","French")
qpf4ratio(data = f2_blocks, 
          pops = pops)

pops <- c("Han", "pan_troglodytes", "Yoruba","AA","French")
qpf4ratio(data = f2_blocks, 
          pops = pops)

lbk_modern_panel <- c("Basque", "Bedouin2", "Druze", "Cypriot", "Tuscan",
  "Sardinian", "French", "Spanish", "Onge", "Han", "Mayan", "Mixe", "Surui")
modern_panel_gt <- merged_gt %>% filter (population %in% lbk_modern_panel)
# remove monomorphic sites (ungrouping the tibble, as we want to use global frequencies
# rather than population frequencies)
modern_panel_gt <- modern_panel_gt %>% ungroup() %>% select_loci_if(loci_maf(genotypes)>0)

modern_panel_gt <- gt_update_backingfile(modern_panel_gt)
# reset the ploidy for this tibble, as it is now all diploids
attr(modern_panel_gt$genotypes, "ploidy") <- 2
modern_panel_gt <- modern_panel_gt %>% gt_impute_simple(method="mode")

modern_panel_gt <- modern_panel_gt %>% select_loci_if(loci_ld_clump(genotypes))

modern_pca <- modern_panel_gt %>% gt_pca_randomSVD()

autoplot(modern_pca, type = "scores") 

library(ggplot2)
autoplot(modern_pca, type = "scores") +
  aes(color = modern_panel_gt$population) +
  labs(color = "population")


lbk_gt <- merged_gt %>% filter(id == "LBK")
lbk_pca_scores <- predict(modern_pca, new_data = lbk_gt, project_method = "least_square")
lbk_pca_scores

autoplot(modern_pca, type = "scores") +
  aes(color = modern_panel_gt$population) +
  labs(color = "population")+
  geom_point(data=lbk_pca_scores, mapping=aes(x=.data$.PC1, y=.data$.PC2), col = "black")


autoplot(modern_pca, type = "scores") +
  aes(color = modern_panel_gt$population) +
  labs(color = "population") +
  geom_point(data=lbk_pca_scores, mapping=aes(x=.data$.PC1, y=.data$.PC2), col = "black") +
  lims(x=c(30, 70), y = c(-10, 15))

