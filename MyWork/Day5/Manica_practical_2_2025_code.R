library(admixtools)
library(tidypopgen)

setwd("/home/patel/embo_popgen_2025/!_MyWork/Day5")
dir("./data")

f2_blocks = f2_from_precomp("./data/f2_tidypopgen", verbose = FALSE)

neand_euras_wave <- qpwave(data = f2_blocks,
      left = c("French","Spanish","Tuscan"),
      right = c("AltaiNea","Mota", "Yoruba", "Denisova", "Mbuti")
)
neand_euras_wave

neand_euras_wave <- qpwave(data = f2_blocks,
      left = c("French","Spanish","Tuscan", "Han", "Onge"),
      right = c("AltaiNea","Mota", "Yoruba", "Denisova", "Mbuti")
)
neand_euras_wave

french_adm <- qpadm(data = f2_blocks,
      left = c("Loschbour", "LBK", "Yamnaya"),
      right = c("Mbuti", "Mota", "Dinka", "Yoruba", "Han"),
      target= "French")
french_adm$popdrop


basque_adm <- qpadm(data = f2_blocks,
      left = c("Loschbour", "LBK", "Yamnaya"),
      right = c("Mbuti", "Mota", "Dinka", "Yoruba", "Han"),
      target= "Basque")
basque_adm$popdrop

base_edges <- matrix(
  c("R",	"Mbuti",
    "R", "eAfr",
    "eAfr",	"Dinka",
    "eAfr",	"outAfrica",
    "outAfrica",	"Han",
    "outAfrica",	"Loschbour"),
  ncol=2,byrow = TRUE,
  dimnames=list(NULL, c("from","to")))

png("plots/graph_humans.png")
base_edges %>% edges_to_igraph() %>% plot_graph()
dev.off()

base_edges <- matrix(
  c("R",	"Mbuti",
    "R", "eAfr",
    "eAfr",	"Dinka",
    "eAfr",	"outAfrica",
    "outAfrica",	"Han",
    "outAfrica",	"Loschbour"),
  ncol=2,
  byrow = TRUE,
  dimnames=list(NULL, c("from","to")))

base_edges

base_igraph <- base_edges %>% edges_to_igraph()

is_valid(base_igraph)

png("plots/igraph_hamns.png")
base_igraph %>% plot_graph()
dev.off()

base_igraph %>% plotly_graph()

base_qpgraph <- qpgraph(data = f2_blocks, graph = base_igraph)

base_qpgraph$f3

base_qpgraph$f3 %>% filter(abs(z)>2)

png("plots/graph_withedges.png")
base_qpgraph$edges %>% plot_graph()
dev.off()

# base swapped 
base_edges_swapped <- matrix(
  c("R",	"Dinka",
    "R", "eAfr",
    "eAfr",	"Mbuti",
    "eAfr",	"outAfrica",
    "outAfrica",	"Han",
    "outAfrica",	"Loschbour"),
  ncol=2,
  byrow = TRUE,
  dimnames=list(NULL, c("from","to")))

base_edges_swapped

base_igraph <- base_edges_swapped %>% edges_to_igraph()

is_valid(base_igraph)

png("plots/igraph_swapped.png")
base_igraph %>% plot_graph()
dev.off()

base_igraph %>% plotly_graph()

base_qpgraph_swapped <- qpgraph(data = f2_blocks, graph = base_igraph)

base_qpgraph_swapped$f3

fits = qpgraph_resample_multi(f2_blocks, 
graphlist = list(base_qpgraph[[1]], base_qpgraph_swapped[[1]]), 
nboot = 100)
compare_fits(fits[[1]]$score_test, fits[[2]]$score_test)

# yamanaya graph
# base swapped 
base_edges_y <- matrix(
  c("R",	"Mbuti",
    "R", "eAfr",
    "eAfr",	"Dinka",
    "eAfr",	"outAfrica",
    "outAfrica",	"Han",
    "outAfrica", "wEurasian",
    "wEurasian", "Yamnaya",
    "wEurasian", "Loschbour"),
  ncol=2,
  byrow = TRUE,
  dimnames=list(NULL, c("from","to")))

base_edges_y

base_igraph <- base_edges_y %>% edges_to_igraph()

is_valid(base_igraph)

png("plots/igraph_yamnaya.png")
base_igraph %>% plot_graph()
dev.off()

base_igraph %>% plotly_graph()

base_qpgraph_y <- qpgraph(data = f2_blocks, graph = base_igraph)

base_qpgraph_y$f3

png("plots/fitted_y.png")
base_igraph %>% plot_graph(highlight_unidentifiable = TRUE)
dev.off()

fits = qpgraph_resample_multi(f2_blocks, 
graphlist = list(base_qpgraph[[1]], base_qpgraph_y[[1]]), 
nboot = 100)


yamnaya_edges <- matrix(
  c(
    "R",	"Mbuti",
    "R",	"eAfr",
    "eAfr", "Dinka",
    "eAfr", "outAfrica",
    "outAfrica", "Han",
    "outAfrica", "wEurasian",
    "wEurasian", "Yamnaya",
    "wEurasian", "Loschbour"),
  ncol = 2,
  byrow = TRUE,
  dimnames = list(NULL, c("from", "to")))
yamnaya_igraph <- yamnaya_edges %>% edges_to_igraph()
yamnaya_igraph %>% plot_graph()

lbk_extra_edges <- matrix(
  c(
    "R",	"Mbuti",
    "R",	"eAfr",
    "eAfr", "pBasalEurasian",
    "eAfr", "Dinka",
    "pBasalEurasian", "BasalEurasian",
    "pBasalEurasian","outAfrica",
    "outAfrica", "Han",
    "outAfrica","wEurasian",
    "wEurasian", "Yamnaya",
    "wEurasian", "pLoschbour",
    "pLoschbour", "Loschbour",
    "pLoschbour","WHG",
    "BasalEurasian", "pLBK",
    "WHG", "pLBK",
    "pLBK","LBK"),
  ncol = 2,
  byrow = TRUE,
  dimnames = list(NULL, c("from", "to")))
lbk_extra_igraph <- lbk_extra_edges %>% edges_to_igraph()
lbk_extra_igraph %>% plot_graph()

is_valid(lbk_extra_igraph)

pdf("plots/lbk_extra.pdf")
lbk_extra_igraph %>% plot_graph(highlight_unidentifiable = TRUE)
dev.off()

lbk_extra_qpgraph <- qpgraph(data = f2_blocks, graph = lbk_extra_igraph)

pdf("plots/fitted_lbk.pdf")
lbk_extra_qpgraph$edges %>% plot_graph()
dev.off()



lbk_extra_edges <- matrix(
  c(
    "R",	"Mbuti",
    "R",	"eAfr",
    "eAfr", "pBasalEurasian",
    "eAfr", "Dinka",
    "pBasalEurasian", "BasalEurasian",
    "pBasalEurasian","outAfrica",
    "outAfrica", "Han",
    "outAfrica","wEurasian",
    "wEurasian", "Yamnaya",
    "wEurasian", "pLoschbour",
    "pLoschbour", "Loschbour",
    "pLoschbour","WHG",
    "BasalEurasian", "pLBK",
    "WHG", "pLBK",
    "pLBK","LBK",
    "pLBK", "Sardinian",
    "Yamnaya", "Sardinian"),
  ncol = 2,
  byrow = TRUE,
  dimnames = list(NULL, c("from", "to")))
lbk_extra_igraph <- lbk_extra_edges %>% edges_to_igraph()
lbk_extra_igraph %>% plot_graph()

is_valid(lbk_extra_igraph)

pdf("plots/lbk_extra_sard.pdf")
lbk_extra_igraph %>% plot_graph(highlight_unidentifiable = TRUE)
dev.off()

lbk_extra_qpgraph <- qpgraph(data = f2_blocks, graph = lbk_extra_igraph)

pdf("plots/fitted_lbk_sard.pdf")
lbk_extra_qpgraph$edges %>% plot_graph()
dev.off()


lbk_extra_edges <- matrix(
  c(
    "R",	"Mbuti",
    "R",	"eAfr",
    "eAfr", "pBasalEurasian",
    "eAfr", "Dinka",
    "pBasalEurasian", "BasalEurasian",
    "pBasalEurasian","outAfrica",
    "outAfrica", "Han",
    "outAfrica","wEurasian",
    "wEurasian", "Yamnaya",
    "wEurasian", "pLoschbour",
    "pLoschbour", "Loschbour",
    "pLoschbour","WHG",
    "BasalEurasian", "pLBK",
    "WHG", "pLBK",
    "pLBK","LBK",
    "pLBK", "Sardinian"),
  ncol = 2,
  byrow = TRUE,
  dimnames = list(NULL, c("from", "to")))
lbk_extra_igraph <- lbk_extra_edges %>% edges_to_igraph()
lbk_extra_igraph %>% plot_graph()

is_valid(lbk_extra_igraph)

pdf("plots/lbk_extra_sard_noad.pdf")
lbk_extra_igraph %>% plot_graph(highlight_unidentifiable = TRUE)
dev.off()

lbk_extra_qpgraph <- qpgraph(data = f2_blocks, graph = lbk_extra_igraph)

pdf("plots/fitted_lbk_sard_noad.pdf")
lbk_extra_qpgraph$edges %>% plot_graph()
dev.off()
