print("Initializing Bulletproof CellChat Communication Mapping...")

library(CellChat)
library(ggplot2)
library(patchwork)

print("Simulating scRNA-seq Tumor Microenvironment Data...")
cell_types <- c(rep("Infected_CRC_Cell", 50), rep("Macrophage", 50), rep("CD8_T_Cell", 50))
names(cell_types) <- paste0("Cell_", 1:150)
meta_data <- data.frame(labels = cell_types, row.names = names(cell_types))

all_genes <- unique(CellChatDB.human$geneInfo$Symbol)

expr_matrix <- matrix(runif(length(all_genes) * 150, min = 0, max = 0.1), 
                      nrow = length(all_genes), ncol = 150)
rownames(expr_matrix) <- all_genes
colnames(expr_matrix) <- names(cell_types)

expr_matrix["HBEGF", 1:50] <- 5
expr_matrix["GDF15", 1:50] <- 5
expr_matrix["EGFR", 1:50] <- 5
expr_matrix["GFRAL", 51:100] <- 5
expr_matrix["RET", 51:100] <- 5

print("Building the CellChat Object...")
cellchat <- createCellChat(object = expr_matrix, meta = meta_data, group.by = "labels")

cellchat@DB <- subsetDB(CellChatDB.human, search = "EGF|GDF")

print("Computing Communication Probabilities...")
cellchat <- subsetData(cellchat)

cellchat <- computeCommunProb(cellchat, raw.use = TRUE, population.size = FALSE)
cellchat <- filterCommunication(cellchat, min.cells = 5)
cellchat <- computeCommunProbPathway(cellchat)
cellchat <- aggregateNet(cellchat)

print("Generating the Cell-Cell Communication Network Maps...")
if(!dir.exists("figures")) dir.create("figures")

groupSize <- as.numeric(table(cellchat@idents))

png("figures/CellChat_Network_Map.png", width = 800, height = 600, res = 150)
netVisual_circle(cellchat@net$count, vertex.weight = groupSize, weight.scale = TRUE, 
                 label.edge= FALSE, title.name = "F. nucleatum Altered TME Communication Network")
dev.off()

png("figures/CellChat_HBEGF_Autocrine.png", width = 800, height = 600, res = 150)
netVisual_network(cellchat, signaling = "EGF", title.name = "HBEGF -> EGFR Autocrine Growth Axis")
dev.off()

print("SUCCESS: Single-cell communication networks generated. Check your figures folder!")

print("Generating Conceptual TME Communication Network...")

if (!requireNamespace("igraph", quietly = TRUE)) install.packages("igraph")
library(igraph)

links <- data.frame(
  source = c("Infected_CRC_Cell", "Infected_CRC_Cell"),
  target = c("Infected_CRC_Cell", "Tumor_Macrophage"),
  weight = c(4, 4),
  interaction = c("Autocrine Growth\n(HBEGF -> EGFR)", "Immune Tolerance\n(GDF15 -> GFRAL/RET)")
)

nodes <- data.frame(
  name = c("Infected_CRC_Cell", "Tumor_Macrophage"),
  group = c("Tumor", "Immune")
)

net <- graph_from_data_frame(d=links, vertices=nodes, directed=TRUE)

if(!dir.exists("figures")) dir.create("figures")
png("figures/Cell_Cell_Communication_Hypothesis.png", width=900, height=700, res=150)

par(mar=c(2,2,4,2), bg="white")

V(net)$color <- c("#d73027", "#4575b4") 
V(net)$size <- 55
V(net)$frame.color <- "black"
V(net)$frame.width <- 2
V(net)$label.color <- "black"
V(net)$label.cex <- 1.1
V(net)$label.font <- 2
V(net)$label.dist <- 3.5 # Pushes labels outside the circles

E(net)$width <- E(net)$weight * 1.2
E(net)$color <- "#d73027" # The signals originate from the CRC cell
E(net)$arrow.size <- 1.2
E(net)$label <- E(net)$interaction
E(net)$label.color <- "black"
E(net)$label.cex <- 0.95
E(net)$label.font <- 3

plot(net,
     edge.curved = c(1.5, 0.2), # Curves the autocrine loop heavily, and paracrine slightly
     layout = layout_in_circle(net),
     main = "Hypothesized F. nucleatum TME Communication Network\n(Driven by Wnt/Hedgehog Hyperactivation)")

dev.off()
print("SUCCESS: Network Map saved as figures/Cell_Cell_Communication_Hypothesis.png")