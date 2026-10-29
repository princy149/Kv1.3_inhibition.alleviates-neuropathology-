library(CellChat)
library(dplyr)
library(Seurat)
library(patchwork)
library(ComplexHeatmap)
library(RColorBrewer)
library(ggplot2)
library(compiler)
library(presto)
cellchat.1 <- readRDS("~/Desktop/NeuroChat.feb22/March14.2024/Cellchat/cellchat.1.rds")
cellchat.2 <- readRDS("~/Desktop/NeuroChat.feb22/March14.2024/Cellchat/cellchat.2.rds")
# use a subset of CellChatDB for cell-cell communication analysis
CellChatDB <- CellChatDB.mouse
CellChatDB
CellChatDB <- CellChatDB.mouse # use CellChatDB.human if running on mouse data
showDatabaseCategory(CellChatDB)
# use a subset of CellChatDB for cell-cell communication analysis
CellChatDB.Secreted.Cell.Cell.Contact  <- subsetDB(CellChatDB, search = "Cell-Cell Contact")
# set the used database in the object
cellchat.1@DB <- CellChatDB.Secreted.Cell.Cell.Contact
cellchat.2@DB <- CellChatDB.Secreted.Cell.Cell.Contact
cellchat.1 <- subsetData(cellchat.1) # subset the expression data of signaling genes
cellchat.2 <- subsetData(cellchat.2) # subset the expression data of signaling genes
#Compute the communication probability and infer cellular communication network
cellchat.1 <- identifyOverExpressedGenes(cellchat.1)
cellchat.1 <- identifyOverExpressedInteractions(cellchat.1)
cellchat.1 <- projectData(cellchat.1, PPI.mouse )
cellchat.1 <- computeCommunProb(cellchat.1, raw.use = TRUE)
# Filter out the cell-cell communication if there are only few number of cells in certain cell groups
cellchat.1 <- filterCommunication(cellchat.1, min.cells = 10)
cellchat.1 <- computeCommunProbPathway(cellchat.1)
cellchat.1 <- aggregateNet(cellchat.1)
#Compute the communication probability and infer cellular communication network
cellchat.2 <- identifyOverExpressedGenes(cellchat.2)
cellchat.2 <- identifyOverExpressedInteractions(cellchat.2)
cellchat.2 <- projectData(cellchat.2, PPI.mouse )
cellchat.2 <- computeCommunProb(cellchat.2, raw.use = TRUE)
# Filter out the cell-cell communication if there are only few number of cells in certain cell groups
cellchat.2  <- filterCommunication(cellchat.2, min.cells = 10)
cellchat.2  <- computeCommunProbPathway(cellchat.2)
cellchat.2  <- aggregateNet(cellchat.2)
object.list <- list(Pap1 = cellchat.1, Control = cellchat.2)
cellchat <- mergeCellChat(object.list, add.names = names(object.list))
cellchat
table(cellchat@meta$condition)
getwd()
saveRDS('Pap1.vs.Control.cell.cell.contact.merged.rds')
saveRDS(cellchat, file='Pap1.vs.Control.cell.cell.contact.merged.rds')
object.list
cellchat
save(object.list, file = "cellchat_object.list_Pap1.Control.cell.cell.contact.April.RData")
#Part I: Identify altered interactions and cell populations
#Compare the total number of interactions and interaction strength
ptm = Sys.time()
gg1 <- compareInteractions(cellchat, show.legend = F, group = c(1,2))
gg2 <- compareInteractions(cellchat, show.legend = F, group = c(1,2), measure = "weight")
gg1 + gg2
savehistory('Pap1.vs.Control.cell.cell.contact.r')
gg1 + gg2
dev.off()
##
num.link <- sapply(object.list, function(x) {rowSums(x@net$count) + colSums(x@net$count)-diag(x@net$count)})
weight.MinMax <- c(min(num.link), max(num.link)) # control the dot size in the different datasets
num.link
weight.MinMax
gg <- list()
for (i in 1:length(object.list)) {
gg[[i]] <- netAnalysis_signalingRole_scatter(object.list[[i]], title = names(object.list)[i], weight.MinMax = weight.MinMax)
}
gg1 <- rankNet(cellchat, mode = "comparison", measure = "weight", sources.use = NULL, targets.use = NULL, stacked = T, do.stat = TRUE)
gg2 <- rankNet(cellchat, mode = "comparison", measure = "weight", sources.use = NULL, targets.use = NULL, stacked = F, do.stat = TRUE)
gg1 + gg2
dev.off()
gg1 <- netVisual_bubble(cellchat, sources.use = c(1,5), targets.use = c(5:11),  comparison = c(1, 2), max.dataset = 2, title.name = "Increased signaling", angle.x = 45, remove.isolate = T)
gg2 <- netVisual_bubble(cellchat, sources.use = c(1,5), targets.use = c(5:11),  comparison = c(1, 2), max.dataset = 1, title.name = "Decreased signaling", angle.x = 45, remove.isolate = T)
gg1 + gg2
dev.off()
# Chord diagram
par(mfrow = c(1,2), xpd=TRUE)
netVisual_chord_gene(object.list[[2]], sources.use = 4, targets.use = c(5:11), slot.name = 'net', net = net.up, lab.cex = 0.8, small.gap = 3.5, title.name = paste0("Up-regulated signaling in ", names(object.list)[2]))
Pap1.signaling.Net.mat.count <- cellchat@net$Pap1$count
Pap1.signaling.Net.mat.count
Pap1.signaling.Net.mat.Strength <- cellchat@net$Pap1$weight
Pap1.signaling.Net.mat.Strength
write.csv(Pap1.signaling.Net.mat.count, file = 'Pap1.signaling.Net.mat.count.csv')
write.csv(Pap1.signaling.Net.mat.Strength, file = 'Pap1.signaling.Net.mat.Strength.csv')
Pap1.Cell.cell.contact.Net.mat.count <- cellchat@net$Pap1$count
Pap1.Cell.cell.contact.Net.mat.count
Pap1.Cell.cell.contact.Net.mat.Strength <- cellchat@net$Pap1$weight
Pap1.Cell.cell.contact.Net.mat.Strength
write.csv(Pap1.Cell.cell.contact.Net.mat.count, file = 'Pap1.Cell.cell.contact.Net.mat.count.csv')
write.csv(Pap1.Cell.cell.contact.Net.mat.Strength, file = 'Pap1.Cell.cell.contact.Net.mat.Strength.csv')
Control.Cell.cell.contact.Net.mat.count <- cellchat@net$Control$count
Control.Cell.cell.contact.Net.mat.count
Control.Cell.cell.contact.Net.mat.Strength <- cellchat@net$Control$weight
Control.Cell.cell.contact.Net.mat.Strength
write.csv(Control.Cell.cell.contact.Net.mat.count, file = 'Control.Cell.cell.contact.Net.mat.count.csv')
write.csv(Control.Cell.cell.contact.Net.mat.Strength, file = 'Control.Cell.cell.contact.Net.mat.Strength.csv')
savehistory('Pap1.vs.Control.cell.cell.contact.r')
# Signaling role analysis on the aggregated cell-cell communication network from all signaling pathways
ht1 <- netAnalysis_signalingRole_heatmap(cellchat, pattern = "outgoing")
x <- cellchat
x <- netAnalysis_computeCentrality(cellchat)
p1 <- netAnalysis_computeCentrality(cellchat.1)
cellchat.1 <- netAnalysis_computeCentrality(cellchat.1)
cellchat.2 <- netAnalysis_computeCentrality(cellchat.2)
object.list <- list(Pap1 = cellchat.1, Control = cellchat.2)
cellchat <- mergeCellChat(object.list, add.names = names(object.list))
cellchat
table(cellchat@meta$condition)
savehistory('Pap1.vs.Control.cell.cell.contact.r')
