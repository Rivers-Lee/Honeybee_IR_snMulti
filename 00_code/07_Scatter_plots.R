# FigS2D  SCATTER PLOTS (Or2 Vs potential Orco)
# 读取必须的包
library(Matrix)
packageVersion("Matrix") # 确认版本
library(ggplot2)
library(dplyr)
library(patchwork)
library(irlba)
library(Signac)
library(Seurat)
library(SeuratObject)

mytsnecolors <- c(
  "#E41A1C", "#E95C59", "#FB8072", "#FF9999", "#FFBCBC", "#F3B1A0", "#FF9A76",
  "#E63863", "#F781BF", "#E59CC4", "#9C27B0", "#AB3282", "#8C549C", "#5F3D69", 
  "#BC80BD", "#BEBADA", "#FFB897", "#E4D457", "#DDA03C", "#DCC1DD", "#FCCDE5", 
  "#B53E2B", "#712820", "#A65628", "#FF7F00", "#E39A35", "#F1BB72", "#FDB462", 
  "#FFCC1D", "#AA9A59", "#BD956A", "#968175", "#A19882", "#678983", "#D9D9D9",
  "#9FA3A8", "#585658", "#E0D4CA", "#B3DE69", "#D6E7A3", "#4DAF4A", "#53A85F", 
  "#68A180", "#116530", "#23452F", "#3A6963", "#24A19C", "#8DD3C7", "#CCEBC5",
  "#57C3F3", "#C1E6F3", "#CCE0F5", "#58A4C3", "#80B1D3", "#377EB8", "#344CB7", 
  "#6778AE", "#625D9E", "#161853", "#476D87"
)

ORN <- readRDS("/data/R04/liyh526/project/01_Honeybee_IR/00_data/00_rds_data/Validate_OSN_merge_Orco_neg_neuron_tsne_recallpeak_remix.rds")
levels(ORN)
DefaultAssay(ORN) <- "RNA"
scatterColors <- c('#A06CB4', '#DF6C78', '#911A2E', '#CD9139', '#B4B4B6', '#21918c', '#3b528b', '#440154')
### Plotting scatter plot: 1. single-cell level, 2. cluster level
ORN$Or2_UMIs <- ORN@assays$SCT@counts['Or2',]
ORN$LOC552552_UMIs <- ORN@assays$SCT@counts['LOC552552',]
ORN$LOC726019_UMIs <- ORN@assays$SCT@counts['LOC726019',]
ORN$LOC551704_UMIs <- ORN@assays$SCT@counts['LOC551704',]
ORN$LOC102653640_UMIs <- ORN@assays$SCT@counts['LOC102653640',]
ORN$LOC408645_UMIs <- ORN@assays$SCT@counts['LOC408645',]
ORN$LOC412818_UMIs <- ORN@assays$SCT@counts['LOC412818',]
ORN$LOC412949_UMIs <- ORN@assays$SCT@counts['LOC412949',]
ORN$LOC409777_UMIs <- ORN@assays$SCT@counts['LOC409777',]
ORN$LOC102656088_UMIs <- ORN@assays$SCT@counts['LOC102656088',]
ORN$LOC725724_UMIs <- ORN@assays$SCT@counts['LOC725724',]
ORN$LOC412993_UMIs <- ORN@assays$SCT@counts['LOC412993',]
ORN$LOC100576097_UMIs <- ORN@assays$SCT@counts['LOC100576097',]

ORN$Or2_Exp <- ORN@assays$RNA@data['Or2',]
ORN$LOC552552_Exp <- ORN@assays$RNA@data['LOC552552',]
ORN$LOC726019_Exp <- ORN@assays$RNA@data['LOC726019',]
ORN$LOC551704_Exp <- ORN@assays$RNA@data['LOC551704',]
ORN$LOC102653640_Exp <- ORN@assays$SCT@data['LOC102653640',]
ORN$LOC408645_Exp <- ORN@assays$SCT@data['LOC408645',]
ORN$LOC412818_Exp <- ORN@assays$SCT@data['LOC412818',]
ORN$LOC412949_Exp <- ORN@assays$SCT@data['LOC412949',]
ORN$LOC409777_Exp <- ORN@assays$SCT@data['LOC409777',]
ORN$LOC102656088_Exp <- ORN@assays$SCT@data['LOC102656088',]
ORN$LOC725724_Exp <- ORN@assays$SCT@data['LOC725724',]
ORN$LOC412993_Exp <- ORN@assays$SCT@data['LOC412993',]
ORN$LOC100576097_Exp <- ORN@assays$SCT@data['LOC100576097',]

colors_list <- mytsnecolors

library(tidyr)
library(ggrepel)
library(cowplot)
# 2. in cluster level
gene_order_final <- c("Or2",
                      "LOC552552","LOC726019","LOC551704","LOC102653640","LOC408645","LOC412818","LOC412949", # Universal
                      "LOC409777","LOC412993", # Broad
                      "LOC102656088","LOC725724","LOC100576097" # Specific
)

DefaultAssay(ORN)<-"SCT"
p<-DotPlot(ORN,features = gene_order_final) +  xlab('') + ylab('') + theme(axis.text.x = element_text(angle=90, hjust=1, vjust=.5, size = 9)) 
dotplot_data<-p$data;

# LOC552552
plot_avgExp.clusters.df <- 
  dotplot_data %>%
  filter(features.plot %in% c('Or2', 'LOC552552')) %>%
  select(-pct.exp, -avg.exp.scaled) %>%
  pivot_wider(names_from = features.plot, values_from = avg.exp) %>% as.data.frame() 

pmain = plot_avgExp.clusters.df %>%
  ggplot(aes(Or2, LOC552552, color=id)) + 
  geom_point() +
  scale_color_manual(values=colors_list)+
  geom_text_repel(aes(label=id), hjust=0) + 
  theme(legend.position = "none") +
  scale_x_log10() + scale_y_log10()+
  xlab('Average expression of Or2') + ylab('Average expression of LOC552552')
xdens <- axis_canvas(pmain, axis = "x")+
  geom_density(data = plot_avgExp.clusters.df, aes(x = Or2), fill=scatterColors[1])+
  scale_x_log10()
ydens <- axis_canvas(pmain, axis = "y", coord_flip = TRUE)+
  geom_density(data = plot_avgExp.clusters.df, aes(x = LOC552552), fill=scatterColors[2]) + 
  scale_x_log10()+
  coord_flip()
p1 <- insert_xaxis_grob(pmain, xdens, grid::unit(.2, "null"), position = "top")
p2 <- insert_yaxis_grob(p1, ydens, grid::unit(.2, "null"), position = "right")
p3 <- ggdraw(p2)
p3

filename <- paste0("/data/R04/liyh526/project/01_Honeybee_IR/03_fig/neuronFigures_Or2-v-LOC552552_clusters.pdf")
ggsave(filename, limitsize = FALSE, units = "px", width = 2000, height = 2000,p3)

# LOC726019
plot_avgExp.clusters.df <- 
  dotplot_data %>%
  filter(features.plot %in% c('Or2', 'LOC726019')) %>%
  select(-pct.exp, -avg.exp.scaled) %>%
  pivot_wider(names_from = features.plot, values_from = avg.exp) %>% as.data.frame() 
pmain=  plot_avgExp.clusters.df %>%
  ggplot(aes(Or2, LOC726019, color=id)) + 
  geom_point() +scale_color_manual(values=colors_list)+
  geom_text_repel(aes(label=id), hjust=0) + 
  theme(legend.position = "none") +
  scale_x_log10() + scale_y_log10()+
  xlab('Average expression of Or2') + ylab('Average expression of LOC726019')
xdens <- axis_canvas(pmain, axis = "x")+
  geom_density(data = plot_avgExp.clusters.df, aes(x = Or2), fill=scatterColors[1])+
  scale_x_log10()
ydens <- axis_canvas(pmain, axis = "y", coord_flip = TRUE)+
  geom_density(data = plot_avgExp.clusters.df, aes(x = LOC726019), fill=scatterColors[2]) + 
  scale_x_log10()+
  coord_flip()
p1 <- insert_xaxis_grob(pmain, xdens, grid::unit(.2, "null"), position = "top")
p2 <- insert_yaxis_grob(p1, ydens, grid::unit(.2, "null"), position = "right")
p3 <- ggdraw(p2)
p3

filename <- paste0("/data/R04/liyh526/project/01_Honeybee_IR/03_fig/neuronFigures_Or2-v-LOC726019_clusters.pdf")
ggsave(filename, limitsize = FALSE, units = "px", width = 2000, height = 2000,p3)

# LOC551704
plot_avgExp.clusters.df <- 
  dotplot_data %>%
  filter(features.plot %in% c('Or2', 'LOC551704')) %>%
  select(-pct.exp, -avg.exp.scaled) %>%
  pivot_wider(names_from = features.plot, values_from = avg.exp) %>% as.data.frame() 

pmain=  plot_avgExp.clusters.df %>%
  ggplot(aes(Or2, LOC551704, color=id)) + 
  geom_point() +scale_color_manual(values=colors_list)+
  geom_text_repel(aes(label=id), hjust=0) + 
  theme(legend.position = "none") +
  scale_x_log10() + scale_y_log10()+
  xlab('Average expression of Or2') + ylab('Average expression of LOC551704')
xdens <- axis_canvas(pmain, axis = "x")+
  geom_density(data = plot_avgExp.clusters.df, aes(x = Or2), fill=scatterColors[1])+
  scale_x_log10()
ydens <- axis_canvas(pmain, axis = "y", coord_flip = TRUE)+
  geom_density(data = plot_avgExp.clusters.df, aes(x = LOC551704), fill=scatterColors[2]) + 
  scale_x_log10()+
  coord_flip()
p1 <- insert_xaxis_grob(pmain, xdens, grid::unit(.2, "null"), position = "top")
p2 <- insert_yaxis_grob(p1, ydens, grid::unit(.2, "null"), position = "right")
p3 <- ggdraw(p2)
p3

filename <- paste0("/data/R04/liyh526/project/01_Honeybee_IR/03_fig/neuronFigures_Or2-v-LOC551704_clusters.pdf")
ggsave(filename, limitsize = FALSE, units = "px", width = 2000, height = 2000,p3)

# LOC102653640
plot_avgExp.clusters.df <- 
  dotplot_data %>%
  filter(features.plot %in% c('Or2', 'LOC102653640')) %>%
  select(-pct.exp, -avg.exp.scaled) %>%
  pivot_wider(names_from = features.plot, values_from = avg.exp) %>% as.data.frame() 

pmain=  plot_avgExp.clusters.df %>%
  ggplot(aes(Or2, LOC102653640, color=id)) + 
  geom_point() +scale_color_manual(values=colors_list)+
  geom_text_repel(aes(label=id), hjust=0) + 
  theme(legend.position = "none") +
  scale_x_log10() + scale_y_log10()+
  xlab('Average expression of Or2') + ylab('Average expression of LOC102653640')
xdens <- axis_canvas(pmain, axis = "x")+
  geom_density(data = plot_avgExp.clusters.df, aes(x = Or2), fill=scatterColors[1])+
  scale_x_log10()
ydens <- axis_canvas(pmain, axis = "y", coord_flip = TRUE)+
  geom_density(data = plot_avgExp.clusters.df, aes(x = LOC102653640), fill=scatterColors[2]) + 
  scale_x_log10()+
  coord_flip()
p1 <- insert_xaxis_grob(pmain, xdens, grid::unit(.2, "null"), position = "top")
p2 <- insert_yaxis_grob(p1, ydens, grid::unit(.2, "null"), position = "right")
p3 <- ggdraw(p2)
p3

filename <- paste0("/data/R04/liyh526/project/01_Honeybee_IR/03_fig/neuronFigures_Or2-v-LOC102653640_clusters.pdf")
ggsave(filename, limitsize = FALSE, units = "px", width = 2000, height = 2000,p3)

# LOC408645
plot_avgExp.clusters.df <- 
  dotplot_data %>%
  filter(features.plot %in% c('Or2', 'LOC408645')) %>%
  select(-pct.exp, -avg.exp.scaled) %>%
  pivot_wider(names_from = features.plot, values_from = avg.exp) %>% as.data.frame() 

pmain=  plot_avgExp.clusters.df %>%
  ggplot(aes(Or2, LOC408645, color=id)) + 
  geom_point() +scale_color_manual(values=colors_list)+
  geom_text_repel(aes(label=id), hjust=0) + 
  theme(legend.position = "none") +
  scale_x_log10() + scale_y_log10()+
  xlab('Average expression of Or2') + ylab('Average expression of LOC408645')
xdens <- axis_canvas(pmain, axis = "x")+
  geom_density(data = plot_avgExp.clusters.df, aes(x = Or2), fill=scatterColors[1])+
  scale_x_log10()
ydens <- axis_canvas(pmain, axis = "y", coord_flip = TRUE)+
  geom_density(data = plot_avgExp.clusters.df, aes(x = LOC408645), fill=scatterColors[2]) + 
  scale_x_log10()+
  coord_flip()
p1 <- insert_xaxis_grob(pmain, xdens, grid::unit(.2, "null"), position = "top")
p2 <- insert_yaxis_grob(p1, ydens, grid::unit(.2, "null"), position = "right")
p3 <- ggdraw(p2)
p3

filename <- paste0("/data/R04/liyh526/project/01_Honeybee_IR/03_fig/neuronFigures_Or2-v-LOC408645_clusters.pdf")
ggsave(filename, limitsize = FALSE, units = "px", width = 2000, height = 2000,p3)

# LOC412818
plot_avgExp.clusters.df <- 
  dotplot_data %>%
  filter(features.plot %in% c('Or2', 'LOC412818')) %>%
  select(-pct.exp, -avg.exp.scaled) %>%
  pivot_wider(names_from = features.plot, values_from = avg.exp) %>% as.data.frame() 

pmain=  plot_avgExp.clusters.df %>%
  ggplot(aes(Or2, LOC412818, color=id)) + 
  geom_point() +scale_color_manual(values=colors_list)+
  geom_text_repel(aes(label=id), hjust=0) + 
  theme(legend.position = "none") +
  scale_x_log10() + scale_y_log10()+
  xlab('Average expression of Or2') + ylab('Average expression of LOC412818')
xdens <- axis_canvas(pmain, axis = "x")+
  geom_density(data = plot_avgExp.clusters.df, aes(x = Or2), fill=scatterColors[1])+
  scale_x_log10()
ydens <- axis_canvas(pmain, axis = "y", coord_flip = TRUE)+
  geom_density(data = plot_avgExp.clusters.df, aes(x = LOC412818), fill=scatterColors[2]) + 
  scale_x_log10()+
  coord_flip()
p1 <- insert_xaxis_grob(pmain, xdens, grid::unit(.2, "null"), position = "top")
p2 <- insert_yaxis_grob(p1, ydens, grid::unit(.2, "null"), position = "right")
p3 <- ggdraw(p2)
p3

filename <- paste0("/data/R04/liyh526/project/01_Honeybee_IR/03_fig/neuronFigures_Or2-v-LOC412818_clusters.pdf")
ggsave(filename, limitsize = FALSE, units = "px", width = 2000, height = 2000,p3)

# LOC412949
plot_avgExp.clusters.df <- 
  dotplot_data %>%
  filter(features.plot %in% c('Or2', 'LOC412949')) %>%
  select(-pct.exp, -avg.exp.scaled) %>%
  pivot_wider(names_from = features.plot, values_from = avg.exp) %>% as.data.frame() 

pmain=  plot_avgExp.clusters.df %>%
  ggplot(aes(Or2, LOC412949, color=id)) + 
  geom_point() +scale_color_manual(values=colors_list)+
  geom_text_repel(aes(label=id), hjust=0) + 
  theme(legend.position = "none") +
  scale_x_log10() + scale_y_log10()+
  xlab('Average expression of Or2') + ylab('Average expression of LOC412949')
xdens <- axis_canvas(pmain, axis = "x")+
  geom_density(data = plot_avgExp.clusters.df, aes(x = Or2), fill=scatterColors[1])+
  scale_x_log10()
ydens <- axis_canvas(pmain, axis = "y", coord_flip = TRUE)+
  geom_density(data = plot_avgExp.clusters.df, aes(x = LOC412949), fill=scatterColors[2]) + 
  scale_x_log10()+
  coord_flip()
p1 <- insert_xaxis_grob(pmain, xdens, grid::unit(.2, "null"), position = "top")
p2 <- insert_yaxis_grob(p1, ydens, grid::unit(.2, "null"), position = "right")
p3 <- ggdraw(p2)
p3

filename <- paste0("/data/R04/liyh526/project/01_Honeybee_IR/03_fig/neuronFigures_Or2-v-LOC412949_clusters.pdf")
ggsave(filename, limitsize = FALSE, units = "px", width = 2000, height = 2000,p3)

# LOC409777
plot_avgExp.clusters.df <- 
  dotplot_data %>%
  filter(features.plot %in% c('Or2', 'LOC409777')) %>%
  select(-pct.exp, -avg.exp.scaled) %>%
  pivot_wider(names_from = features.plot, values_from = avg.exp) %>% as.data.frame() 

pmain=  plot_avgExp.clusters.df %>%
  ggplot(aes(Or2, LOC409777, color=id)) + 
  geom_point() +scale_color_manual(values=colors_list)+
  geom_text_repel(aes(label=id), hjust=0) + 
  theme(legend.position = "none") +
  scale_x_log10() + scale_y_log10()+
  xlab('Average expression of Or2') + ylab('Average expression of LOC409777')
xdens <- axis_canvas(pmain, axis = "x")+
  geom_density(data = plot_avgExp.clusters.df, aes(x = Or2), fill=scatterColors[1])+
  scale_x_log10()
ydens <- axis_canvas(pmain, axis = "y", coord_flip = TRUE)+
  geom_density(data = plot_avgExp.clusters.df, aes(x = LOC409777), fill=scatterColors[2]) + 
  scale_x_log10()+
  coord_flip()
p1 <- insert_xaxis_grob(pmain, xdens, grid::unit(.2, "null"), position = "top")
p2 <- insert_yaxis_grob(p1, ydens, grid::unit(.2, "null"), position = "right")
p3 <- ggdraw(p2)
p3

filename <- paste0("/data/R04/liyh526/project/01_Honeybee_IR/03_fig/neuronFigures_Or2-v-LOC409777_clusters.pdf")
ggsave(filename, limitsize = FALSE, units = "px", width = 2000, height = 2000,p3)

# LOC102656088
plot_avgExp.clusters.df <- 
  dotplot_data %>%
  filter(features.plot %in% c('Or2', 'LOC102656088')) %>%
  select(-pct.exp, -avg.exp.scaled) %>%
  pivot_wider(names_from = features.plot, values_from = avg.exp) %>% as.data.frame() 

pmain=  plot_avgExp.clusters.df %>%
  ggplot(aes(Or2, LOC102656088, color=id)) + 
  geom_point() +scale_color_manual(values=colors_list)+
  geom_text_repel(aes(label=id), hjust=0) + 
  theme(legend.position = "none") +
  scale_x_log10() + scale_y_log10()+
  xlab('Average expression of Or2') + ylab('Average expression of LOC102656088')
xdens <- axis_canvas(pmain, axis = "x")+
  geom_density(data = plot_avgExp.clusters.df, aes(x = Or2), fill=scatterColors[1])+
  scale_x_log10()
ydens <- axis_canvas(pmain, axis = "y", coord_flip = TRUE)+
  geom_density(data = plot_avgExp.clusters.df, aes(x = LOC102656088), fill=scatterColors[2]) + 
  scale_x_log10()+
  coord_flip()
p1 <- insert_xaxis_grob(pmain, xdens, grid::unit(.2, "null"), position = "top")
p2 <- insert_yaxis_grob(p1, ydens, grid::unit(.2, "null"), position = "right")
p3 <- ggdraw(p2)
p3

filename <- paste0("/data/R04/liyh526/project/01_Honeybee_IR/03_fig/neuronFigures_Or2-v-LOC102656088_clusters.pdf")
ggsave(filename, limitsize = FALSE, units = "px", width = 2000, height = 2000,p3)

# LOC725724
plot_avgExp.clusters.df <- 
  dotplot_data %>%
  filter(features.plot %in% c('Or2', 'LOC725724')) %>%
  select(-pct.exp, -avg.exp.scaled) %>%
  pivot_wider(names_from = features.plot, values_from = avg.exp) %>% as.data.frame() 

pmain=  plot_avgExp.clusters.df %>%
  ggplot(aes(Or2, LOC725724, color=id)) + 
  geom_point() +scale_color_manual(values=colors_list)+
  geom_text_repel(aes(label=id), hjust=0) + 
  theme(legend.position = "none") +
  scale_x_log10() + scale_y_log10()+
  xlab('Average expression of Or2') + ylab('Average expression of LOC725724')
xdens <- axis_canvas(pmain, axis = "x")+
  geom_density(data = plot_avgExp.clusters.df, aes(x = Or2), fill=scatterColors[1])+
  scale_x_log10()
ydens <- axis_canvas(pmain, axis = "y", coord_flip = TRUE)+
  geom_density(data = plot_avgExp.clusters.df, aes(x = LOC725724), fill=scatterColors[2]) + 
  scale_x_log10()+
  coord_flip()
p1 <- insert_xaxis_grob(pmain, xdens, grid::unit(.2, "null"), position = "top")
p2 <- insert_yaxis_grob(p1, ydens, grid::unit(.2, "null"), position = "right")
p3 <- ggdraw(p2)
p3

filename <- paste0("/data/R04/liyh526/project/01_Honeybee_IR/03_fig/neuronFigures_Or2-v-LOC725724_clusters.pdf")
ggsave(filename, limitsize = FALSE, units = "px", width = 2000, height = 2000,p3)

# LOC412993
plot_avgExp.clusters.df <- 
  dotplot_data %>%
  filter(features.plot %in% c('Or2', 'LOC412993')) %>%
  select(-pct.exp, -avg.exp.scaled) %>%
  pivot_wider(names_from = features.plot, values_from = avg.exp) %>% as.data.frame() 

pmain=  plot_avgExp.clusters.df %>%
  ggplot(aes(Or2, LOC412993, color=id)) + 
  geom_point() +scale_color_manual(values=colors_list)+
  geom_text_repel(aes(label=id), hjust=0) + 
  theme(legend.position = "none") +
  scale_x_log10() + scale_y_log10()+
  xlab('Average expression of Or2') + ylab('Average expression of LOC412993')
xdens <- axis_canvas(pmain, axis = "x")+
  geom_density(data = plot_avgExp.clusters.df, aes(x = Or2), fill=scatterColors[1])+
  scale_x_log10()
ydens <- axis_canvas(pmain, axis = "y", coord_flip = TRUE)+
  geom_density(data = plot_avgExp.clusters.df, aes(x = LOC412993), fill=scatterColors[2]) + 
  scale_x_log10()+
  coord_flip()
p1 <- insert_xaxis_grob(pmain, xdens, grid::unit(.2, "null"), position = "top")
p2 <- insert_yaxis_grob(p1, ydens, grid::unit(.2, "null"), position = "right")
p3 <- ggdraw(p2)
p3

filename <- paste0("/data/R04/liyh526/project/01_Honeybee_IR/03_fig/neuronFigures_Or2-v-LOC412993_clusters.pdf")
ggsave(filename, limitsize = FALSE, units = "px", width = 2000, height = 2000,p3)

# LOC100576097
plot_avgExp.clusters.df <- 
  dotplot_data %>%
  filter(features.plot %in% c('Or2', 'LOC100576097')) %>%
  select(-pct.exp, -avg.exp.scaled) %>%
  pivot_wider(names_from = features.plot, values_from = avg.exp) %>% as.data.frame() 

pmain=  plot_avgExp.clusters.df %>%
  ggplot(aes(Or2, LOC100576097, color=id)) + 
  geom_point() +scale_color_manual(values=colors_list)+
  geom_text_repel(aes(label=id), hjust=0) + 
  theme(legend.position = "none") +
  scale_x_log10() + scale_y_log10()+
  xlab('Average expression of Or2') + ylab('Average expression of LOC100576097')
xdens <- axis_canvas(pmain, axis = "x")+
  geom_density(data = plot_avgExp.clusters.df, aes(x = Or2), fill=scatterColors[1])+
  scale_x_log10()
ydens <- axis_canvas(pmain, axis = "y", coord_flip = TRUE)+
  geom_density(data = plot_avgExp.clusters.df, aes(x = LOC100576097), fill=scatterColors[2]) + 
  scale_x_log10()+
  coord_flip()
p1 <- insert_xaxis_grob(pmain, xdens, grid::unit(.2, "null"), position = "top")
p2 <- insert_yaxis_grob(p1, ydens, grid::unit(.2, "null"), position = "right")
p3 <- ggdraw(p2)
p3

filename <- paste0("/data/R04/liyh526/project/01_Honeybee_IR/03_fig/neuronFigures_Or2-v-LOC100576097_clusters.pdf")
ggsave(filename, limitsize = FALSE, units = "px", width = 2000, height = 2000,p3)
