# 读取必须的包
library(Matrix)
packageVersion("Matrix") # 确认版本
library(BSgenome.Amel.HAv3.1.update.chemoreceptor)
library(ggplot2)
library(dplyr)
library(patchwork)
library(irlba)
library(Signac)
library(Seurat)
library(SeuratObject)

ORN <- readRDS("/data/R04/liyh526/project/01_Honeybee_IR/00_data/00_rds_data/Unsupervised_ORN_remove_nopower_modify_the_tsne_recall_peak_latest.rds")
levels(ORN)
All_Neuron <- readRDS("/data/R04/liyh526/project/01_Honeybee_IR/00_data/00_rds_data/WNN_Neuron_integrated.rds")
levels(All_Neuron)

# 注释转换
head(ORN$cell_group)
head(All_Neuron$Annotation)

# 将Annotation转换为字符向量以避免因子问题
All_Neuron$Annotation <- as.character(All_Neuron$Annotation)

# 获取所有Orco+Neuron的细胞ID
orc_pos_cells <- names(All_Neuron$Annotation)[All_Neuron$Annotation == "Orco+Neuron"]

# 检查是否存在Orco+Neuron细胞不在ORN中（可选，给出警告）
missing_cells <- setdiff(orc_pos_cells, names(ORN$cell_group))
if (length(missing_cells) > 0) {
  warning(paste("Some Orco+Neuron cells not found in ORN$cell_group:", length(missing_cells)))
}

# 提取共同细胞ID并获取对应的cell_group（转换为字符）
common_cells <- intersect(orc_pos_cells, names(ORN$cell_group))
new_annotations <- as.character(ORN$cell_group[common_cells])

# 更新All_Neuron中的注释
All_Neuron$Annotation[common_cells] <- new_annotations

Idents(All_Neuron) <- All_Neuron$Annotation

mytsnecolors <- c(
  '#53A85F', '#F1BB72', '#F3B1A0', '#D6E7A3', '#57C3F3', '#476D87',
  '#E95C59', '#E59CC4', '#AB3282', '#23452F', '#BD956A', '#8C549C', '#585658',
  '#9FA3A8', '#E0D4CA', '#5F3D69', '#58A4C3', '#AA9A59', '#E63863', '#E39A35', 
  '#C1E6F3', '#6778AE', '#B53E2B', '#712820', '#DCC1DD', '#CCE0F5', '#625D9E', 
  '#68A180', '#3A6963', '#968175', '#161853', '#FF9999', '#344CB7', '#FFCC1D', 
  '#116530', '#678983', '#A19882', '#FFBCBC', '#24A19C', '#FF9A76', "#8DD3C7",
  "#FFFFB3", "#BEBADA", "#FB8072", "#80B1D3", "#FDB462", "#B3DE69", "#FCCDE5", 
  "#D9D9D9", "#BC80BD", "#CCEBC5", "#FFED6F", "#E41A1C", "#377EB8", "#4DAF4A", 
  "#F1BB72", "#AA9A59", "#A65628", "#F781BF",
  
  # 新增的5个颜色 (保持视觉区分度)
  "#9C27B0",  # 鲜艳紫罗兰 (补充紫色系)
  "#00BCD4",  # 亮青色 (补充蓝绿色系)
  "#795548",  # 深咖啡色 (补充棕色系)
  "#CDDC39",  # 酸柠檬黄 (补充亮黄色)
  "#7B1FA2"   # 深紫罗兰 (增强紫色深度)
)

DimPlot(All_Neuron, reduction = "wnn.umap", cols=mytsnecolors, label = TRUE)

# 直接提取子集
All_valided_Neuron <- subset(All_Neuron,idents=c(seq(1:59), "Orco-Neuron"))

DimPlot(All_valided_Neuron, reduction = "wnn.umap", cols = mytsnecolors, label = TRUE)
DimPlot(All_valided_Neuron, group.by = "orig.ident", reduction = "wnn.umap", cols=mytsnecolors, label = TRUE)


# 重新聚类分群 
All_valided_Neuron.list <- SplitObject(All_valided_Neuron, split.by = "orig.ident");
for (i in 1:length(All_valided_Neuron.list)) {
  DefaultAssay(All_valided_Neuron.list[[i]]) <- "raw_RNA"
  All_valided_Neuron.list[[i]] <- SCTransform(All_valided_Neuron.list[[i]],return.only.var.genes = FALSE,variable.features.n=10000,verbose = FALSE)
}
All_valided_Neuron.features <- SelectIntegrationFeatures(object.list = All_valided_Neuron.list, nfeatures = 4000);
All_valided_Neuron.list <- PrepSCTIntegration(object.list = All_valided_Neuron.list, anchor.features = All_valided_Neuron.features)
#integrate RNA using rpca
All_valided_Neuron.list <- lapply(
  X = All_valided_Neuron.list,
  FUN = ScaleData,
  features = All_valided_Neuron.features,
  verbose = FALSE
)
All_valided_Neuron.list <- lapply(
  X = All_valided_Neuron.list,
  FUN = RunPCA,
  features = All_valided_Neuron.features,
  verbose = FALSE
)

lapply(X = All_valided_Neuron.list, FUN = Key)

integration_anchors <- FindIntegrationAnchors(
  object.list = All_valided_Neuron.list,
  normalization.method = "SCT",
  anchor.features = All_valided_Neuron.features,
  dims = 1:30,
  reduction = "rpca",
  k.anchor = 20,
)
all_features <- lapply(All_valided_Neuron.list, row.names) %>% Reduce(intersect,.) 

All_valided_Neuron <- IntegrateData(
  anchorset = integration_anchors,
  normalization.method = "SCT",
  new.assay.name = "integratedRNA_ORN",
  features.to.integrate = all_features,
  dims = 1:30
)

ORN <- All_valided_Neuron
#run LSI on new seurat object with integrated RNA assay
DefaultAssay(ORN) <- "ATAC"
Fragments(ORN)
ORN <- RunTFIDF(ORN)
ORN <- FindTopFeatures(ORN, min.cutoff = "q25")
ORN <- RunSVD(ORN,n=25)
#integrate embeddings and output new object to prevent overwriting integrated RNA
#I can't quite remember why I did this the first time, but it was the work around that I needed
ORN_atac <- IntegrateEmbeddings(
  anchorset = integration_anchors,
  new.reduction.name = "integratedLSI_ORN",
  reductions = ORN@reductions$lsi
)

#copy integrated LSI from duplicate seurat object to original object
ORN@reductions$integratedLSI_ORN <- ORN_atac@reductions$integratedLSI_ORN
#####done integrate ATAC and RNA ################

# RNA analysis
DefaultAssay(ORN) <- "integratedRNA_ORN"
ORN <- ScaleData(ORN,features=rownames(ORN))
ORN <- RunPCA(ORN)
ORN <- RunTSNE(ORN,dims = 1:50, min.dist = 0.001, check_duplicates = FALSE, reduction.name = 'tsne.rna', reduction.key = 'rnaTSNE_')

ORN <- FindNeighbors(object = ORN, reduction = 'pca', dims = 1:50)
ORN <- FindClusters(object = ORN, verbose = FALSE,  
                    resolution = c(3, 4, 5), 
                    algorithm = 3,
                    random.seed = 345)
table(ORN$integratedRNA_ORN_snn_res.3)
table(ORN$integratedRNA_ORN_snn_res.4)
table(ORN$integratedRNA_ORN_snn_res.5)
table(ORN$integratedRNA_ORN_snn_res.6)

ORN$seurat_clusters <- ORN$integratedRNA_ORN_snn_res.5

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

# reorder the level
Idents(ORN) <- ORN$seurat_clusters
table(Idents(ORN))
names(ORN@graphs)

# Find SubCluster
ORN <- FindSubCluster(ORN, cluster = "15", graph.name = "integratedRNA_ORN_snn", algorithm = 3, resolution = 0.5)
table(ORN$sub.cluster)
Idents(ORN) <- "sub.cluster"
ORN <- FindSubCluster(ORN, cluster = "26", graph.name = "integratedRNA_ORN_snn", algorithm = 3, resolution = 0.5)
table(ORN$sub.cluster)
Idents(ORN) <- "sub.cluster"
DimPlot(ORN, cols = mytsnecolors, label = TRUE, repel = TRUE, pt.size=0.3, reduction = "tsne.rna", group.by = "sub.cluster")

ORN <- RenameIdents(ORN, '15_0' = '15', '15_1' = '34', '26_0' = '26', '26_1' = '35')
DimPlot(ORN, cols = mytsnecolors, label = TRUE, repel = TRUE, pt.size=0.3, reduction = "tsne.rna")
ORN$seurat_clusters <- Idents(ORN)

ORN$orig.ident<-factor(ORN$orig.ident,levels=c("NE","Nurse","Forager"))

pdf("/data/R04/liyh526/project/01_Honeybee_IR/03_fig/01_cluster/raw_ORN_cluster.pdf",width=6,height=5)
DimPlot(ORN, cols = mytsnecolors, label = TRUE, repel = TRUE, pt.size=0.5, reduction = "tsne.rna", group.by = "seurat_clusters")
DimPlot(ORN, cols = mytsnecolors, label = TRUE, repel = TRUE, pt.size=0.5, reduction = "tsne.rna", group.by = "orig.ident")
DimPlot(ORN, cols = mytsnecolors, label = F, repel = F, pt.size=0.5, reduction = "tsne.rna", group.by = "Annotation")
dev.off()

saveRDS(ORN,"/data/R04/liyh526/project/01_Honeybee_IR/00_data/00_rds_data/Validate_OSN_merge_Orco_neg_neuron_tsne.rds")

# 修改reduction
ORN <- readRDS("/data/R04/liyh526/project/01_Honeybee_IR/00_data/00_rds_data/Validate_OSN_merge_Orco_neg_neuron_tsne_recallpeak.rds")
ORN_only <- readRDS("/data/R04/liyh526/project/01_Honeybee_IR/00_data/00_rds_data/Unsupervised_ORN_remove_nopower_modify_the_tsne_recall_peak_latest.rds")

head(ORN@reductions[["tsne.rna"]]@cell.embeddings)
dim(ORN@reductions[["tsne.rna"]]@cell.embeddings)
head(ORN_only@reductions[["tsne.rna"]]@cell.embeddings)
dim(ORN_only@reductions[["tsne.rna"]]@cell.embeddings)

# 提取两个对象的嵌入矩阵
orn_embeddings <- ORN@reductions[["tsne.rna"]]@cell.embeddings
orn_only_embeddings <- ORN_only@reductions[["tsne.rna"]]@cell.embeddings

# 获取细胞barcode
orn_barcodes <- rownames(orn_embeddings)
orn_only_barcodes <- rownames(orn_only_embeddings)

# 找出独有细胞和共有细胞
unique_to_orn <- setdiff(orn_barcodes, orn_only_barcodes)  # ORN中独有的55个细胞
common_barcodes <- intersect(orn_barcodes, orn_only_barcodes)  # 共有的4120个细胞

# 提取非重叠细胞的嵌入（55个细胞）并改动
unique_embeddings <- orn_embeddings[unique_to_orn, ]
unique_embeddings[,1] <- unique_embeddings[,1] + 30
unique_embeddings[,2] <- unique_embeddings[,2] - 88

# 用ORN_only的嵌入替换共同细胞的嵌入
orn_embeddings[common_barcodes, ] <- orn_only_embeddings[common_barcodes, ]
orn_embeddings[unique_to_orn, ] <- unique_embeddings

# 将更新后的嵌入矩阵重新赋值给ORN对象
ORN@reductions[["tsne.rna"]]@cell.embeddings <- orn_embeddings

DimPlot(ORN, cols = mytsnecolors, label = F, repel = F, pt.size=0.5, reduction = "tsne.rna", group.by = "Annotation")
DimPlot(ORN_only, cols = mytsnecolors, label = F, repel = F, pt.size=0.5, reduction = "tsne.rna", group.by = "cell_group")

levels(ORN$Annotation)
levels(ORN$Annotation)[levels(ORN$Annotation) == "Orco-Neuron"] <- "0"
Idents(ORN) <- ORN$Annotation

saveRDS(ORN,"/data/R04/liyh526/project/01_Honeybee_IR/00_data/00_rds_data/Validate_OSN_merge_Orco_neg_neuron_tsne_recallpeak_remix.rds")



# 绘制cluster数目的柱状图
# 1. 绘制DimPlot（假设seurat_object是您的Seurat对象）
p <- DimPlot(ORN, cols = mytsnecolors, label = F, repel = F, pt.size=0.5, reduction = "tsne.rna", group.by = "Annotation")

# 2. 提取ggplot构建数据
pb <- ggplot2::ggplot_build(p)

# 3. 获取颜色映射表
color_data <- unique(pb$data[[1]][c("group", "colour")])
color_data <- color_data[order(color_data$group), ]

# 4. 关联cluster名称
clusters <- levels(p$data$Annotation)  # 获取有序cluster名称
color_table <- data.frame(
  cluster = clusters,
  color_hex = color_data$colour
)

# 查看颜色表
print(color_table)

# 3. 创建频数表并转换Annotation类型
freq_table <- table(ORN$Annotation)
df <- data.frame(
  Annotation = names(freq_table) %>% as.character(),  # 转为字符型
  Frequency = as.numeric(freq_table)
) %>% 
  arrange(desc(Frequency))  # 按频数降序

# 4. 匹配颜色
df <- left_join(df, color_table, by = c("Annotation" = "cluster"))

# 按频数排序并转换为因子（确保绘图顺序）
df <- df %>%
  arrange(desc(Frequency)) %>%
  mutate(Annotation = factor(Annotation, levels = Annotation))

# 5. 创建美观的柱状图
pdf("/data/R04/liyh526/project/01_Honeybee_IR/06_final/Fig1/Cluster_number.pdf",width=4,height=10)
ggplot(df, aes(x = Annotation, y = Frequency, fill = Annotation)) +
  geom_bar(stat = "identity", width = 0.8, color = "white", linewidth = 0.3) +
  geom_text(aes(label = Frequency), 
            vjust = 0.5,
            hjust = -0.25,
            size = 3.5, 
            color = "black",
            fontface = "bold") +
  scale_fill_manual(values = setNames(df$color_hex, df$Annotation)) +
  labs(title = "",
       x = "",
       y = "Number of Cells",
       caption = paste("Total cells:", sum(df$Frequency))) +
  theme_classic(base_size = 12) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
    axis.title = element_text(face = "bold"),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.y = element_blank(),
    legend.position = "none",
    plot.caption = element_text(hjust = 1, face = "italic", size = 9)
  ) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1))) + # 为顶部标签留空间
  coord_flip()
dev.off()