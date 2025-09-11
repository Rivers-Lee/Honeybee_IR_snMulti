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

ORN <- readRDS("/data/R04/liyh526/project/01_Honeybee_IR/00_data/00_rds_data/Validate_OSN_merge_Orco_neg_neuron_tsne_recallpeak_remix.rds")
levels(ORN)

DefaultAssay(ORN)<-"RNA"

# 1. 确保 Annotation 列是因子类型
if (!is.factor(ORN$Annotation)) {
  ORN$Annotation <- factor(ORN$Annotation)
}

levels(ORN$Annotation)

# 创建LOC到基因符号的映射
loc_to_symbol <- c(
  "Or2" = "Or2",
  "LOC552552" = "Nmdar2-Like",
  "LOC726019" = "Ir75d",
  "LOC551704" = "Ir76b",
  "LOC102653640" = "Ir75n1",
  "LOC408645" = "KaiR1D",
  "LOC412818" = "Nmdar2",
  "LOC412949" = "Ir25a",
  "LOC412993" = "Ekar",
  "LOC409777" = "Ir8a",
  "LOC411220" = "GluRI-Like",
  "LOC102656088" = "Ir75f",
  "LOC725724" = "Ir76a",
  "LOC100576097" = "Ir75n3"
)

# 修改后的绘图函数
plot_sorted_violin <- function(seurat_obj, features, colors, gene_symbol_map, log_scale = TRUE) {
  # 创建结果列表存储每个基因的绘图
  plot_list <- list()
  
  # 获取当前分组信息
  clusters <- seurat_obj$Annotation
  original_levels <- levels(clusters)
  
  for (feature in features) {
    # 提取表达数据
    expr_data <- GetAssayData(seurat_obj, assay = "RNA", slot = "data")[feature, ]
    
    # 计算每个cluster的平均表达
    avg_expr <- tapply(expr_data, clusters, mean)
    
    # 按平均表达降序排序
    sorted_clusters <- names(avg_expr)[order(avg_expr, decreasing = TRUE)]
    
    # 创建临时对象避免修改原始数据
    temp_obj <- seurat_obj
    
    # 更新因子水平（按表达排序）
    temp_obj$Annotation <- factor(
      as.character(temp_obj$Annotation),
      levels = sorted_clusters
    )
    
    # 调整颜色顺序以匹配新的聚类顺序
    new_colors <- colors[match(sorted_clusters, original_levels)]
    
    # 生成小提琴图
    p <- VlnPlot(
      temp_obj,
      features = feature,
      group.by = "Annotation",
      cols = new_colors,
      log = log_scale,
      add.noise = F,
      stack = F, 
      pt.size = 0
    ) + 
      theme() +
      NoLegend()
    
    # 替换基因名（如果存在于映射表中）
    if (feature %in% names(gene_symbol_map)) {
      new_name <- gene_symbol_map[feature]
      p <- p + ggtitle(new_name)  # 修改标题
      p$plot_env$feature.plot <- new_name  # 修改内部存储的基因名
    }
    
    plot_list[[feature]] <- p
  }
  
  # 组合所有图形
  wrap_plots(plot_list, ncol = 1)
}

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

# 使用映射表绘制图形
pdf("/data/R04/liyh526/project/01_Honeybee_IR/06_final/Fig2/Orco_vln.pdf", width=22, height=14)
Orco <- c("Or2", "LOC726019", "LOC551704", "LOC102653640")
print(plot_sorted_violin(ORN, Orco, mytsnecolors, loc_to_symbol))
dev.off()

pdf("/data/R04/liyh526/project/01_Honeybee_IR/06_final/2025_09_01/Fig3/Broad_vln.pdf", width=22, height=4)
Broad <- c("LOC412949", "LOC409777")
print(plot_sorted_violin(ORN, Broad, mytsnecolors, loc_to_symbol))
dev.off()

pdf("/data/R04/liyh526/project/01_Honeybee_IR/06_final/2025_09_01/Fig3/Specifc_vln.pdf", width=22, height=10)
Specifc <- c("LOC725724", "LOC102656088", "LOC100576097")
print(plot_sorted_violin(ORN, Specifc, mytsnecolors, loc_to_symbol))
dev.off()


# 第一个函数：按表达量排序并添加UMI过滤
plot_sorted_violin_UMI <- function(seurat_obj, features, colors, gene_symbol_map, log_scale = TRUE) {
  # 创建结果列表存储每个基因的绘图
  plot_list <- list()
  
  # 获取当前分组信息
  clusters <- seurat_obj$Annotation
  original_levels <- levels(clusters)
  
  # 获取标准化后的表达数据和原始UMI计数数据
  norm_data <- GetAssayData(seurat_obj, assay = "RNA", slot = "data")[features, ]
  umi_data <- GetAssayData(seurat_obj, assay = "RNA", slot = "counts")[features, ]
  
  for (feature in features) {
    # 创建临时对象避免修改原始数据
    temp_obj <- seurat_obj
    
    # 复制表达数据
    temp_expr_data <- norm_data[feature, ]
    
    # 应用UMI过滤：如果UMI < 4，则表达设为0
    for (cluster in original_levels) {
      cluster_cells <- WhichCells(temp_obj, idents = cluster)
      total_umi <- sum(umi_data[feature, cluster_cells])
      
      if (total_umi < 4) {
        temp_expr_data[cluster_cells] <- 0
      }
    }
    
    # 更新临时对象的表达数据
    temp_obj[["RNA"]]@data[feature, ] <- temp_expr_data
    
    # 计算每个cluster的平均表达
    avg_expr <- tapply(temp_expr_data, clusters, mean)
    
    # 按平均表达降序排序
    sorted_clusters <- names(avg_expr)[order(avg_expr, decreasing = TRUE)]
    
    # 更新因子水平（按表达排序）
    temp_obj$Annotation <- factor(
      as.character(temp_obj$Annotation),
      levels = sorted_clusters
    )
    
    # 调整颜色顺序以匹配新的聚类顺序
    new_colors <- colors[match(sorted_clusters, original_levels)]
    
    # 生成小提琴图
    p <- VlnPlot(
      temp_obj,
      features = feature,
      group.by = "Annotation",
      cols = new_colors,
      log = log_scale,
      add.noise = F,
      stack = F, 
      pt.size = 0
    ) + 
      theme() +
      NoLegend()
    
    # 替换基因名（如果存在于映射表中）
    if (feature %in% names(gene_symbol_map)) {
      new_name <- gene_symbol_map[feature]
      p <- p + ggtitle(new_name)
      p$plot_env$feature.plot <- new_name
    }
    
    plot_list[[feature]] <- p
  }
  
  # 组合所有图形
  wrap_plots(plot_list, ncol = 1)
}

plot_unsorted_violin_UMI <- function(seurat_obj, features, colors, gene_symbol_map, log_scale = TRUE) {
  # 创建结果列表存储每个基因的绘图
  plot_list <- list()
  
  # 获取当前分组信息
  clusters <- seurat_obj$Annotation
  original_levels <- levels(clusters)
  
  # 获取标准化后的表达数据和原始UMI计数数据
  norm_data <- GetAssayData(seurat_obj, assay = "RNA", slot = "data")[features, ]
  umi_data <- GetAssayData(seurat_obj, assay = "RNA", slot = "counts")[features, ]
  
  # 找出在所有基因中都无表达的cluster
  # 首先确定每个基因在每个cluster中的表达状态
  expr_status <- list()
  for (feature in features) {
    expr_status[[feature]] <- sapply(original_levels, function(cluster) {
      cluster_cells <- WhichCells(seurat_obj, idents = cluster)
      total_umi <- sum(umi_data[feature, cluster_cells])
      return(total_umi >= 4) # 只有UMI >= 4才认为有表达
    })
  }
  
  # 找出在所有基因中都无表达的cluster
  all_zero_clusters <- original_levels[sapply(original_levels, function(cluster) {
    all(sapply(expr_status, function(x) !x[cluster]))
  })]
  
  # 找出至少在一个基因中有表达的cluster
  expressed_clusters <- setdiff(original_levels, all_zero_clusters)
  
  # 确保0 cluster在所有无表达cluster的最右侧
  if ("0" %in% all_zero_clusters) {
    all_zero_clusters <- c(setdiff(all_zero_clusters, "0"), "0")
  }
  
  # 确定统一的顺序：先是有表达的cluster（按原始顺序），然后是无表达的cluster
  unified_order <- c(expressed_clusters, all_zero_clusters)
  
  # 打印调试信息
  cat("统一顺序:", paste(unified_order, collapse = ", "), "\n")
  cat("无表达cluster:", paste(all_zero_clusters, collapse = ", "), "\n")
  
  for (feature in features) {
    # 创建临时对象避免修改原始数据
    temp_obj <- seurat_obj
    
    # 复制表达数据
    temp_expr_data <- norm_data[feature, ]
    
    # 应用UMI过滤：如果UMI < 4，则表达设为0
    for (cluster in original_levels) {
      cluster_cells <- WhichCells(temp_obj, idents = cluster)
      total_umi <- sum(umi_data[feature, cluster_cells])
      
      if (total_umi < 4) {
        temp_expr_data[cluster_cells] <- 0
      }
    }
    
    # 更新临时对象的表达数据
    temp_obj[["RNA"]]@data[feature, ] <- temp_expr_data
    
    # 使用统一的cluster顺序
    temp_obj$Annotation <- factor(
      as.character(temp_obj$Annotation),
      levels = unified_order
    )
    
    # 调整颜色顺序以匹配统一的聚类顺序
    new_colors <- colors[match(unified_order, original_levels)]
    
    # 生成小提琴图
    p <- VlnPlot(
      temp_obj,
      features = feature,
      group.by = "Annotation",
      cols = new_colors,
      log = log_scale,
      add.noise = FALSE,
      stack = FALSE, 
      pt.size = 0
    ) + 
      theme() +
      NoLegend()
    
    # 替换基因名（如果存在于映射表中）
    if (feature %in% names(gene_symbol_map)) {
      new_name <- gene_symbol_map[feature]
      p <- p + ggtitle(new_name)
      p$plot_env$feature.plot <- new_name
    }
    
    # 如果不是最后一个图，移除X轴标签
    if (feature != features[length(features)]) {
      p <- p + theme(axis.title.x = element_blank(), 
                     axis.text.x = element_blank(),
                     axis.ticks.x = element_blank())
    } else {
      # 最后一个图，确保X轴标签可读
      p <- p + theme(axis.text.x = element_text(angle = 45, hjust = 1))
    }
    
    plot_list[[feature]] <- p
  }
  
  # 组合所有图形
  wrap_plots(plot_list, ncol = 1)
}


# 使用映射表绘制图形
pdf("/data/R04/liyh526/project/01_Honeybee_IR/06_final/2025_09_01/Fig2/Orco_vln.pdf", width=22, height=8)
Orco <- c("Or2", "LOC726019", "LOC551704", "LOC102653640")
print(plot_unsorted_violin_UMI(ORN, Orco, mytsnecolors, loc_to_symbol))
dev.off()

pdf("/data/R04/liyh526/project/01_Honeybee_IR/06_final/2025_09_01/Fig3/Broad_vln.pdf", width=22, height=4)
Broad <- c("LOC412949", "LOC409777")
print(plot_sorted_violin_UMI(ORN, Broad, mytsnecolors, loc_to_symbol))
dev.off()

pdf("/data/R04/liyh526/project/01_Honeybee_IR/06_final/2025_09_01/Fig3/Specifc_vln.pdf", width=22, height=6)
Specifc <- c("LOC725724", "LOC102656088", "LOC100576097")
print(plot_sorted_violin_UMI(ORN, Specifc, mytsnecolors, loc_to_symbol))
dev.off()








# IR feature plot
library(ggplot2)
library(tidyverse)
DefaultAssay(ORN) <- "RNA"

# 更新 Seurat 对象中的 tsne 降维结果
pdf("/data/R04/liyh526/project/01_Honeybee_IR/06_final/FigS1/IR_iGluR_All_Neuron_featureplot.pdf")

FeaturePlot(ORN,features = "LOC409777",min.cutoff = 0,max.cutoff = 2, pt.size = 0.7, reduction = "tsne.rna", order = TRUE)

FeaturePlot(ORN,features = "LOC412949",min.cutoff = 0,max.cutoff = 2, pt.size = 0.7, reduction = "tsne.rna", order = TRUE)

FeaturePlot(ORN,features = "LOC551704",min.cutoff = 0,max.cutoff = 2, pt.size = 0.7, reduction = "tsne.rna", order = TRUE)

FeaturePlot(ORN,features = "Or2",min.cutoff = 0,max.cutoff = 2, pt.size = 0.7, reduction = "tsne.rna")

FeaturePlot(ORN,features = "LOC726019",min.cutoff = 0,max.cutoff = 2, pt.size = 0.7, reduction = "tsne.rna", order = TRUE)

FeaturePlot(ORN,features = "LOC552552",min.cutoff = 0,max.cutoff = 2, pt.size = 0.7, reduction = "tsne.rna", order = TRUE)

FeaturePlot(ORN,features = "LOC408645",min.cutoff = 0,max.cutoff = 2, pt.size = 0.7, reduction = "tsne.rna", order = TRUE)

FeaturePlot(ORN,features = "LOC102653640",min.cutoff = 0,max.cutoff = 2, pt.size = 0.7, reduction = "tsne.rna", order = TRUE)

FeaturePlot(ORN,features = "LOC412818",min.cutoff = 0,max.cutoff = 2, pt.size = 0.7, reduction = "tsne.rna", order = TRUE)

FeaturePlot(ORN,features = "LOC411220",min.cutoff = 0,max.cutoff = 2, pt.size = 0.7, reduction = "tsne.rna", order = TRUE)

FeaturePlot(ORN,features = "LOC412993",min.cutoff = 0,max.cutoff = 2, pt.size = 0.7, reduction = "tsne.rna", order = TRUE)

FeaturePlot(ORN,features = "LOC100576097",min.cutoff = 0,max.cutoff = 2, pt.size = 0.7, reduction = "tsne.rna", order = TRUE)

FeaturePlot(ORN,features = "LOC102656088",min.cutoff = 0,max.cutoff = 2, pt.size = 0.7, reduction = "tsne.rna", order = TRUE)

FeaturePlot(ORN,features = "LOC725724",min.cutoff = 0,max.cutoff = 2, pt.size = 0.7, reduction = "tsne.rna", order = TRUE)

dev.off()

pdf("/data/R04/liyh526/project/01_Honeybee_IR/03_fig/00_All_IR_iGluR_All_Neuron_Dimplot.pdf")
DimPlot(ORN, cols = mytsnecolors, label = F, repel = F, pt.size=0.7, reduction = "tsne.rna", group.by = "Annotation") + NoLegend()
dev.off()
