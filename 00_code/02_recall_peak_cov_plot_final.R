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

# call peak
DefaultAssay(ORN)<-"ATAC"
Fragments(ORN)
peak<-CallPeaks(
  ORN,
  group.by = "seurat_clusters",
  macs2.path = "/public/home/nieyg/biosoft/conda/bin/macs2",
  broad = FALSE,
  format = "BED",
  fragment.tempdir = tempdir(),
  effective.genome.size = 2.7e+09,
  outdir="/data/R04/liyh526/project/01_Honeybee_IR/04_call_peak",
  combine.peaks=TRUE
)

macs2_counts <- FeatureMatrix(
  fragments = Fragments(ORN),
  features = peak,
  cells = colnames(ORN)
)     

# create a new assay using the MACS2 peak set and add it to the Seurat object
ORN[["peaks_ORN"]] <- CreateChromatinAssay(
  counts = macs2_counts,
  fragments = Fragments(ORN),
  annotation = Annotation(ORN)
)

saveRDS(ORN,"/data/R04/liyh526/project/01_Honeybee_IR/00_data/00_rds_data/Validate_OSN_merge_Orco_neg_neuron_tsne_recallpeak.rds")

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

# 绘制covplot
# 预处理
ORN <- readRDS("/data/R04/liyh526/project/01_Honeybee_IR/00_data/00_rds_data/Validate_OSN_merge_Orco_neg_neuron_tsne_recallpeak_remix.rds")
levels(ORN)
DefaultAssay(ORN) <- "peaks_ORN"
# first compute the GC content for each peak

Signac::Annotation(ORN)$tx_id <- Signac::Annotation(ORN)$gene_name
# link peaks to genes

ORN <- RegionStats(ORN, genome = BSgenome.Amel.HAv3.1.update.chemoreceptor)
######Visulize track######
plot_sorted_coverage <- function(seurat_obj, region_str, assay_name = "peaks_ORN", color_palette) {
  # 保存原始id顺序
  original_idents <- Idents(seurat_obj)
  original_levels <- levels(original_idents)
  
  # 创建颜色映射表（确保颜色与cluster绑定）
  color_map <- setNames(color_palette[1:length(original_levels)], original_levels)
  
  # 将区域字符串转换为GRanges对象
  region_gr <- Signac::StringToGRanges(region_str)
  
  # 获取fragment数据
  frags <- Fragments(seurat_obj[[assay_name]])
  if (length(frags) == 0) stop("No fragment data found in the assay")
  
  # 计算目标区域的覆盖深度
  cov_data <- Signac::CoveragePlot(
    object = seurat_obj,
    region = region_str,
    annotation = FALSE,
    peaks = FALSE,
    links = FALSE,
    return_data = TRUE  # 关键：返回原始数据而非绘图
  )
  
  # 提取覆盖深度数据
  coverage_df <- cov_data[[1]]$data
  
  # 计算每个cluster的最大覆盖深度
  group_max_coverage <- coverage_df %>%
    dplyr::group_by(group) %>%
    dplyr::summarise(max_coverage = max(coverage)) %>%
    ungroup()
  
  # 按最大覆盖深度降序排序
  sorted_clusters <- group_max_coverage %>%
    arrange(desc(max_coverage)) %>%
    pull(group)
  
  # 设置新的idents顺序
  new_idents <- factor(original_idents, levels = sorted_clusters)
  Idents(seurat_obj) <- new_idents
  
  # 绘制CoveragePlot
  p <- CoveragePlot(
    object = seurat_obj,
    region = region_str,
    annotation = TRUE,
    peaks = FALSE,
    links = FALSE
  ) & scale_fill_manual(values = color_map)  # 使用命名向量确保颜色正确绑定
  
  # 恢复原始idents顺序
  Idents(seurat_obj) <- original_idents
  seurat_obj@active.ident <- factor(Idents(seurat_obj), levels = original_levels)
  
  print(p)
}

# 使用新函数绘制所有图形
# 绘制特定IR基因的peak track (Fig3)
# Ir76a
p7 <- plot_sorted_coverage(ORN, "Group1-24543400-24546800", color_palette = mytsnecolors)
# Ir75f
p8 <- plot_sorted_coverage(ORN, "Group6-1669000-1675000", color_palette = mytsnecolors)
# Ir75n3
p9 <- plot_sorted_coverage(ORN, "Group15-7742000-7746000", color_palette = mytsnecolors)

pdf("/data/R04/liyh526/project/01_Honeybee_IR/06_final/2025_09_01/Fig3/Separate_peak/Fig3-Specific-IR-peaktrack.pdf", height=20, width=12)
p7 | p8 | p9
dev.off()

# 绘制Broad-IR基因的peak track (Fig2)
# Ir25a
p5 <- plot_sorted_coverage(ORN, "Group1-15746000-15750000", color_palette = mytsnecolors)
# Ir8a
p6 <- plot_sorted_coverage(ORN, "Group16-4598000-4605000", color_palette = mytsnecolors)

pdf("/data/R04/liyh526/project/01_Honeybee_IR/06_final/2025_09_01/Fig3/Separate_peak/Fig3-Broad-IR-peaktrack.pdf", height=20, width=8)
p5 | p6
dev.off()


# 修改后的函数：统一cluster顺序，0 cluster在最下部，优化多图显示
plot_unsorted_coverage <- function(seurat_obj, region_str, assay_name = "peaks_ORN", color_palette) {
  # 保存原始id顺序
  original_idents <- Idents(seurat_obj)
  original_levels <- levels(original_idents)
  
  # 确保0 cluster在最后（最下部）
  if ("0" %in% original_levels) {
    # 将0 cluster移到levels的最后
    non_zero_levels <- setdiff(original_levels, "0")
    fixed_levels <- c(non_zero_levels, "0")
    
    # 更新factor levels
    seurat_obj$Annotation <- factor(
      as.character(seurat_obj$Annotation),
      levels = fixed_levels
    )
    Idents(seurat_obj) <- seurat_obj$Annotation
  }
  
  # 创建颜色映射表（确保颜色与cluster绑定）
  color_map <- setNames(color_palette[1:length(levels(Idents(seurat_obj)))], levels(Idents(seurat_obj)))
  
  # 将区域字符串转换为GRanges对象
  region_gr <- Signac::StringToGRanges(region_str)
  
  # 获取fragment数据
  frags <- Fragments(seurat_obj[[assay_name]])
  if (length(frags) == 0) stop("No fragment data found in the assay")
  
  # 绘制CoveragePlot
  p <- CoveragePlot(
    object = seurat_obj,
    region = region_str,
    annotation = TRUE,
    peaks = FALSE,
    links = FALSE
  ) 
  
  # 应用颜色映射
  p <- p & scale_fill_manual(values = color_map)
  
  # 恢复原始idents顺序
  Idents(seurat_obj) <- original_idents
  seurat_obj@active.ident <- factor(Idents(seurat_obj), levels = original_levels)
  seurat_obj$Annotation <- factor(as.character(seurat_obj$Annotation), levels = original_levels)
  
  return(p)
}

# 辅助函数：移除Y轴标签和标题
remove_y_axis <- function(plot) {
  plot <- plot + 
    theme(
      axis.title.y = element_blank(),
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank()
    )
  return(plot)
}

# 绘制Orco-IR基因的peak track (Fig2)
# Or2
p1 <- plot_unsorted_coverage(ORN, "Group1-5723250-5725250", color_palette = mytsnecolors)
# Ir75d
p2 <- plot_unsorted_coverage(ORN, "Group11-16143000-16145000", color_palette = mytsnecolors)
# Ir76b
p3 <- plot_unsorted_coverage(ORN, "Group14-4996500-4997500", color_palette = mytsnecolors)
# Ir75n1
p4 <- plot_unsorted_coverage(ORN, "Group2-11821000-11824000", color_palette = mytsnecolors)

pdf("/data/R04/liyh526/project/01_Honeybee_IR/06_final/2025_09_01/Fig2/Fig2-Orco-IR-peaktrack.pdf", height=40, width=16)
p1 | p2 | p3 | p4
dev.off()




# 合并Track
# 读取数据
ORN <- readRDS("/data/R04/liyh526/project/01_Honeybee_IR/00_data/00_rds_data/Validate_OSN_merge_Orco_neg_neuron_tsne_recallpeak_remix.rds")

# 准备绘图数据
plot_data <- fread("/data/R04/liyh526/project/01_Honeybee_IR/06_final/Fig1/Fig1_Dotplot/Dotplot_with_accessibility_table.tsv", header = T)

# 从plot_data中提取每个基因在每个cluster的表达状态
get_expr_status_map <- function(plot_data) {
  expr_map <- list()
  
  # 遍历每个基因
  for (symbol in unique(plot_data$gene)) {
    gene_data <- plot_data[plot_data$gene == symbol, ]
    status_map <- setNames(gene_data$expr_status, as.character(gene_data$cluster))
    expr_map[[symbol]] <- status_map
  }
  
  return(expr_map)
}

# 根据表达状态合并cluster
merge_clusters_by_expr_status <- function(seurat_obj, gene_name, expr_status_map) {
  # 获取该基因的表达状态映射
  gene_status_map <- expr_status_map[[gene_name]]
  
  if (is.null(gene_status_map)) {
    warning(paste("No expression status data for gene", gene_name))
    return(seurat_obj)
  }
  
  # 创建新的分组列
  new_groups <- sapply(as.character(seurat_obj$Annotation), function(cluster) {
    if (cluster %in% names(gene_status_map)) {
      status <- gene_status_map[cluster]
      
      # 仅当状态为"No expression"时归类为无表达组
      if (status == "No expression") {
        paste0("NoExpr_", gene_name)
      } else {
        # 所有其他状态都归类为表达组
        paste0("Expr_", gene_name)
      }
    } else {
      # 如果cluster不在映射中，默认为无表达组
      paste0("NoExpr_", gene_name)
    }
  })
  
  # 添加到metadata
  seurat_obj@meta.data[[paste0("expr_group_", gene_name)]] <- factor(
    new_groups,
    levels = c(paste0("Expr_", gene_name), paste0("NoExpr_", gene_name))
  )
  
  return(seurat_obj)
}

# 修改Coverage Plot函数
plot_coverage_by_expr_status <- function(seurat_obj, region_str, gene_name, 
                                         assay_name = "peaks_ORN", 
                                         color_expr = "#dd6e42", 
                                         color_noexpr = "#e8dab2",
                                         window_size = 200) { # 添加平滑窗口参数
  
  group_col <- paste0("expr_group_", gene_name)
  
  if (!(group_col %in% colnames(seurat_obj@meta.data))) {
    stop(paste("Expression group column", group_col, "not found in metadata"))
  }
  
  Idents(seurat_obj) <- seurat_obj@meta.data[[group_col]]
  
  # 添加 window 参数进行平滑
  p <- CoveragePlot(
    object = seurat_obj,
    region = region_str,
    group.by = group_col,
    annotation = TRUE,
    peaks = FALSE,
    links = FALSE,
    window = window_size # 控制平滑程度，值越大越平滑
  ) + 
    ggtitle(paste(gene_name, "-", region_str)) +
    theme(plot.title = element_text(hjust = 0.5, size = 12))
  
  p <- p & scale_fill_manual(values = c(color_expr, color_noexpr))
  return(p)
}

# 从plot_data创建表达状态映射
expr_status_map <- get_expr_status_map(plot_data)

# 3. 为每个目标基因创建表达分组
genes_of_interest <- c("Ir25a", "Ir8a", 
                       "Ir75f", "Ir76a", "Ir75n3"
)

for (gene in genes_of_interest) {
  if (gene %in% names(expr_status_map)) {
    ORN <- merge_clusters_by_expr_status(ORN, gene, expr_status_map)
  } else {
    warning(paste("Gene", gene, "not found in expression status map. Skipping."))
  }
}

# 4. 绘制特定基因的Coverage Plot
# 定义基因与区域的对应关系
gene_region_map <- list(
  "Ir25a" = "Group1-15746000-15750000",
  "Ir8a" = "Group16-4598270-4605000",
  "Ir76a" = "Group1-24543400-24546800",
  "Ir75f" = "Group6-1669000-1675000",
  "Ir75n3" = "Group15-7742000-7746000"
)

# 绘制所有组合
# 在绘图循环中指定平滑窗口大小
plots <- list()
DefaultAssay(ORN) <- "peaks_ORN"
for (gene in names(gene_region_map)) {
  if (gene %in% genes_of_interest) {
    region <- gene_region_map[[gene]]
    
    # 设置平滑窗口
    window_size <- 700
    
    p <- plot_coverage_by_expr_status(ORN, region, gene,
                                      color_expr = "#dd6e42", 
                                      color_noexpr = "#e8dab2",
                                      window_size = window_size) # 传入平滑参数
    plots[[gene]] <- p
  }
}

# 5. 按基因类型分组组合图形
# 定义基因类型
broad_genes <- c("Ir25a", "Ir8a")
specific_genes <- c("Ir75f", "Ir76a", "Ir75n3")

# 创建分组绘图列表
grouped_plots <- list(
  broad = plots[names(plots) %in% broad_genes],
  specific = plots[names(plots) %in% specific_genes]
)

# 6. 保存结果
# 保存Broad基因
pdf("/data/R04/liyh526/project/01_Honeybee_IR/06_final/2025_09_01/Fig3/Broad_Genes_Coverage_Plots.pdf", height = 6, width = 12)
wrap_plots(grouped_plots$broad, ncol = 2)
dev.off()

# 保存Specific基因
pdf("/data/R04/liyh526/project/01_Honeybee_IR/06_final/2025_09_01/Fig3/Specific_Genes_Coverage_Plots.pdf", height = 6, width = 18)
wrap_plots(grouped_plots$specific, ncol = 3)
dev.off()