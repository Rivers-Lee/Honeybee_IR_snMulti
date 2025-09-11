library(Seurat)
library(ggplot2)
library(dplyr)
library(cowplot)
library(tidyr)
library(pheatmap)
library(lsa)
library(GenomicRanges)
library(Signac)

# 基因列表
IR_gene_list <- c("LOC412949", "LOC408645", "LOC410788", "LOC100577496", 
                  "LOC102653640", "LOC412993", "LOC411220", "Nmdar1",
                  "LOC552552", "LOC412818", "LOC726019", "LOC551704",
                  "LOC100576097", "LOC409777", "LOC102656088", "LOC410623", 
                  "LOC727346", "LOC725724", "LOC102653800", "LOC411050")


ORN <- readRDS("/data/R04/liyh526/project/01_Honeybee_IR/00_data/00_rds_data/Validate_OSN_merge_Orco_neg_neuron_tsne_recallpeak_remix.rds")
DefaultAssay(ORN) <- "RNA"
expr_threshold <- 0.1

# 计算每个cluster的平均表达量
avg_expr <- AverageExpression(ORN, 
                              assays = "RNA",
                              features = IR_gene_list,
                              group.by = "Annotation")[[1]]

Idents(ORN) <- ORN$Annotation
# 计算每个cluster中表达基因的细胞比例
prop_matrix <- matrix(0, nrow = length(IR_gene_list), ncol = length(unique(ORN$Annotation)))
rownames(prop_matrix) <- IR_gene_list
colnames(prop_matrix) <- unique(ORN$Annotation)

exp_data <- GetAssayData(ORN, slot = "data")[IR_gene_list, ]

for (cluster in unique(ORN$Annotation)) {
  cells_in_cluster <- WhichCells(ORN, idents = cluster)
  cluster_data <- exp_data[, cells_in_cluster, drop = FALSE]
  
  # 计算该cluster中每个基因的表达比例
  prop_exp <- rowMeans(cluster_data > expr_threshold)
  prop_matrix[, cluster] <- prop_exp * 100  # 转换为百分比
}

# 筛选满足条件的基因
filtered_genes <- IR_gene_list[
  sapply(IR_gene_list, function(gene) {
    any(avg_expr[gene, ] > expr_threshold) |   # 表达量满足
      any(prop_matrix[gene, ] > 10)            # 或表达比例满足
  })
]

# 输出结果
print(paste0("筛选后基因数量: ", length(filtered_genes)))
print(filtered_genes)

# 剔除iGluR
filtered_genes <- c(
  "LOC726019", "LOC551704", "LOC102653640", "LOC412949",
  "LOC409777", "LOC725724", "LOC102656088", "LOC100576097"
)

# 创建LOC到基因符号的映射
loc_to_symbol <- c(
  "LOC726019" = "Ir75d",
  "LOC551704" = "Ir76b",
  "LOC102653640" = "Ir75n1",
  "LOC412949" = "Ir25a",
  "LOC409777" = "Ir8a",
  "LOC102656088" = "Ir75f",
  "LOC725724" = "Ir76a",
  "LOC100576097" = "Ir75n3"
)
symbol_to_loc <- setNames(names(loc_to_symbol), loc_to_symbol)

# 获取最终基因顺序（使用基因符号）
gene_order_final <- c(
  "Ir75d", "Ir76b", 
  "Ir75n1", "Ir25a", 
  "Ir8a", "Ir76a", 
  "Ir75f", "Ir75n3"
)

# 提取最终表达
avg_expr <- AverageExpression(ORN, 
                              assays = "RNA",
                              features = filtered_genes,  # 使用原始LOC名称提取数据
                              group.by = "Annotation")[[1]]

# 将表达矩阵的行名转换为基因符号
rownames(avg_expr) <- loc_to_symbol[rownames(avg_expr)]


# 启动子区域映射
gene_region_map <- list(
  "LOC726019" = "Group11-16143000-16145000",
  "LOC551704" = "Group14-4996500-4997500",
  "LOC102653640" = "Group2-11821000-11824000",
  "LOC412949" = "Group1-15746000-15750000",
  "LOC409777" = "Group16-4598000-4600000",
  "LOC725724" = "Group1-24545600-24546800",
  "LOC102656088" = "Group6-1669000-1673500",
  "LOC100576097" = "Group15-7744500-7745500"
)

# 创建启动子区域的GRanges对象
promoter_list <- strsplit(unlist(gene_region_map), "-")

# 启动子区域定义
promoter_gr <- GRanges(
  seqnames = sapply(promoter_list, `[`, 1),
  ranges = IRanges(
    start = as.numeric(sapply(promoter_list, `[`, 2)),
    end = as.numeric(sapply(promoter_list, `[`, 3))
  ))
names(promoter_gr) <- names(gene_region_map)

# 提取区域信号矩阵 (使用Signac高效方法)
region_signal <- FeatureMatrix(
  fragments = Fragments(ORN[["peaks_ORN"]]),
  features = promoter_gr,
  cells = colnames(ORN),
  process_n = 2000,  # 并行处理
  sep = c("-", "-")
)

# 深度归一化 (消除技术偏差)
depth_factor <- ORN$nCount_peaks_ORN / median(ORN$nCount_peaks_ORN)
norm_signal <- t(t(region_signal) / depth_factor)

# 计算cluster水平最大信号
cluster_max_signal <- matrix(0, nrow = nrow(norm_signal), 
                             ncol = length(unique(ORN$Annotation)))
rownames(cluster_max_signal) <- names(refined_promoter_gr)
colnames(cluster_max_signal) <- unique(ORN$Annotation)

for (cluster in unique(ORN$Annotation)) {
  cluster_cells <- colnames(ORN)[ORN$Annotation == cluster]
  if (length(cluster_cells) > 0) {
    cluster_data <- norm_signal[, cluster_cells, drop = FALSE]
    cluster_max_signal[, cluster] <- apply(cluster_data, 1, max)
  }
}

# 动态阈值设定 (基于基因特异性)
promoter_prop_matrix <- matrix(0, nrow = nrow(cluster_max_signal), 
                               ncol = ncol(cluster_max_signal))
rownames(promoter_prop_matrix) <- rownames(cluster_max_signal)
colnames(promoter_prop_matrix) <- colnames(cluster_max_signal)

for (gene in rownames(cluster_max_signal)) {
  max_val <- max(cluster_max_signal[gene, ])
  threshold <- max_val / 6  # 使用最高信号的1/3作为阈值
  promoter_prop_matrix[gene, ] <- ifelse(
    cluster_max_signal[gene, ] >= threshold, 1, 0
  )
}

# 转换为基因符号
rownames(promoter_prop_matrix) <- loc_to_symbol[rownames(promoter_prop_matrix)]

# 检查结果
print(promoter_prop_matrix)
  
# 准备绘图数据
gene_order_final <- c("Ir75d", "Ir76b", "Ir75n1", "Ir25a", "Ir8a", "Ir76a", "Ir75f", "Ir75n3")
  
plot_data <- expand.grid(
  gene = gene_order_final,
  cluster = unique(ORN$Annotation)
) %>%
  mutate(
    avg_expr = sapply(1:nrow(.), function(i) {
      avg_expr[as.character(.$gene[i]), as.character(.$cluster[i])]
    }),
    prop = sapply(1:nrow(.), function(i) {
      # 确保prop_matrix的行名也已转换为基因符号
      prop_matrix_symbols <- prop_matrix
      rownames(prop_matrix_symbols) <- loc_to_symbol[rownames(prop_matrix)]
      prop_matrix_symbols[as.character(.$gene[i]), as.character(.$cluster[i])]
    }),
    promoter_prop = sapply(1:nrow(.), function(i) {
      promoter_prop_matrix[as.character(.$gene[i]), as.character(.$cluster[i])]
    }),
    # 定义四种表达状态
    expr_status = case_when(
      # 通用规则
      avg_expr > expr_threshold & prop > 10 ~ "High expression",
      avg_expr > expr_threshold & prop <= 10 ~ "Lack expression",
      avg_expr <= expr_threshold & prop > 10 ~ "Low expression",
      
      # 特殊规则：只对Ir25a和Ir8a应用Accessible only
      gene %in% c("Ir25a", "Ir8a") & 
        avg_expr <= expr_threshold & 
        prop <= 10 & 
        promoter_prop == 1 ~ "Accessible only",
      
      # 其他情况都是No expression
      TRUE ~ "No expression"
    )
  )

fwrite(plot_data, "/data/R04/liyh526/project/01_Honeybee_IR/06_final/Fig1/Fig1_Dotplot/Dotplot_with_accessibility_table.tsv", 
       sep = "\t", quote = FALSE)

# 排序基因和cluster
# 计算每个基因的"No expression"数量并排序基因
no_expr_count <- sapply(gene_order_final, function(g) {
  sum(plot_data$gene == g & plot_data$expr_status == "No expression")
})
gene_order_by_no <- names(sort(no_expr_count))
  
# 创建表达状态矩阵
expr_status_matrix <- dcast(plot_data, cluster ~ gene, value.var = "expr_status")
rownames(expr_status_matrix) <- expr_status_matrix$cluster
expr_status_matrix$cluster <- NULL
  
cluster_order <- c("0", "5", "47", "43", "3", "15", 
                   "16", "49", "9", "34", "45", "32", 
                   "40", "37", "39", "18", "27", "29", 
                   "42", "25", "50", "6", "31", "46", 
                   "51", "58", "55", "1", "38", "57", 
                   "33", "20","30", "19", "17", "28", 
                   "24", "56", "59", "54", "13", "53", 
                   "11", "35", "44", "52", "26", "14", 
                   "23", "48", "8", "12", "21", "10", 
                   "4", "7", "22", "41", "36", "2")
  
# 定义基因分组
universal_genes <- c("Ir75d", "Ir76b", "Ir75n1")
broad_genes <- c("Ir25a", "Ir8a")
specific_genes <- c("Ir75f", "Ir76a", "Ir75n3")

status_levels <- c("No expression", "Accessible only", 
                   "Lack expression", "Low expression", "High expression")

# 创建排序数据框
df_order <- data.frame(cluster = cluster_order)

# 添加所有基因列（保持原始顺序）
for (g in gene_order_by_no) {
  df_order[[g]] <- factor(
    expr_status_matrix[as.character(df_order$cluster), g],
    levels = status_levels
  )
}

# 1. 对Universal/Broad基因进行块内排序
# 初始化排序序列
current_order <- 1:nrow(df_order)

# 按顺序处理每个universal/broad基因
for (g in c(universal_genes, broad_genes)) {
  # 仅在当前基因存在时排序
  if (g %in% colnames(df_order)) {
    # 获取当前基因的状态值
    status_values <- as.numeric(df_order[[g]])
    
    # 在当前块内排序（保持前次排序的相对顺序）
    new_order <- order(current_order, status_values)
    
    # 更新排序
    df_order <- df_order[new_order, ]
    current_order <- order(new_order)
  }
}

# 2. 对Specific基因进行反向排序
for (g in specific_genes) {
  # 仅在当前基因存在时排序
  if (g %in% colnames(df_order)) {
    # 创建反向排序指标：High expression优先在底部
    reverse_score <- ifelse(
      df_order[[g]] == "High expression", 0,  # 高表达得最低分
      ifelse(
        df_order[[g]] == "Low expression", 1,
        ifelse(
          df_order[[g]] == "Lack expression", 2,
          ifelse(
            df_order[[g]] == "Accessible only", 3,
            4  # No expression得最高分
          )
        )
      )
    )
    
    # 在当前块内反向排序
    new_order <- order(current_order, reverse_score)
    
    # 更新排序
    df_order <- df_order[new_order, ]
    current_order <- order(new_order)
  }
}

# 获取最终聚类顺序
new_cluster_order <- df_order$cluster

new_cluster_order <- c("0", "5", "8", "43", "3", "47", "37", "16", "31", "25", 
                       "15", "42", "29", "27", "32", "45", "6", "51", "55", "40", 
                       "34", "49", "50", "58", "1", "18", "46", "9", "39", "59", 
                       "54", "57", "53", "20", "12", "41", "19", "13", "35", "30", 
                       "11", "2", "48", "56", "33", "14", "38", "23", "28", "17", 
                       "36", "26", "44", "10", "4", "52", "21", "22", "24", "7")
  
# 更新绘图数据中的因子水平
plot_data$cluster <- factor(plot_data$cluster, levels = new_cluster_order)
plot_data$gene <- factor(plot_data$gene, levels = gene_order_by_no)
  
# 计算分组线位置（基于新基因顺序）
universal_genes <- c("Ir75d", "Ir76b", "Ir75n1") 
broad_genes <- c("Ir25a", "Ir8a")
specific_genes <- c("Ir75f", "Ir76a", "Ir75n3")
  
universal_end <- max(which(levels(plot_data$gene) %in% universal_genes)) + 0.5
broad_end <- max(which(levels(plot_data$gene) %in% broad_genes)) + 0.5
  
# 绘图
point_colors <- c("High expression" = "#4d8076", 
                  "Lack expression" = "#4d8076",
                  "Low expression" = "#70a19e",
                  "Accessible only" = "white")
  
pdf("/data/R04/liyh526/project/01_Honeybee_IR/06_final/Fig1/Fig1_Dotplot/Only_IR_Neuron_recluster_Dotplot_Annotation.pdf", width=10, height=16)
ggplot(plot_data, aes(x = cluster, y = gene)) +
  geom_tile(fill = "white", color = "grey90", linewidth = 0.5) +
    
  # 绘制不同表达状态的点
  geom_point(
    data = subset(plot_data, expr_status == "Low expression"),
    aes(fill = expr_status),
    size = 6, shape = 21, color = "grey60", alpha = 0.5
  ) +
  geom_point(
    data = subset(plot_data, expr_status == "Lack expression"),
    aes(fill = expr_status),
    size = 3, shape = 21, color = "grey30"
  ) +
  geom_point(
    data = subset(plot_data, expr_status == "High expression"),
    aes(fill = expr_status),
    size = 6, shape = 21, color = "grey30"
  ) +
  geom_point(
    data = subset(plot_data, expr_status == "Accessible only"),
    shape = 21, size = 3, color = "grey30", fill = "white"
  ) +
  geom_point(
    data = subset(plot_data, expr_status == "No expression"),
    shape = 4, size = 2.5, color = "grey50"
  ) +
    
  # 设置颜色和主题
  scale_fill_manual(values = point_colors) +
  geom_hline(
    yintercept = c(universal_end, broad_end), 
    color = "grey30", linetype = "dashed", linewidth = 0.6
  ) +
  theme_minimal(base_size = 14) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, face = "bold"),
    axis.text.y = element_text(face = "italic"),
    axis.title = element_blank(),
    panel.grid = element_blank(),
    panel.border = element_rect(fill = NA, color = "black", linewidth = 1),
    legend.position = "none",
    plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
    plot.subtitle = element_text(hjust = 0.5, color = "grey40")
  ) +
  # 坐标轴设置
  scale_x_discrete(expand = c(0, 0)) +
  scale_y_discrete(expand = c(0, 0)) +
  labs(
    title = "Ionotropic Receptor Expression Pattern",
    subtitle = "Expression threshold: avg > 0.1 and > 10% cells"
  ) +
  coord_flip(clip = "off") +
    
  # 添加分组标签（调整到图上方）
  annotate(
    "text", 
    x = length(unique(ORN$Annotation)) + 1.5,  # 放在cluster轴上方
    y = c(mean(which(levels(plot_data$gene) %in% universal_genes)),
          mean(which(levels(plot_data$gene) %in% broad_genes)),
          mean(which(levels(plot_data$gene) %in% specific_genes))),
    label = c("Universal", "Broad", "Specific"),
    hjust = 0.5, fontface = "bold", size = 5,
    color = c("#1f77b4", "#9365b4", "#c65293")
  )
dev.off()