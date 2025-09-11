library(chorddiag)
library(Matrix)

ORN <- readRDS("/data/R04/liyh526/project/01_Honeybee_IR/00_data/00_rds_data/Validate_OSN_merge_Orco_neg_neuron_tsne.rds")

gene_order_final <- c("Or2", # Orco
                      #"LOC552552",
                      "LOC726019", "LOC551704",
                      #"LOC408645",
                      "LOC102653640",
                      #"LOC412818", # Universal
                      "LOC412949", "LOC409777",
                      #"LOC412993", # Broad
                      "LOC725724","LOC102656088","LOC100576097" # Specific
)

highReceptorRNA = ORN@assays$raw_RNA@counts[gene_order_final,]
highReceptorRNA.data = ORN@assays$raw_RNA@data[gene_order_final,]

highReceptorRNA.expCol = highReceptorRNA[,colSums(highReceptorRNA) != 0]
highReceptorRNA.data.expCol = highReceptorRNA.data[,colSums(highReceptorRNA.data) != 0]

coexp.df <- data.frame(matrix(ncol = 9, nrow = 9))
colnames(coexp.df) <- rownames(highReceptorRNA.expCol)
rownames(coexp.df) <- rownames(highReceptorRNA.expCol)

coexp.df2 <- data.frame(matrix(ncol = 9, nrow = 9))
colnames(coexp.df2) <- rownames(highReceptorRNA.expCol)
rownames(coexp.df2) <- rownames(highReceptorRNA.expCol)

coexp.df3 <- data.frame(matrix(ncol = 9, nrow = 9))
colnames(coexp.df3) <- rownames(highReceptorRNA.expCol)
rownames(coexp.df3) <- rownames(highReceptorRNA.expCol)

# 在计算完协同表达矩阵后添加名称映射
loc_to_symbol <- c(
  "Or2" = "Orco",
  "LOC552552" = "Nmdar2-Like",
  "LOC726019" = "Ir75d",
  "LOC551704" = "Ir76b",
  "LOC102653640" = "Ir75n1",
  "LOC408645" = "KaiR1D",
  "LOC412818" = "Nmdar2",
  "LOC412949" = "Ir25a",
  "LOC412993" = "Ekar",
  "LOC409777" = "Ir8a",
  "LOC725724" = "Ir76a",
  "LOC102656088" = "Ir75f",
  "LOC100576097" = "Ir75n3"
)

# 对三个矩阵应用相同的名称转换
rename_matrix <- function(mat, mapping) {
  new_names <- sapply(rownames(mat), function(x) {
    ifelse(x %in% names(mapping), mapping[x], x)
  })
  rownames(mat) <- new_names
  colnames(mat) <- new_names
  return(mat)
}

coexp.df <- rename_matrix(coexp.df, loc_to_symbol)
coexp.df2 <- rename_matrix(coexp.df2, loc_to_symbol)
coexp.df3 <- rename_matrix(coexp.df3, loc_to_symbol)

DefaultAssay(ORN) <- "RNA"
p <- DotPlot(ORN, 
        features = gene_order_final,
        group.by = "Annotation",
        cols = c("lightgrey", "blue"),
        dot.scale = 6) +
  xlab('') + ylab('Cluster') +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 9),
    legend.position = "right",
    plot.margin = margin(t = 40, b = 5, l = 5, r = 20)
  )

dotplot_data <- p$data

# normalized expression
for (i in 1:9) {
  for (j in 1:9) {
    if (i==j) { coexp.df[i,j] <- 0 }
    else {
      lapply(1:length(ncol(highReceptorRNA.data.expCol)), function(x){
        print(c('i = ', i))
        print(c('j = ', j))
        iExpressedCells.mx = highReceptorRNA.data.expCol[,highReceptorRNA.data.expCol[i,] >= 1]
        jExpressedN=NCOL(iExpressedCells.mx[, iExpressedCells.mx[j,] >= 1])
        # print(jExpressedN)
        gene1<- gene_order_final[i]
        gene2<- gene_order_final[j]
        cluster_id<-dotplot_data[which(dotplot_data$features.plot%in%c(gene1,gene2)),]$id;
        if(length(which(duplicated(cluster_id)))){
          coexp.df[i,j]<<-jExpressedN
          iExpressedCells.mx2 = highReceptorRNA.data.expCol[,highReceptorRNA.data.expCol[i,] >= 2]
          jExpressedN2=NCOL(iExpressedCells.mx2[, iExpressedCells.mx2[j,] >= 2])
          # print(jExpressedN)
          coexp.df2[i,j]<<-jExpressedN2
          
          iExpressedCells.mx3 = highReceptorRNA.data.expCol[,highReceptorRNA.data.expCol[i,] >= 1.5]
          jExpressedN3=NCOL(iExpressedCells.mx3[, iExpressedCells.mx3[j,] >= 1.5])
          # print(jExpressedN)
          coexp.df3[i,j]<<-jExpressedN3;}
        else{coexp.df[i,j]<<-0;
        coexp.df2[i,j]<<-0;
        coexp.df3[i,j]<<-0;
        }
      })
    }
  }
}

mytsnecolors <- c(
  "#005d56", "#476D87", "#BD956A", "#968175", "#DDA03C", "#DCC1DD",
  "#FFBCBC", "#F3B1A0", "#D9D9D9", "#9FA3A8", "#F1BB72", "#FDB462", 
  "#23452F", "#3A6963", "#FF9A76", "#E63863", "#FFCC1D", "#FF9999", 
  "#E41A1C", "#E95C59", "#9C27B0", "#AB3282", "#8C549C", "#5F3D69", 
  "#BC80BD", "#BEBADA", "#FFB897", "#E4D457", "#FB8072", "#E39A35",  
  "#FCCDE5", "#B53E2B", "#712820", "#678983", "#A19882", "#FF7F00",  
  "#585658", "#E0D4CA", "#B3DE69", "#D6E7A3", "#4DAF4A", "#53A85F", 
  "#68A180", "#116530", "#24A19C", "#8DD3C7", "#F781BF", "#E59CC4", 
  "#CCEBC5", "#57C3F3", "#C1E6F3", "#CCE0F5", "#58A4C3", "#80B1D3", 
  "#6778AE", "#625D9E", "#377EB8", "#344CB7", "#161853", "#AA9A59"
)
hex_codes1<- mytsnecolors

# 创建弦图对象
chord_plot <- chorddiag(
  as.matrix(coexp.df),
  groupColors = mytsnecolors,
  groupnamePadding = 12,
  groupnameFontsize = 12,
  showTicks = FALSE
)

htmlwidgets::saveWidget(
  chord_plot,
  file = "/data/R04/liyh526/project/01_Honeybee_IR/06_final/Fig1/Fig1-coexp_IR_chord_plot.html",
  selfcontained = TRUE
)

#install.packages("webshot")
webshot::install_phantomjs()

# 在保存HTML后添加以下代码
library(webshot)
webshot(
  url = "/data/R04/liyh526/project/01_Honeybee_IR/06_final/Fig1/Fig1-coexp_IR_chord_plot.html",
  file = "/data/R04/liyh526/project/01_Honeybee_IR/06_final/Fig1/Fig1-coexp_IR_chord_plot.pdf",
  delay = 10,       # 等待10秒确保页面加载完成
  vwidth = 1200,    # 设置宽度
  vheight = 1200    # 设置高度
)
