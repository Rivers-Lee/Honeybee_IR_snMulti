#install.packages("RIdeogram")
library(RIdeogram)

AmelChromosome <- read.table(
  "/data/R04/liyh526/project/01_Honeybee_IR/01_genome_info/Honeybee/Honeybee_chrom_length.txt",
  header = TRUE,
  sep = "\t",
  colClasses = c("character", "integer", "integer")  # 明确指定列类型
)
AmelChromosome
sapply(AmelChromosome, class)

gtf_data <- read.table(
  "/data/R04/liyh526/project/01_Honeybee_IR/01_genome_info/Honeybee/GCF_003254395.2_Amel_HAv3.1_genomic_UPDATED_chemoreceptors_CORRECTED_geneid_withoutGroup.gtf",
  sep = "\t", header = FALSE, comment.char = "#"
)
unique_gtf_chr <- unique(gtf_data$V1)
print(unique_gtf_chr)

# 检查每条染色体的最大坐标
for(chr in unique(gtf_data$V1)) {
  max_pos <- max(gtf_data[gtf_data$V1 == chr, "V4"], na.rm = TRUE)
  chr_length <- AmelChromosome[AmelChromosome$Chr == chr, "End"]
  cat(sprintf("Chromosome %s: GTF max=%d, Karyotype length=%d\n", 
              chr, max_pos, chr_length))
}

# 计算每条染色体上的基因密度
library(dplyr)

# 提取基因位置（GTF第3列为"transcript"）
genes <- gtf_data[gtf_data$V3 == "transcript", c("V1", "V4")]
colnames(genes) <- c("Chr", "Position")

# 创建滑动窗口
window_size <- 200000
density_data <- data.frame()

for (chr in unique(AmelChromosome$Chr)) {
  chr_length <- AmelChromosome[AmelChromosome$Chr == chr, "End"]
  starts <- seq(0, chr_length, by = window_size)
  ends <- pmin(starts + window_size, chr_length)
  
  # 计算每个窗口的基因数量
  gene_count <- sapply(seq_along(starts), function(i) {
    sum(genes$Chr == chr & 
          genes$Position >= starts[i] & 
          genes$Position <= ends[i])
  })
  
  # 添加到结果
  density_data <- rbind(density_data, data.frame(
    Chr = chr,
    Start = starts,
    End = ends,
    Value = gene_count
  ))
}

density_data

gene_list <- read.table(
  "/data/R04/liyh526/project/01_Honeybee_IR/01_genome_info/Honeybee/Honeybee_IR_gene_position.txt",
  header = T,
  sep = "\t")

setwd("/data/R04/liyh526/project/01_Honeybee_IR/06_final/FigS1/Chromesome_Ideogram")
ideogram(karyotype=AmelChromosome, #染色体信息
         overlaid=density_data, #基因密度
         label=gene_list, #目标基因标签
         label_type="marker",width=135)
convertSVG("/data/R04/liyh526/project/01_Honeybee_IR/06_final/FigS1/Chromesome_Ideogram/chromosome.svg",device="png")
