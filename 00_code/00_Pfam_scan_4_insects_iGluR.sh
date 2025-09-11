#!/bin/bash
# 设置工作目录和物种列表
species=(
  "Apis_mellifera" 
  "Drosophila_melanogaster"
  "Aedes_aegypti"
  "Ooceraea_biroi"
)

# 定义iGluR必需的Pfam结构域名称（使用hmm_name列）
required_domains=("ANF_receptor" "Lig_chan")

# 设置Pfam数据库路径
pfam_db="/data/R04/liyh526/project/01_Honeybee_IR/05_homology_find_in_4/input/pfam"  # 替换为实际Pfam数据库路径

# 创建输出目录
mkdir -p iGluR_results

# 循环处理每个物种
for sp in "${species[@]}"; do
  # 输入文件路径
  fasta_file="${sp}.faa"
  
  # 检查输入文件是否存在
  if [[ ! -f "$fasta_file" ]]; then
    echo "Error: Missing input file - $fasta_file"
    continue
  fi
  
  # Step 1: 运行pfam_scan.py
  pfam_out="iGluR_results/${sp}_pfam_scan.csv"
  
  echo "Processing $sp..."
  /data/R04/liyh526/project/01_Honeybee_IR/05_homology_find_in_4/input/pfam_scan/pfam_scan.py "$fasta_file" "$pfam_db" -out "$pfam_out" -cpu 16 -evalue 1e-5
  
  # Step 2: 筛选具有完整iGluR结构域的蛋白
  # 使用awk高效处理CSV文件，匹配hmm_name列
  awk -F, -v domains="${required_domains[*]}" '
  BEGIN {
    split(domains, req_dom, " ")
    for (i in req_dom) required[req_dom[i]] = 1
  }
  NR > 1 {  # 跳过标题行
    seq_id = $1
    hmm_name = $7
    
    # 去除可能的引号
    gsub(/"/, "", seq_id)
    gsub(/"/, "", hmm_name)
    
    # 记录每个序列的所有结构域名称
    domains_all[seq_id][hmm_name] = 1
  }
  END {
    # 检查每个序列是否包含所有必需结构域
    for (seq in domains_all) {
      missing = 0
      for (dom in required) {
        if (!(dom in domains_all[seq])) {
          missing = 1
          break
        }
      }
      if (!missing) print seq
    }
  }' "$pfam_out" > "iGluR_results/${sp}_iGluR_ids.txt"
  
  # Step 3: 提取iGluR蛋白序列
  if [[ -s "iGluR_results/${sp}_iGluR_ids.txt" ]]; then
    count=$(wc -l < "iGluR_results/${sp}_iGluR_ids.txt")
    echo "Found $count iGluR proteins for $sp"
    
    # 使用seqkit提取序列
    seqkit grep -f "iGluR_results/${sp}_iGluR_ids.txt" \
      "$fasta_file" > "iGluR_results/${sp}_iGluR.faa"
  else
    echo "No complete iGluR proteins found for $sp"
    # 创建空结果文件
    touch "iGluR_results/${sp}_iGluR.faa"
  fi
done

echo "Processing complete. Results saved in iGluR_results/"