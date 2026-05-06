library(data.table)
library(rtracklayer)
library(GenomicRanges)

setwd("/home/rstudio/2. Validation test/HN_DNA_CytoSNP_202603/20260223/")

dir.create(path = "ASCAT_input",showWarnings = F)
cyto_SNP_files <- list.files(pattern = "FinalReport.\\.txt$", full.names = TRUE)
print(cyto_SNP_files)

target_files <- cyto_SNP_files 

# 2. 定義需要的欄位
target_columns <- c("Sample ID", "Position", "Chr", "SNP Name", "Log R Ratio", "B Allele Freq")

# 3. 批次讀取並過濾
data_list <- lapply(target_files, function(x) {
  df <- fread(x, skip = "Sample ID", select = target_columns)
  return(df)
})

final_df <- rbindlist(data_list, use.names = TRUE, fill = TRUE)
final_df <- final_df[!(Chr %in% c("0", "XY"))]
# final_df[(final_df$`Sample ID` != c("H0148N","H0201N")),]
# filtered_df <- final_df[!( `Sample ID` %in% c("H0148N", "H0201N") )]

# 4. 設定名稱 (確保與讀取的檔案清單對應)
names(data_list) <- basename(target_files)

# Long to Wide (log R)
short_df <- dcast(final_df, `SNP Name` + Chr + Position ~ `Sample ID`, 
                  value.var = "Log R Ratio")
# Long to Wide (BAF)
short_df_baf <- dcast(final_df, `SNP Name` + Chr + Position ~ `Sample ID`,
                      value.var = "B Allele Freq")

# for GC and rep
cols <- c("SNP Name", "Chr", "Position")
selected_df <- short_df[, cols, with = FALSE]
# fwrite(selected_df, "./ASCAT_input/20260428_SNP_pos_hg38.txt", sep = "\t", quote = FALSE)

## for repication hg38 to hg19
# 1. 讀取連鎖檔 (Chain file)
# 請修改路徑指向您下載的 chain 檔案
chain <- import.chain("/home/rstudio/hg382hg19_chain/hg38ToHg19.over.chain")

# 2. 將 data.table 轉換為 GRanges 物件
# 確保 Chr 有 "chr" 前綴，否則與 chain 檔對不起來
if (!grepl("chr", selected_df$Chr[1])) {
  selected_df$Chr <- paste0("chr", selected_df$Chr)
}

gr_hg38 <- makeGRangesFromDataFrame(selected_df,
                                    keep.extra.columns = TRUE,
                                    ignore.strand = TRUE,
                                    seqnames.field = "Chr",
                                    start.field = "Position",
                                    end.field = "Position")

# 3. 執行 LiftOver
cur_hg19_list <- liftOver(gr_hg38, chain)
gr_hg19 <- unlist(cur_hg19_list)

# 4. 轉回 data.table 並整理格式
hg19_df <- as.data.table(gr_hg19)

# 重新命名與挑選欄位，以符合您原本的格式
final_hg19_pos <- hg19_df[, .(`SNP Name` = SNP.Name, 
                              Chr = seqnames, 
                              Position = start)]
# 1. 以 SNP Name 為基準進行交集 (Inner Join)
# 這會確保兩邊都存在的 SNP 才保留，並同時擁有 hg38 與 hg19 的座標
combined_pos <- merge(selected_df, final_hg19_pos, 
                      by = "SNP Name", 
                      suffixes = c("_hg38", "_hg19"))

# 2. 檢查交集後的數量變化
cat("原始 hg38 數量:", nrow(selected_df), "\n")
cat("轉換成功 hg19 數量:", nrow(final_hg19_pos), "\n")
cat("交集後共同數量:", nrow(combined_pos), "\n")

# 定義正確的染色體順序 (1~22, X, Y, MT 視需求保留)
chr_order <- paste0("chr", c(1:22, "X", "Y"))

# 將 Chr 欄位轉換為 Factor，這樣 R 就會知道要按照我們給的順序排，而不是字母順序
combined_pos[, Chr_hg38 := factor(Chr_hg38, levels = chr_order)]
combined_pos[, Chr_hg19 := factor(Chr_hg19, levels = chr_order)]

# 執行排序：先排染色體 (Chr_hg38)，再排物理位置 (Position_hg38)
# 這裡統一以 hg38 的座標為基準進行排序
setorder(combined_pos, Chr_hg38, Position_hg38)

# (選用) 排序後，如果有未能對應到 chr_order 的奇異染色體會變成 NA，可以過濾掉
combined_pos <- combined_pos[!is.na(Chr_hg38) & !is.na(Chr_hg19)]

# 3. 分別存出「乾淨且一致」的座標檔
# 存出給 hg38 流程用的 for createGCcontentFile.R
fwrite(combined_pos[, .(`SNP Name`, Chr = Chr_hg38, Position = Position_hg38)], 
       "/home/rstudio/2. Validation test/HN_DNA_CytoSNP_202603/20260223/ASCAT_input/20260428_SNP_pos_hg38_final.txt", 
       sep = "\t", quote = FALSE)

# 存出給 hg19 流程用的 for replication without "chr"
fwrite(combined_pos[, .(`SNP Name`, Chr = gsub("chr", "", Chr_hg19), Position = Position_hg19)], 
       "/home/rstudio/2. Validation test/HN_DNA_CytoSNP_202603/20260223/ASCAT_input/20260428_SNP_pos_hg19_final.txt",
       sep = "\t", quote = FALSE)

## Go to run
# Rscript createReplicTimingFile.R '/home/rstudio/2. Validation test/HN_DNA_CytoSNP_202603/20260223/ASCAT_input/20260428_SNP_pos_hg19.txt'
# Rscript createGCcontentFile.R "/home/rstudio/2. Validation test/HN_DNA_CytoSNP_202603/20260223/ASCAT_input/20260428_SNP_pos_hg38_final.txt" 24 "/home/rstudio/2. Validation test/bwa_ref/GRCh38.d1.vd1.fa"

## rep hg19 to hg38 name and pos
rt_data <- fread("/home/rstudio/ReplicationTiming_SNPloci.txt")
colnames(rt_data)[colnames(rt_data) == "V1"] <- "SNP Name"

# 2. 與我們先前建立的 master list (combined_pos) 做交集
# combined_pos 裡已經有 hg38 和 hg19 的座標
final_master <- merge(combined_pos, rt_data, by = "SNP Name")

# 3. 再次執行生物學排序 (chr1 -> chrX，然後依 Position_hg38 排序)
chr_order <- paste0("chr", c(1:22, "X", "Y"))
final_master[, Chr_hg38 := factor(Chr_hg38, levels = chr_order)]
setorder(final_master, Chr_hg38, Position_hg38)

# 4. 根據這份「最終交集名單」過濾並排序 LRR 與 BAF 矩陣
valid_snps <- final_master$`SNP Name`

# 過濾 LRR 矩陣
short_df_final <- short_df[`SNP Name` %in% valid_snps]
short_df_final <- short_df_final[match(valid_snps, `SNP Name`)]

# 過濾 BAF 矩陣
short_df_baf_final <- short_df_baf[`SNP Name` %in% valid_snps]
short_df_baf_final <- short_df_baf_final[match(valid_snps, `SNP Name`)]

# 5. 輸出所有 ASCAT 最終輸入檔 (全部統一為 hg38 座標)

cell_line_cols <- names(final_master)[8:ncol(final_master)]
# B. 最終 RT 檔 (對齊 hg38 順序)
rt_output <- final_master[, .SD, .SDcols = c("SNP Name","Chr_hg38","Position_hg38", cell_line_cols)]
colnames(rt_output)[colnames(rt_output) == "Chr_hg38"] <- "Chr"
colnames(rt_output)[colnames(rt_output) == "Position_hg38"] <- "Position"
rt_output$Chr = gsub("chr", "", rt_output$Chr)
row.names(rt_output) = rt_output$`SNP Name`
rt_output$`SNP Name` <- NULL
write.table(
  rt_output,
  "./ASCAT_input/20260428_RT_hg38_final.txt",
  sep = "\t",
  quote = FALSE,
  col.names = NA
)


test<-fread("./ASCAT_input/20260428_RT_hg38_final.txt")
# LRR
fwrite(short_df_final, "./ASCAT_input/20260428_LRR_hg38_final.txt", sep = "\t", quote = FALSE)

# BAF
fwrite(short_df_baf_final, "./ASCAT_input/20260428_BAF_hg38_final.txt", sep = "\t", quote = FALSE)

# GC remove chr
GC_content <- fread("/home/rstudio/GCcontent_SNPloci.txt")
GC_content$Chr = gsub("chr", "", GC_content$Chr)
row.names(GC_content) = GC_content$V1
GC_content$V1 <- NULL
write.table(
  GC_content,
  "./ASCAT_input/20260428_GC_hg38_final.txt",
  sep = "\t",
  quote = FALSE,
  col.names = NA
)


test<-fread("./ASCAT_input/20260428_GC_hg38_final.txt")

cat("最終 SNP 總數：", nrow(final_master), "\n")
