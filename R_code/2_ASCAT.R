library(ASCAT)
library(data.table)

setwd("/home/rstudio/2. Validation test/HN_DNA_CytoSNP_202603/20260223/ASCAT_input/")

ascat.bc <- ascat.loadData(
  Tumor_LogR_file = "20260428_LRR_hg38_final.txt",
  Tumor_BAF_file  = "20260428_BAF_hg38_final.txt",
  gender          = rep("XY",8),
  genomeVersion   = "hg38"
)
ascat.plotRawData(ascat.bc, img.prefix = "Before_correction_")

ascat.bc = ascat.correctLogR(ascat.bc,
                             GCcontentfile = "20260428_GC_hg38_final.txt",
                             replictimingfile = "20260428_RT_hg38_final.txt")
ascat.plotRawData(ascat.bc, img.prefix = "After_correction_")
# 使用平台 preset / germline genotype prediction
gg = ascat.predictGermlineGenotypes(
  ascat.bc,
  platform = "IlluminaCytoSNP850k"
)

ascat.bc = ascat.aspcf(ascat.bc, ascat.gg=gg)

ascat.plotSegmentedData(ascat.bc)

ascat.output = ascat.runAscat(ascat.bc, write_segments = TRUE)

QC = ascat.metrics(ascat.bc,ascat.output)
save(ascat.bc, ascat.output, QC, file = 'ASCAT_objects.Rdata')


# Read fiel and save ploidy and purity
seg.files <- list.files(
  pattern = "^.*segments\\.txt$",
  full.names = TRUE
)

purity <- ascat.output[["purity"]]
ploidy <- ascat.output[["ploidy"]]

for (file in seg.files) {
  
  seg <- fread(file)
  
  # 從檔名抓 sample name
  # 例如 ./H0050I.segments.txt -> H0050I
  sample_id <- sub("\\.segments\\.txt$", "", basename(file))
  
  seg[, purity := purity[sample_id]]
  seg[, ploidy := ploidy[sample_id]]
  
  # 覆蓋原檔
  fwrite(seg, file, sep = "\t")
}
