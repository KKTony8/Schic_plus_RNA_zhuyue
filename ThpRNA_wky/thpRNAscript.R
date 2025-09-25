getwd()
setwd("C:/Users/tommyw/Desktop/zhuyuethpRNA/thpRNAQC")
.libPaths()
library(Seurat)
library(dplyr)
library(ggplot2)
# 加载 Seurat 包
library(Seurat)
library(dplyr)
library(ggplot2)

### 1.数据预处理部分
# 设置文件路径
file_path <- "你的桌面路径/zhuyuethpRNA/thpsthpRNAcount/333.txt"
# 读取数据
data <- read.table("你的桌面路径/zhuyuethpRNA/thpsthpRNAcount/333.txt",
                   header = TRUE, sep = "\t", check.names = FALSE)

# 查看列名
colnames(data)

# 创建一个映射表，将旧列名映射为新列名
old_names <- c(
  "JZ24177048-RNA-1-RNA-1_combined_R1",
  "JZ24177051-RNA-4-RNA-4_combined_R1",
  "JZ24182315-RNA-5-RNA-5_combined_R1",
  "JZ24182316-RNA-6-RNA-6_combined_R1",
  "R54_1",
  "R56_1"
)

new_names <- c(
  "thp_RNA_1",
  "thp_RNA_4",
  "thp_RNA_5",
  "thp_RNA_6",
  "thp_RNA_54",
  "thp_RNA_56"
)

# 替换列名
colnames(data)[colnames(data) %in% old_names] <- new_names

# 检查修改后的列名
colnames(data)

# 保存修改后的数据
write.table(data, "你的桌面路径/zhuyuethpRNA/thpsthpRNAcount/333_renamed.txt",
            sep = "\t", quote = FALSE, row.names = FALSE)

# 读取数据
data <- read.table("你的桌面路径/zhuyuethpRNA/thpsthpRNAcount/333_renamed.txt",
                   header = TRUE, sep = "\t", check.names = FALSE)

# 查看原始列名
colnames(data)

# 定义函数处理列名
clean_names <- function(names_vec) {
  names_vec <- gsub("^SThp1_RNA_", "sthp_RNA_", names_vec)  # SThp1 -> sthp
  names_vec <- gsub("^Thp1_RNA_", "thp_RNA_", names_vec)    # Thp1 -> thp
  names_vec <- gsub("^sthp_RNA_", "sthp_RNA_", names_vec)   # 小写 sthp 保留
  names_vec <- gsub("^thp-1RNA-", "thp_RNA_", names_vec)    # thp-1RNA- -> thp_RNA_
  names_vec <- gsub("_1_val$", "", names_vec)               # 去掉尾部 _1_val
  names_vec <- gsub("_SThp1_RNA_.*$", "", names_vec)        # 处理重复复杂列名
  return(names_vec)
}

# 应用函数
colnames(data) <- clean_names(colnames(data))

# 检查修改后的列名
colnames(data)

# 保存修改后的数据
write.table(data, "你的桌面路径/zhuyuethpRNA/thpsthpRNAcount/333_cleaned.txt",
            sep = "\t", quote = FALSE, row.names = FALSE)

ncol(data)







# 一.数据处理部分#### 


# 一. 20250903重新分析####        
library(Seurat)
library(dplyr)
library(ggplot2)

# 1. 读取合并矩阵
matrix_file <- "你的桌面路径/zhuyuethpRNA/thpsthpRNAcount/merged_final_matrix.txt"
data <- read.table(matrix_file, header = TRUE, row.names = 1, sep = "\t", check.names = FALSE)

# 2. 创建 Seurat 对象
seurat_obj <- CreateSeuratObject(counts = data, project = "thpRNA")

# 3. 添加质控指标
seurat_obj <- PercentageFeatureSet(seurat_obj, pattern = "^MT-", col.name = "percent.mt")  # 可选线粒体比例
seurat_obj
# 4. 根据列名添加分组信息
cell_names <- colnames(seurat_obj)
group <- ifelse(grepl("^thp", cell_names, ignore.case = TRUE), "THP", "STHP")
seurat_obj$Group <- group
seurat_obj$Group
group
table(seurat_obj$Group)
# 5. 绘制按组的小提琴图
VlnPlot(
  seurat_obj,
  features = c("nFeature_RNA", "nCount_RNA"),
  group.by = "Group",   # 按 Group 分组
  pt.size = 0.1,
  ncol = 3
) + 
  scale_fill_manual(values = c("THP" = "#1f78b4", "STHP" = "#33a02c"))  # 自定义颜色


# 二. 20250903重新分析，阈值去掉样本####  


filtered_seurat <- subset(seurat_obj, subset = nFeature_RNA > 5000 & nCount_RNA > 1e6 & nCount_RNA < 1.5e7)

# 查看过滤后的细胞数
filtered_cells_count <- ncol(filtered_seurat)
cat("过滤后剩余细胞数:", filtered_cells_count, "\n")

# 查看按 Group 分组的数量
table(filtered_seurat$Group)

VlnPlot(
  filtered_seurat,
  features = c("nFeature_RNA", "nCount_RNA"),
  group.by = "Group",
  pt.size = 0.1,
  ncol = 3
) + 
  scale_fill_manual(values = c("THP" = "#1f78b4", "STHP" = "#33a02c"))


# 三. 20250903重新分析，看看要不要去批次#### 


# 1. 标准化数据
filtered_seurat <- NormalizeData(filtered_seurat)

# 2. 找高变基因
filtered_seurat <- FindVariableFeatures(filtered_seurat, selection.method = "vst", nfeatures = 2000)

# 3. 缩放数据
all_genes <- rownames(filtered_seurat)
filtered_seurat <- ScaleData(filtered_seurat, features = all_genes)

# 4. PCA
filtered_seurat <- RunPCA(filtered_seurat, features = VariableFeatures(object = filtered_seurat))

# 5. 查看前几个 PC 的重要性
print(filtered_seurat[["pca"]], dims = 1:5, nfeatures = 5)

# 6. PCA 可视化
DimPlot(filtered_seurat, reduction = "pca", group.by = "Group", pt.size = 1) +
  scale_color_manual(values = c("THP" = "#1f78b4", "STHP" = "#33a02c")) +
  ggtitle("PCA of filtered cells")

# 5. 构建邻接图
filtered_seurat <- FindNeighbors(filtered_seurat, dims = 1:20)

# 6. UMAP
filtered_seurat <- RunUMAP(filtered_seurat, dims = 1:20)

# 7. 按 Group 可视化
DimPlot(filtered_seurat, reduction = "umap", group.by = "Group", pt.size = 1) +
  scale_color_manual(values = c("THP" = "#1f78b4", "STHP" = "#33a02c")) +
  ggtitle("UMAP by Group")


# 三. 20250903重新分析，不用去批次直接分析####  


# 1. 运行 PCA 后，可以直接用 ElbowPlot
# 将 ident 中的 "SThp1" 改成 "sthp"
levels(filtered_seurat)  # 查看当前 level
filtered_seurat <- RenameIdents(filtered_seurat, "SThp1" = "sthp")

# 再次检查
table(Idents(filtered_seurat))

ElbowPlot(filtered_seurat, ndims = 50) + 
  ggtitle("Elbow Plot to determine number of PCs")

# 5. 运行 UMAP
# 通常我们使用前10个主成分（PCs），可以根据需要调整
filtered_seurat <- RunUMAP(filtered_seurat, dims = 1:10)

# 6. 可视化 UMAP 图
DimPlot(filtered_seurat, reduction = "umap", label = TRUE, pt.size = 0.5) +
  ggtitle("UMAP of filtered Seurat object")






# 四. 20250903重新分析，细胞注释#### 

seurat.markers <- FindAllMarkers(filtered_seurat, only.pos = TRUE,
                                 min.pct = 0.25, logfc.threshold = 0.25)

# 查看每个簇前几个 marker
top_markers <- seurat.markers %>%
  group_by(cluster) %>%
  slice_max(n = 2, order_by = avg_log2FC)
top_markers




# 3. 在 UMAP 上显示细胞类型注释
DimPlot(filtered_seurat, reduction = "umap", label = TRUE, pt.size = 0.5) +
  ggtitle("UMAP of filtered Seurat object with cell type annotation")

# --------------------------
# 4. 可视化 marker 基因在 UMAP 上的表达（FeaturePlot）
FeaturePlot(filtered_seurat, features = c("LYZ", "PLB1", "CST7", "SLC7A11"))


# 五. 20250903重新分析，差异基因热图#### 


top10 <- seurat.markers %>%
  group_by(cluster) %>%
  slice_max(n = 10, order_by = avg_log2FC) %>%
  pull(gene)
top10
# --------------------------
# 2. 绘制热图
DoHeatmap(filtered_seurat, features = top10, size = 3) +
  ggtitle("Heatmap of top marker genes by cluster")


#六. 20250903重新分析，富集通路#### 


# 安装必要包（如果还没安装）
# BiocManager::install("clusterProfiler")
# BiocManager::install("org.Hs.eg.db")

library(clusterProfiler)
library(org.Hs.eg.db)
library(dplyr)
library(enrichplot)

# --------------------------
# 1. 准备分组的差异基因
# --------------------------
# 假设 seurat.markers 已经有 cluster = "thp" / "sthp"
# 取 avg_log2FC > 0.25
markers_filtered <- seurat.markers %>%
  filter(avg_log2FC > 0.25)

# 分组提取基因
thp_genes <- markers_filtered %>%
  filter(cluster == "thp") %>%
  pull(gene)

sthp_genes <- markers_filtered %>%
  filter(cluster == "sthp") %>%
  pull(gene)

# --------------------------
# 2. 转换为 Entrez ID
# --------------------------
thp_entrez <- bitr(thp_genes, fromType="SYMBOL", toType="ENTREZID", OrgDb=org.Hs.eg.db)$ENTREZID
sthp_entrez <- bitr(sthp_genes, fromType="SYMBOL", toType="ENTREZID", OrgDb=org.Hs.eg.db)$ENTREZID

# --------------------------
# 3. GO 富集分析
# --------------------------
go_thp <- enrichGO(gene=thp_entrez, OrgDb=org.Hs.eg.db, keyType="ENTREZID",
                   ont="BP", pAdjustMethod="BH", pvalueCutoff=0.05, qvalueCutoff=0.05)
go_sthp <- enrichGO(gene=sthp_entrez, OrgDb=org.Hs.eg.db, keyType="ENTREZID",
                    ont="BP", pAdjustMethod="BH", pvalueCutoff=0.05, qvalueCutoff=0.05)

# --------------------------
# 4. KEGG 富集分析
# --------------------------
kegg_thp <- enrichKEGG(gene=thp_entrez, organism="hsa", pvalueCutoff=0.05)
kegg_sthp <- enrichKEGG(gene=sthp_entrez, organism="hsa", pvalueCutoff=0.05)

# --------------------------
# 5.提取前20个通路的名字
# --------------------------

go_thp_names   <- head(as.data.frame(go_thp)$Description, 20)
go_sthp_names  <- head(as.data.frame(go_sthp)$Description, 20)
kegg_thp_names <- head(as.data.frame(kegg_thp)$Description, 20)
kegg_sthp_names<- head(as.data.frame(kegg_sthp)$Description, 20)

# 打印结果
cat("\n=== GO (BP) - thp 前20通路 ===\n")
print(go_thp_names)

cat("\n=== GO (BP) - sthp 前20通路 ===\n")
print(go_sthp_names)

cat("\n=== KEGG - thp 前20通路 ===\n")
print(kegg_thp_names)

cat("\n=== KEGG - sthp 前20通路 ===\n")
print(kegg_sthp_names)


# --------------------------
# 6. 可视化示例
# --------------------------
dotplot(go_thp, showCategory=20) + ggtitle("GO BP enrichment - thp")
dotplot(kegg_thp, showCategory=20) + ggtitle("KEGG enrichment - thp")

dotplot(go_sthp, showCategory=20) + ggtitle("GO BP enrichment - sthp")
dotplot(kegg_sthp, showCategory=20) + ggtitle("KEGG enrichment - sthp")



### 七.看差异基因染色体上位置


# 安装和加载 biomaRt
if(!require(biomaRt)) install.packages("biomaRt")
library(biomaRt)

# 选择 Ensembl GRCh37 (hg19)
ensembl <- useEnsembl(biomart = "genes", 
                      dataset = "hsapiens_gene_ensembl", 
                      GRCh = 37)

# 你的基因列表
genes <- c("LYZ","PLB1","LINC01220","SERPINB2","EMP2","TMEM255A","MCUB",
           "EHHADH","CCR2","KCTD12","CST7","SLC7A11","TTC39B","ALDH1L2",
           "ASNS","ASS1","VLDLR-AS1","SLC6A9","TUBE1","C5AR1")

# 查询 hg19 基因位置信息
gene_locations <- getBM(attributes = c("hgnc_symbol", 
                                       "chromosome_name", 
                                       "start_position", 
                                       "end_position", 
                                       "strand"),
                        filters = "hgnc_symbol",
                        values = genes,
                        mart = ensembl)

# 打印结果
print(gene_locations)

# 保存成文件
write.table(gene_locations, "gene_locations_hg19.txt", 
            sep="\t", quote=FALSE, row.names=FALSE)








