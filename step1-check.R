
a = read.table('GSE111229_Mammary_Tumor_fibroblasts_768samples_rawCounts.txt.gz', header = T, sep = '\t')

a[1:6, 1:4]

# 筛选表达量合格的行
dat = a[apply(a, 1, function(x) {sum(x>1) > floor(ncol(a)/50)}), ]

dat[1:4, 1:4]

# sum(dat[,3])
# log2( 18*1000000/sum(dat[,3]) + 1)
# 去除文库大小差异, CPM归一化
dat = log2(edgeR::cpm(dat)+1)
dat[1:4, 1:4]


x = 1:10
y = 2*x
z = rnorm(10)
tmp = data.frame(x, y, z)

# 层次聚类
hc = hclust(dist(t(dat)))

clus = cutree(hc, 4)
group_list = as.factor(clus)
table(group_list)

colnames(dat)
library(stringr)
plate = str_split(colnames(dat), '_', simplify = T)[,3]
table(plate)

# 统计每个样本有表达的多少行
n_g = apply(a, 2, function(x) sum(x>1))

df = data.frame(g=group_list, plate=plate, n_g=n_g)

df$all = 'all'

save(a, dat, df, file = 'input.Rdata')

# 再对rpkm的矩阵操作，因为我们需要复现文章的图片
if(F) {
    a = read.table('GSE111229_Mammary_Tumor_fibroblasts_768samples_rpkmNormalized.txt.gz', header = T, sep = '\t')
    dat = a[apply(a, 1, function(x) sum(x>0) > floor(ncol(a)/50)), ]
    hc = hclust(dist(t(dat))) ## 样本间层次聚类
    plot(hc)
    
    clus = cutree(hc, 4)
    group_list = as.factor(clus)
    table(group_list)
    
    # 提取批次信息
    colnames(dat)
    library(stringr)
    plate = str_split(colnames(dat), '_', simplify = T)[,3]
    table(plate)
    
    n_g = apply(a, 2, function(x) sum(x>0))
    df = data.frame(g = group_list, plate = plate, n_g = n_g)    
    df$all = 'all'
    
    save(a, dat, df, file = 'input_rpkm.Rdata')
        
            
}









