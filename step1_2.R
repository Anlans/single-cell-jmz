rm(list = ls())
options(stringsAsFactors = F)

load(file = 'input.Rdata')
a[1:4, 1:4]
head(df)

group_list = df$g
table(group_list)

# 理解热图、PCA图
cg = names(tail(sort(apply(dat, 1, sd)), 1000)) # 取表达量标准差最大的1000行的行名

library(pheatmap)
pheatmap(dat[cg,], show_colnames = F, show_rownames = F, 
         filename = 'all_cells_top_1000_sd.png')
dev.off()

# 对top1000的sd的基因集的表达矩阵，归一化，分组画图

if(T) {
    n = t(scale(t(dat[cg, ]))) # 因为scale默认对列操作，而我们需要对gene操作，所以t(data)将gene变成列进行scale
    n[n>2] = 2
    n[n<-2] = -2
    n[1:4, 1:4]
    
    ac = data.frame(g=group_list)
    rownames(ac) = colnames(n)
    pheatmap(n, show_rownames = F, show_colnames = F,
             annotation_col = ac,
             filename = 'all_cells_top_10000_sd_cutree1.png')
    dev.off()
}













