
if(F) {
    set.seed(0)
    
    library(pheatmap)
    library(Rtsne)
    library(ggfortify)
    library(mvtnorm)
    
    
    # 模拟数据理解tsne
    
    # 创建两个正态分布的随机表达矩阵，这是无法区分开来的
    if(T) {
        ng = 500
        nc = 20
        a1 = rnorm(ng * nc); dim(a1) = c(ng, nc)
        a2 = rnorm(ng * nc); dim(a2) = c(ng, nc)
        a3 = cbind(a1, a2)
        colnames(a3) = c(paste0('cell_01_', 1:nc),
                         paste0('cell_02_', 1:nc))
        
        rownames(a3) = paste('gene_', 1:ng, sep = '')
        pheatmap(a3)
        dist(a3) # 这种默认的用就是行与行之间的，这里就是基因与基因之间的关系。人们关心的是样本之间的
        # 所以我们下面在进行pca时，需要转置一下
        a3 = t(a3); dim(a3)
        pca_dat = prcomp(a3, scale. = T) # prcomp() 主成分分析
        p = autoplot(pca_dat, ) + 
            theme_classic() +
            ggtitle('PCA plot') # 做出PCA可以得到此处随机生成的数据细胞是无法被区分开的
        df = cbind(as.data.frame(a3), group=c(rep('b1', 20), rep('b2', 20)))
        pca_result <- prcomp(df[, 1:(ncol(df) - 1)])
        autoplot(pca_result, data = df, colour = 'group') + 
            theme_bw()
        
        set.seed(42)
        tsne_out = Rtsne(a3, pca = F, perplexity = 10, theta = 0.0)
        tsnes = tsne_out$Y
        colnames(tsnes) = c('tSNE1', 'tSNE2')
        group = c(rep('b1', 20), rep('b2', 20))
        tsnes = as.data.frame(tsnes)
        tsnes$group = group
        ggplot(tsnes, aes(x = tSNE1, y = tSNE2)) + geom_point(aes(col=group))
         
    }
    
    
    # 同样的正态分布随机表达矩阵，但其中部分细胞+3，可以区分开来。
    if(T) {
        ng = 500
        nc = 20
        a1 = rnorm(ng * nc); dim(a1) = c(ng, nc)
        a2 = rnorm(ng * nc)+3; dim(a2) = c(ng, nc)
        a3 = cbind(a1, a2)
        colnames(a3) = c(paste0('cell_01_', 1:nc),
                         paste0('cell_02_', 1:nc))
        
        rownames(a3) = paste('gene_', 1:ng, sep = '')
        pheatmap(a3)
        dist(a3) # 这种默认的用就是行与行之间的，这里就是基因与基因之间的关系。人们关心的是样本之间的
        # 所以我们下面在进行pca时，需要转置一下
        a3 = t(a3); dim(a3)
        pca_dat = prcomp(a3, scale. = T) # prcomp() 主成分分析
        p = autoplot(pca_dat, ) + 
            theme_classic() +
            ggtitle('PCA plot') # 做出PCA可以得到此处随机生成的数据细胞是无法被区分开的
        df = cbind(as.data.frame(a3), group=c(rep('b1', 20), rep('b2', 20)))
        pca_result <- prcomp(df[, 1:(ncol(df) - 1)])
        autoplot(pca_result, data = df, colour = 'group') + 
            theme_bw()
        
        set.seed(42)
        tsne_out = Rtsne(a3, pca = F, perplexity = 10, theta = 0.0)
        tsnes = tsne_out$Y
        colnames(tsnes) = c('tSNE1', 'tSNE2')
        group = c(rep('b1', 20), rep('b2', 20))
        tsnes = as.data.frame(tsnes)
        tsnes$group = group
        ggplot(tsnes, aes(x = tSNE1, y = tSNE2)) + geom_point(aes(col=group))
        
    }
}

# 真实数据演练
rm(list = ls())
load(file = 'input_rpkm.Rdata')

dat_back = dat

dat = t(dat)
dat[1:2, 1:4]
dat = as.data.frame(dat)

colnames(dat)
library(stringr)
plate = str_split(colnames(dat), '_', simplify = T)[,3]

dat = cbind(dat, plate)
dat[1:4, 1:4]
table(dat$plate)

library(FactoMineR)
library(factoextra)

dat.pca = PCA(dat[,-ncol(dat)], graph = F)
fviz_pca_ind(dat.pca, # repel =T,
             geom.ind = "point", # 只显示点，不显示文字
             col.ind = dat$plate, # 按分组上色
             #palette = c("#00AFBB", "#E7B800"),
             addEllipses = TRUE, # 添加晕环
             legend.title = "Groups")
ggsave('all_cells_PCA_by_plate.png')





