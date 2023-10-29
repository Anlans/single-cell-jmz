rm(list = ls())
options(stringsAsFactors = F)


load('input_rpkm.Rdata')
a[1:4, 1:4]
head(df)

dat[1:4, 1:4]
exprSet = dat

mean_per_gene = apply(exprSet, 1, mean, na.rm = T)
sd_per_gene = apply(exprSet, 1, sd, na.rm = T)
mad_per_gene = apply(exprSet, 1, mad, na.rm = T)

# 构建一个df存放结果
cv_per_gene = data.frame(mean = mean_per_gene,
                         sd = sd_per_gene,
                         mad = mad_per_gene,
                         cv = sd_per_gene / mean_per_gene)
rownames(cv_per_gene)
head(cv_per_gene)

with(cv_per_gene, plot(log10(mean), log10(cv^2)))

cv_per_gene$log10cv2 = log10(cv_per_gene$cv^2)
cv_per_gene$log10mean = log10(cv_per_gene$mean)

library(ggpubr)
cv_per_gene = cv_per_gene[cv_per_gene$log10mean < 4 & cv_per_gene$log10mean > 0, ]
ggscatter(cv_per_gene, x = 'log10mean', y = 'log10cv2',
          color = 'black', shape = 16, size = 1, add = 'loess',
          xlab = 'log10(mean)RPKM', ylab = 'log10(cv^2)')

if(F) {
    
}











