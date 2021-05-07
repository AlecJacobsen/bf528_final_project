install.packages('ggplot2')
library(ggplot2)

rm(list = ls())
#load in cufflinks output
genes = data.frame(read.table('/projectnb/bf528/users/lava_lamp/project_4/Alec_final/P0_1_cufflinks/genes.fpkm_tracking', header=T))

pdf('/projectnb/bf528/users/lava_lamp/project_4/Alec_final/FPKM_hist.pdf')
ggplot(genes[genes$FPKM > 0.0001,], aes(log(FPKM))) + geom_histogram(bins = 100, color = 'black') + 
  theme_bw() + ggtitle('Distribution of FPKM for All Genes') + ylab('Count')
dev.off()
