setwd('/Pathopticon_all_analyses/')
cont.tab <- read.csv('contingency_tables_CD.csv')
fisher.df <- data.frame(matrix(ncol = 5, nrow = length(cont.tab$Cell_type)))
colnames(fisher.df) <- c('pval', 'logpval', 'OR', 'OR_CI1', 'OR_CI2')
rownames(fisher.df) <- cont.tab$Cell_type

pvals <- rep(0, length(cont.tab$Cell_type))
ORs <- rep(0, length(cont.tab$Cell_type))
CI1s <- rep(0, length(cont.tab$Cell_type))
CI2s <- rep(0, length(cont.tab$Cell_type))

for (i in seq(length(cont.tab$Cell_type))) {
  f <- fisher.test(matrix(unlist(cont.tab[i, 2:5]), nrow=2, byrow = TRUE))
  pvals[i] <- f$p.value
  ORs[i] <- f$estimate
  CI1s[i] <- f$conf.int[1]
  CI2s[i] <- f$conf.int[2]
}

fisher.df$pval <- pvals
fisher.df$logpval <- -log10(pvals)
fisher.df$OR <- ORs
fisher.df$OR_CI1 <- CI1s
fisher.df$OR_CI2 <- CI2s

write.csv(fisher.df, '/Pathopticon_all_analyses/Pathopticon_intermediary_outputs/edge_overlap_pvals_CD.csv')