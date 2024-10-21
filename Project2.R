#Jayam Sutariya
#Project 2

library(magrittr)

### Question 1 ###

#Retrieving the data matrix and patient group information 
dataset <- "GSE19804"
gsets <- getGEO(dataset, GSEMatrix = T, getGPL = T)
gset <- gsets[[1]]
expr <- exprs(gset)

pdata <- pData(gset)
control <- rownames(pdata[grep("Lung Normal", pdata$title), ])
cancer <- rownames(pdata[grep("Lung Cancer", pdata$title), ])


### Question 2 ###

#Calculating difference in mean
cal_mean_diff <- function(x, cancer, control) {
  mean(x[cancer]) - mean(x[control])
}

#Calculating p-value
cal_p_value <- function(x, cancer, control) {
  t.test(x[cancer], x[control])$p.value
}

#Calculating t-score
cal_t_score <- function(x, cancer, control) {
  t.test(x[cancer], x[control])$statistic
}

#Calling the previous functions on the data
logFC <- apply(expr, 1, cal_mean_diff, cancer, control)
PValue <- apply(expr, MARGIN = 1, FUN = cal_p_value, cancer, control)
TScore <- apply(expr, MARGIN = 1, FUN = cal_t_score, cancer, control)

#Creating data frame and storing it in file
geneIds <- rownames(expr)

df <- data.frame(
  row.names = NULL,
  "GeneID" = geneIds,
  "PValue" = PValue,
  "TScore" = TScore,
  "LogFC" = logFC
)
saveRDS(df, file = "DE-results.rds")

### Question 3 ###

#Selecting DE genes, plotting as a volcano plot and saving it to file
pdf("volcano.pdf")
plot(
  x = df$LogFC,
  y = -log10(df$PValue),
  xlab = "logFC",
  ylab = "-log10(p-value)",
  main = "Volcano plot",
  col = ifelse(abs(df$LogFC) > 1 & df$PValue < 0.05, "red", "black")
)
abline(h = -log10(0.05), col = "red")
abline(v = -1, col = "blue")
abline(v = 1, col = "blue")
dev.off()


### Question 4 ###

#Calculating the t-SCORE between the groups and storing it in a data frame
cal_t_score <- function(x, cancer, control) {
  t.test(x[cancer], x[control])$statistic
}

TScore <- apply(expr, MARGIN = 1, FUN = cal_t_score, cancer, control)
t_OBSERVED <- data.frame(TScore)
saveRDS(t_OBSERVED, file = "t-observed.rds")


### Question 5 ###

#Calculating the empirical p-values for the t-SCORES and saving it to file
permuteExpr <- expr
t_NULL_DISTRIBUTION <- lapply(c(1:100), function(i) {
  message(i)
  colnames(permuteExpr) <- sample(colnames(expr))
  apply(permuteExpr, MARGIN = 1, FUN = cal_t_score, cancer, control)
}) %>%
  do.call(what = cbind) %>%
  as.data.frame()

pT <- lapply(rownames(t_OBSERVED), function(r) {
  sum(abs(t_NULL_DISTRIBUTION[r, ]) > abs(t_OBSERVED[r, ])) / ncol(t_NULL_DISTRIBUTION)
}) %>% unlist()
saveRDS(pT, file = "p-empirical-t-score.rds")


### Question 6 ###

#Helper function to calculate euclidean distance
cal_e_score <- function(x, cancer, control) {
  dist(rbind(mean(x[cancer]), mean(x[control])), method = "euclidean")
}

#Repeating steps 4 and 5 but with Euclidean distance
eSCORE <- apply(expr, MARGIN = 1, FUN = cal_e_score, cancer, control)
e_OBSERVED <- data.frame(eSCORE)

e_NULL_DISTRIBUTION <- lapply(c(1:100), function(i) {
  message(i)
  samples <- sample(ncol(expr), replace = F)
  controls <- samples[1:length(control)]
  cancers <- samples[(length(control) + 1):length(samples)]
  apply(
    expr,
    MARGIN = 1,
    FUN = cal_e_score,
    cancer = cancers,
    control = controls
  )
}) %>%
  do.call(what = cbind) %>%
  as.data.frame()

pE <- lapply(rownames(e_OBSERVED), function(r) {
  mean(abs(e_NULL_DISTRIBUTION[r, ]) > abs(e_OBSERVED[r, ]))
}) %>% unlist()
saveRDS(pE, file = "p-empirical-euclidean.rds")


### Question 7 ###

#Plotting pT and pE saving them to file
pdf("hist-pT.pdf")
hist(pT, breaks = 100, col = "red", main = "pT")
dev.off()

pdf("hist-pE.pdf")
hist(pE, breaks = 100, col = "red", main = "pE")
dev.off()


### Question 8 ###

#Calculating the correlation between the two p-value vectors
cor(pT, pE)