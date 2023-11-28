library(ranger); library(snpR)
args <- commandArgs(TRUE) 
data <- readRDS(as.character(args[1]))
out <- as.character(args[2])
i <- as.numeric(args[3])

num.trees <- 1000000
set.seed(4+i)

data_sn <- format_snps(data, output="sn")
data_sn <- t(data_sn[,-c(1:3)])
data_sites <- as.factor(sample.meta(data)$Site_Code)


rf_data <- ranger(x=data_sn[-i,], importance="permutation", keep.inbag = TRUE, y=data_sites[-i], num.trees = num.trees, mtry = ncol(data_sn))

eval <- predict(rf_data, data_sn[i,,drop=FALSE], num.trees = num.trees)

err  <- forestError::quantForestError(rf_data, # the forest 
                                      X.train = data_sn[-i,], # the training data
                                      X.test =  data_sn[i,,drop=FALSE], # the test data
                                      Y.train = data_sites[-i])  # the classifications of the training data

fwrite(data.frame(sample.meta(data)$Site_Code[i],eval$predictions,err), file = paste0(out, "_", i, ".txt"), row.names = FALSE, col.names = FALSE, quote = FALSE, sep = "\t") 


### Above run on bell 


### is this different from "random"? 
data <- read.csv("gt_ranger_loo_adult.csv")
addmargins(table(data$site, data$pred))[-3,]

hist(data$mcr)

n = 10000
results <- numeric(n)
for(i in 1:n){
  results[i] <- sum(sample(unique(data$site), nrow(data), replace=TRUE) == data$site)/nrow(data)
}
prop <- sum(data$site == data$pred) / nrow(data)
sum(results >= prop)/n

