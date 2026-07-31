rm(list = ls())

library('data.table')
library('dplyr')
library('clubSandwich')
library('stargazer')
library('Hmisc')
library('readstata13')
library('car')
library('miceadds')
library('multiwayvcov')

globdir <- "C:/Users/asankar/Dropbox (MIT)/HIM DATASETS/Analysis/ArunEstherPaper/"
datadir<-paste0(globdir,"data/")
tablesdir <- paste0(globdir,"tables/")


villagexmonth_level <- fread(paste0(datadir,"haryana_rct_Table3.csv"),header = TRUE, sep = ",", data.table = FALSE)
villagexmonth_level$fes <- group_indices(villagexmonth_level, id_phc, created_year, created_month)

villagexmonth_level$gossipplusSMS33 <- villagexmonth_level$gossip + villagexmonth_level$SMS33
villagexmonth_level$gossipminusSMS33 <- villagexmonth_level$gossip - villagexmonth_level$SMS33

villagexmonth_level$gossipplusSMS66 <- villagexmonth_level$gossip + villagexmonth_level$SMS66
villagexmonth_level$gossipminusSMS66 <- villagexmonth_level$gossip - villagexmonth_level$SMS66

data_level <- villagexmonth_level

outcomes <- c("shot_Penta1", "shot_Penta2", "shot_Penta3", "shot_Measles1", "nchildren")
outcomes_label <- c("Penta1 level", "Penta2 level", "Penta3 level", "Measles1 level", "Number of Children")

model_list <- list()
se_list <- list()
cmean_list <- c()
nna_list <- c()
num_clusters_list <- c()
p_val_gossip_random_list <- c()
p_val_gossip_trusted_list <- c()
p_val_gossip_trustgossip_list <- c()


for (y in outcomes) {
  formule <- as.formula("data_level[,y]~gossip+trusted+trustgossip+highslope+highflat+lowslope+lowflat+factor(fes)")
  
  model <- lm(formule, data = data_level)
  
  coeftests <- coef_test(model, vcov = "CR1", cluster = data_level[, "id_village_grp"], test = "naive-t")
  cluster_robust <- coeftests$SE
  
  
  print(coeftests[1:10,])
  
  nna <- sum(!is.na(data_level[,y]))
  num_clusters <- length(unique(data_level[,"id_village_grp"]))
  cmean <- round(mean(data_level[which(data_level$random == 1), gsub("ln_","",y)], na.rm = TRUE), 2)
  p_val_gossip_random <- round(coeftests['gossip','p_t'],3)
  
  formule_gossiptests <- as.formula("data_level[,y]~random+trusted+trustgossip+highslope+highflat+lowslope+lowflat+factor(fes)") #just to do additional linear hypothesis testing
  model_gossiptests <- lm(formule_gossiptests, data = data_level)
  coeftests_gossiptests <- coef_test(model_gossiptests, vcov = "CR1", cluster = data_level[, "id_village_grp"], test = "naive-t")
  p_val_gossip_trusted <- round(coeftests_gossiptests['trusted','p_t'],3)
  p_val_gossip_trustgossip <- round(coeftests_gossiptests['trustgossip','p_t'],3)
  
  
  
  model_list <- append(model_list, list(model))
  se_list <- append(se_list, list(cluster_robust))
  nna_list <- append(nna_list, nna)
  num_clusters_list <- append(num_clusters_list, num_clusters)
  cmean_list <- append(cmean_list, cmean)
  p_val_gossip_random_list <- append(p_val_gossip_random_list,p_val_gossip_random)
  p_val_gossip_trusted_list <- append(p_val_gossip_trusted_list, p_val_gossip_trusted)
  p_val_gossip_trustgossip_list <- append(p_val_gossip_trustgossip_list, p_val_gossip_trustgossip )
  
}

writeLines(capture.output(stargazer(model_list,
                                    se = se_list,
                                    multicolumn = FALSE,
                                    column.sep.width = "0.25pt",
                                    dep.var.labels=outcomes_label,
                                    font.size = "scriptsize",
                                    omit = c("fes", "Constant", "slope","flat"),
                                    omit.stat = "all",
                                    add.lines = list(c("Observations", nna_list),
                                                     c("Villages.", num_clusters_list),
                                                     c("Mean (Random Seeds)", cmean_list),
                                                     c("Gossip=Random (pval.)", p_val_gossip_random_list),
                                                     c("Gossip=Trusted (pval.)", p_val_gossip_trusted_list ),
                                                     c("Gossip=Trusted Gossip (pval.)", p_val_gossip_trustgossip_list)))
), 
paste0(tablesdir, "Table3_A.tex"))


model_list <- list()
se_list <- list()
cmean_list <- c()
nna_list <- c()
num_clusters_list <- c()
p_val_gossip_SMS33_list <- c()
p_val_gossip_SMS66_list <- c()
p_val_gossip_random_list <- c()
p_val_gossip_trusted_list <- c()
p_val_gossip_trustgossip_list <- c()


for (y in outcomes) {
  formule <- as.formula("data_level[,y]~gossip+trusted+trustgossip+SMS33+SMS66+highslope+highflat+lowslope+lowflat+factor(fes)")
  
  model <- lm(formule, data = data_level)
  
  coeftests <- coef_test(model, vcov = "CR1", cluster = data_level[, "id_village_grp"], test = "naive-t")
  cluster_robust <- coeftests$SE
  
  nna <- sum(!is.na(data_level[,y]))
  num_clusters <- length(unique(data_level[,"id_village_grp"]))
  cmean <- round(mean(data_level[which(data_level$random == 1), gsub("ln_","",y)], na.rm = TRUE), 2)
  p_val_gossip_random <- round(coeftests['gossip','p_t'],3)
  
  formule_gossiptests <- as.formula("data_level[,y]~random+trusted+trustgossip+SMS33+SMS66+highslope+highflat+lowslope+lowflat+factor(fes)") #just to do additional linear hypothesis testing
  model_gossiptests <- lm(formule_gossiptests, data = data_level)
  coeftests_gossiptests <- coef_test(model_gossiptests, vcov = "CR1", cluster = data_level[, "id_village_grp"], test = "naive-t")
  

  p_val_gossip_trusted <- round(coeftests_gossiptests['trusted','p_t'],3)
  p_val_gossip_trustgossip <- round(coeftests_gossiptests['trustgossip','p_t'],3)
  
  formule_gossiptests_SMS33 <-as.formula("data_level[,y]~gossipplusSMS33+gossipminusSMS33+SMS66+trusted+trustgossip+highslope+highflat+lowslope+lowflat+factor(fes)") #just to do additional linear hypothesis testing
  model_gossiptests_SMS33 <- lm(formule_gossiptests_SMS33, data = data_level)
  coeftests_gossiptests_SMS33 <- coef_test(model_gossiptests_SMS33, vcov = "CR1", cluster = data_level[, "id_village_grp"], test = "naive-t")
  
  p_val_gossip_SMS33 <- round(coeftests_gossiptests_SMS33['gossipminusSMS33','p_t'],3)
  
  
  formule_gossiptests_SMS66 <-as.formula("data_level[,y]~gossipplusSMS66+gossipminusSMS66+SMS33+trusted+trustgossip+highslope+highflat+lowslope+lowflat+factor(fes)") #just to do additional linear hypothesis testing
  model_gossiptests_SMS66 <- lm(formule_gossiptests_SMS66, data = data_level)
  coeftests_gossiptests_SMS66 <- coef_test(model_gossiptests_SMS66, vcov = "CR1", cluster = data_level[, "id_village_grp"], test = "naive-t")
  
  p_val_gossip_SMS66 <- round(coeftests_gossiptests_SMS66['gossipminusSMS66','p_t'],3)
  
  
  
  
  model_list <- append(model_list, list(model))
  se_list <- append(se_list, list(cluster_robust))
  nna_list <- append(nna_list, nna)
  num_clusters_list <- append(num_clusters_list, num_clusters)
  cmean_list <- append(cmean_list, cmean)
  p_val_gossip_SMS33_list <- append(p_val_gossip_SMS33_list, p_val_gossip_SMS33)
  p_val_gossip_SMS66_list <- append(p_val_gossip_SMS66_list, p_val_gossip_SMS66)
  p_val_gossip_random_list <- append(p_val_gossip_random_list,p_val_gossip_random)
  p_val_gossip_trusted_list <- append(p_val_gossip_trusted_list, p_val_gossip_trusted)
  p_val_gossip_trustgossip_list <- append(p_val_gossip_trustgossip_list, p_val_gossip_trustgossip )
  
}

writeLines(capture.output(stargazer(model_list,
                                    se = se_list,
                                    multicolumn = FALSE,
                                    column.sep.width = "0.25pt",
                                    dep.var.labels=outcomes_label,
                                    font.size = "scriptsize",
                                    omit = c("fes", "Constant", "slope","flat"),
                                    omit.stat = "all",
                                    add.lines = list(c("Observations", nna_list),
                                                     c("Villages.", num_clusters_list),
                                                     c("Mean (Random Seeds)", cmean_list),
                                                     c("Gossip=SMS33 (pval.)", p_val_gossip_SMS33_list),
                                                     c("Gossip=SMS66 (pval.)", p_val_gossip_SMS66_list ),
                                                     c("Gossip=Random (pval.)", p_val_gossip_random_list),
                                                     c("Gossip=Trusted (pval.)", p_val_gossip_trusted_list ),
                                                     c("Gossip=Trusted Gossip (pval.)", p_val_gossip_trustgossip_list)))
), 
paste0(tablesdir, "Table3_B.tex"))