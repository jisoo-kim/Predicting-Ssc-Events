
##### Packages & Read in Data #####

Packages <- c("MCMCglmm", "splines", "dplyr", "plyr", "magrittr", "lme4", "mice",
              "ggplot2", "tidyr", "matrixcalc", "abind", "gtable", "data.table", "pROC")

lapply(Packages, library, character.only = TRUE)


# Setwd to current directory
setwd("~/Desktop/Git/Data & Preprocessing")

# load data
load("minndat0.Rdata")

##### Filter Data #####

# assign dataset
dat <- minndat %>% filter(YearsSinceOnset >= 0 & YearsSinceOnset <= 40)
pidname <- "Patient.ID"
pidall <- unique(c(t(dat %>% select(pidname))))

variablenames <- c("pFVCt", "pDLCOt", "EFt", "RVSPt") 

# EF cutoffs 
ef50 <- -1.75; ef35 <- -2.43

# RVSP cutoffs 
rvsp50 <- -1.33; rvsp45 <- -1.09

# FVC cutoffs
fvc70 <- -0.41; fvc60 <- -0.91

baselinevarnames <- c("AgeOnset", "race", "Type", "Sex", "ACA", "RNAPol", "SCL70")
datename <- "YearsSinceOnset"
timename <- "YearsSinceFirstObs"

# remove row with missing time & order by time
regdat <- dat %>% filter(complete.cases(!!as.symbol(datename))) %>%
  dplyr::arrange(!!as.symbol(pidname), !!as.symbol(datename)) %>%
  select(c(pidname, datename, timename, variablenames, baselinevarnames)) %>%
  set_colnames(c("pid", "time", "ref.time", variablenames, baselinevarnames))


# select patients with more than th observations
npid <- ddply(regdat, .(pid), summarize,
              nfvc = sum(!is.na(pFVCt)), ndlco = sum(!is.na(pDLCOt)),
              nef =  sum(!is.na(EFt)), nrvsp =  sum(!is.na(RVSPt)))

# for m = 4
th <- 2
filterdat <- npid %>% filter(nfvc >= th & ndlco >= th & nef >= th & nrvsp >= th)
filterid <- filterdat$pid
regdat.d <- regdat %>% filter(pid %in% filterid)

# create dataset with dummy variables
dummyvarnames <- c("race_AA", "Type_Diffuse", "Type_Sine", "Sex_Male", "ACApos", "RNAPolpos", "SCL70pos")
regdat.d <- fastDummies::dummy_cols(regdat.d, remove_most_frequent_dummy = T) %>%
  select(-c(baselinevarnames[-1])) %>% set_colnames(c("pid", "time", "ref.time", variablenames, "AgeOnset", dummyvarnames))

allid <- regdat.d$pid %>% unique

# create spline variables
regdat.d <- regdat.d %>% merge(ddply(regdat.d, .(pid), summarize, last = last(time)), by = "pid") %>%
  mutate(btime = time - last, 
         sp10yr = ifelse(btime + 10 > 0, btime + 10, 0), 
         sp3yr = ifelse(btime + 3 > 0, btime + 3, 0))


##### manipulate data to fit model #####

# remove rows where all 4 variables are missing
regdat.d <- regdat.d %>% filter(!(is.na(pFVCt) & is.na(pDLCOt) & is.na(EFt) & is.na(RVSPt)))

# find the most recent past fvc, dlco, rvsp
d <- data.table(regdat.d)

d[, prevFVC := c(NA, pFVCt[-.N]), by = pid]
d[, prevDLCO := c(NA, pDLCOt[-.N]), by = pid]
d[, prevRVSP := c(NA, RVSPt[-.N]), by = pid]

d[, prevFVCfill := nafill(prevFVC, "locf"), by = pid]
d[, prevDLCOfill := nafill(prevDLCO, "locf"), by = pid]
d[, prevRVSPfill := nafill(prevRVSP, "locf"), by = pid]

d <- data.table(d %>% filter(!is.na(EFt)))

# find the most recent past EF 
d[, prevEF := c(NA, EFt[-.N]), by = pid]
d[, tprev := c(NA, btime[-.N]), by = pid]
d[, gap1 := btime - tprev]
d[gap1 > 3, prevEF := NA]

# find the second to the most recent past EF
d[, prevEF2 := c(NA, prevEF[-.N]), by = pid]
d[, tprev2 := c(NA, tprev[-.N]), by = pid]
d[, gap2 := tprev - tprev2]
d[(gap1 + gap2) > 3, prevEF2 := NA]

d <- d %>% mutate(#IprevEF = ifelse(is.na(prevEF), 0, prevEF), 
                  #IprevEF2 = ifelse(is.na(prevEF2), 0, prevEF2),
                  Ief35 = ifelse(EFt < ef35, 1, 0),
                  Ief50 = ifelse(EFt < ef50, 1, 0))

setDT(d)
d[, count50 := cumsum(Ief50), by = pid]
d[, count35 := cumsum(Ief35), by = pid]

d[, prevcount50 := c(NA, count50[-.N]), by = pid]
d[, prevcount35 := c(NA, count35[-.N]), by = pid]

##### impute missing values #####

vnames <- c("pid", "Ief50", "Ief35", "prevEF", "prevEF2", 
            "prevFVCfill", "prevDLCOfill", "prevRVSPfill",
            "prevcount50", "prevcount35", dummyvarnames)

compd <- d %>% select(all_of(vnames))

# dataset with NAs
impd <- compd %>% select(-c(pid, Ief50, Ief35))

mice_imp <- mice(data = impd, m = 5, seed = 123)
mdat <- complete(mice_imp)

impdat <- data.frame(compd %>% select(pid, Ief50, Ief35), mdat)




##### fit model #####

formula1_50 <- as.formula(paste("Ief50 ~ ns(prevEF, 2) + prevFVCfill + prevDLCOfill + prevRVSPfill +", 
                                paste(dummyvarnames, collapse = " + ")))
formula2_50  <- as.formula(paste("Ief50 ~ ns(prevEF2, 2) + ns(prevEF, 2) + prevFVCfill + 
                                 prevDLCOfill + prevRVSPfill +", paste(dummyvarnames, collapse = " + ")))
formula3_50  <- as.formula(paste("Ief50 ~ prevcount50 + ns(prevEF2, 2) + ns(prevEF, 2) +
                                 prevFVCfill + prevDLCOfill + prevRVSPfill +",
                                 paste(dummyvarnames, collapse = " + ")))

formula1_35 <- as.formula(paste("Ief35 ~ ns(prevEF, 2) + prevFVCfill + prevDLCOfill + prevRVSPfill +",
                                paste(dummyvarnames, collapse = " + ")))
formula2_35 <- as.formula(paste("Ief35 ~ ns(prevEF2, 2) + ns(prevEF, 2) + prevFVCfill + 
                                prevDLCOfill + prevRVSPfill +", paste(dummyvarnames, collapse = " + ")))
formula3_35 <- as.formula(paste("Ief35 ~ prevcount35 + ns(prevEF2, 2) + ns(prevEF, 2) + 
                                prevFVCfill + prevDLCOfill + prevRVSPfill +",
                                paste(dummyvarnames, collapse = " + ")))

#################################################################################

install.packages("glmnet")
library(glmnet)

fit <- glmnet(x = , y = , family = "binomial")
head(test)


# k = 4
test <- pdat.rest
fit <- glm(formula2_35, data = test, family = "binomial")

t1 <- predict(fit)
sum(exp(t1)/(1+exp(t1)))


test$Ief35 %>% sum
#1: glm.fit: fitted probabilities numerically 0 or 1 occurred 

#################################################################################

##### K fold cross validation #####

# patient id
allid <- impdat$pid %>% unique

# Randomly shuffle selected patient id
set.seed(123)
pid.cv <- sample(allid) # Randomly shuffle the data

# Create 5 equally size folds
K <- 5
folds <- cut(1:length(pid.cv), breaks = K, labels = FALSE)

preddat <- NULL

for (k in 1:K){

  # patient ikds in k folds
  pid.cut <- pid.cv[which(folds == k)]
  pid.rest <-  pid.cv[which(folds != k)]
  
  pdat.cut <- impdat %>% filter(pid %in% pid.cut)
  pdat.rest <- impdat %>% filter(pid %in% pid.rest)
  
  fit1_50 <- glm(formula1_50, data = pdat.rest, family = "binomial") 
  fit2_50 <- glm(formula2_50, data = pdat.rest, family = "binomial") 
  fit3_50 <- glm(formula3_50, data = pdat.rest, family = "binomial")
  
  fit1_35 <- glm(formula1_35, data = pdat.rest, family = "binomial")
  fit2_35 <- glm(formula2_35, data = pdat.rest, family = "binomial")
  fit3_35 <- glm(formula3_35, data = pdat.rest, family = "binomial")
  
  newpreddat <-  data.frame(pdat.cut, 
                            predprob1_50 = predict(fit1_50, newdata = pdat.cut, type = "response"),
                            predprob2_50 = predict(fit2_50, newdata = pdat.cut, type = "response"),
                            predprob3_50 = predict(fit3_50, newdata = pdat.cut, type = "response"),
                            
                            predprob1_35 = predict(fit1_35, newdata = pdat.cut, type = "response"),
                            predprob2_35 = predict(fit2_35, newdata = pdat.cut, type = "response"),
                            predprob3_35 = predict(fit3_35, newdata = pdat.cut, type = "response"))
  
  preddat <- rbind(preddat, newpreddat)
  
}

preddat <- preddat %>% arrange(pid)

#save(preddat, file = "efemppred.Rdata")


##### ROC curves & AUC #####

par(mfrow = c(1, 3))
ef50_fit1 <- roc(preddat$Ief50, preddat$predprob1_50,
                 smoothed = TRUE,
                 # arguments for ci
                 ci = TRUE, ci.alpha = 0.95, stratified = FALSE,
                 # arguments for plot
                 plot = TRUE, auc.polygon = T, max.auc.polygon = T, grid=F,
                 print.auc = TRUE, show.thres=TRUE, main = "LM1")

ef50_fit2 <- roc(preddat$Ief50, preddat$predprob2_50,
                 smoothed = TRUE,
                 # arguments for ci
                 ci = TRUE, ci.alpha = 0.95, stratified = FALSE,
                 # arguments for plot
                 plot = TRUE, auc.polygon = T, max.auc.polygon = T, grid=F,
                 print.auc = TRUE, show.thres=TRUE, main = "LM2")


ef50_fit3 <- roc(preddat$Ief50, preddat$predprob3_50,
                 smoothed = TRUE,
                 # arguments for ci
                 ci = TRUE, ci.alpha = 0.95, stratified = FALSE,
                 # arguments for plot
                 plot = TRUE, auc.polygon = T, max.auc.polygon = T, grid=F,
                 print.auc = TRUE, show.thres=TRUE, main = "LM3")


ef35_fit1 <- roc(preddat$Ief35, preddat$predprob1_35,
                 smoothed = TRUE,
                 # arguments for ci
                 ci = TRUE, ci.alpha = 0.95, stratified = FALSE,
                 # arguments for plot
                 plot = TRUE, auc.polygon = T, max.auc.polygon = T, grid=F,
                 print.auc = TRUE, show.thres=TRUE, main = "LM1")

ef35_fit2 <- roc(preddat$Ief35, preddat$predprob2_35,
                 smoothed = TRUE,
                 # arguments for ci
                 ci = TRUE, ci.alpha = 0.95, stratified = FALSE,
                 # arguments for plot
                 plot = TRUE, auc.polygon = T, max.auc.polygon = T, grid=F,
                 print.auc = TRUE, show.thres=TRUE, main = "LM2")

ef35_fit3 <- roc(preddat$Ief35, preddat$predprob3_35,
                 smoothed = TRUE,
                 # arguments for ci
                 ci = TRUE, ci.alpha = 0.95, stratified = FALSE,
                 # arguments for plot
                 plot = TRUE, auc.polygon = T, max.auc.polygon = T, grid=F,
                 print.auc = TRUE, show.thres=TRUE, main = "LM3")
