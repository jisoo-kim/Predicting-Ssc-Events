
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

d <- data.table(regdat.d)

# find the most recent past rvsp, dlco, ef
d[, prevEF := c(NA, EFt[-.N]), by = pid]
d[, prevDLCO := c(NA, pDLCOt[-.N]), by = pid]
d[, prevRVSP := c(NA, RVSPt[-.N]), by = pid]

d[, prevEFfill := nafill(prevEF, "locf"), by = pid]
d[, prevDLCOfill := nafill(prevDLCO, "locf"), by = pid]
d[, prevRVSPfill := nafill(prevRVSP, "locf"), by = pid]

d <- data.table(d %>% filter(!is.na(pFVCt)))

# find the most recent past EF 
d[, prevFVC := c(NA, pFVCt[-.N]), by = pid]
d[, tprev := c(NA, btime[-.N]), by = pid]
d[, gap1 := btime - tprev]
d[gap1 > 3, prevFVC := NA]

# find the second to the most recent past EF
d[, prevFVC2 := c(NA, prevFVC[-.N]), by = pid]
d[, tprev2 := c(NA, tprev[-.N]), by = pid]
d[, gap2 := tprev - tprev2]
d[(gap1 + gap2) > 3, prevFVC2 := NA]

d <- d %>% mutate(Ifvc60 = ifelse(pFVCt <= fvc60, 1, 0),
                  Ifvc70 = ifelse(pFVCt <= fvc70, 1, 0))

setDT(d)
d[, count70 := cumsum(Ifvc70), by = pid]
d[, count60 := cumsum(Ifvc60), by = pid]

d[, prevcount70 := c(NA, count70[-.N]), by = pid]
d[, prevcount60 := c(NA, count60[-.N]), by = pid]



##### impute missing values #####

vnames <- c("pid", "Ifvc70", "Ifvc60", "prevFVC", "prevFVC2", 
            "prevEFfill", "prevDLCOfill", "prevRVSPfill",
            "prevcount70", "prevcount60", dummyvarnames)

compd <- d %>% select(all_of(vnames))

# dataset with NAs
impd <- compd %>% select(-c(pid, Ifvc70, Ifvc60))

mice_imp <- mice(data = impd, m = 5, seed = 123)
mdat <- complete(mice_imp)

impdat <- data.frame(compd %>% select(pid, Ifvc70, Ifvc60), mdat)



##### fit model #####

formula1_70 <- as.formula(paste("Ifvc70 ~ ns(prevFVC, 2) + prevEFfill + prevDLCOfill + prevRVSPfill +", 
                                paste(dummyvarnames, collapse = " + ")))
formula2_70  <- as.formula(paste("Ifvc70 ~ ns(prevFVC2, 2) + ns(prevFVC, 2) + prevEFfill + 
                                 prevDLCOfill + prevRVSPfill +", paste(dummyvarnames, collapse = " + ")))
formula3_70  <- as.formula(paste("Ifvc70 ~ prevcount70 + ns(prevFVC2, 2) + ns(prevFVC, 2) +
                                 prevEFfill + prevDLCOfill + prevRVSPfill +",
                                 paste(dummyvarnames, collapse = " + ")))

formula1_60 <- as.formula(paste("Ifvc60 ~ ns(prevFVC, 2) + prevEFfill + prevDLCOfill + prevRVSPfill +",
                                paste(dummyvarnames, collapse = " + ")))
formula2_60 <- as.formula(paste("Ifvc60 ~ ns(prevFVC2, 2) + ns(prevFVC, 2) + prevEFfill + 
                                prevDLCOfill + prevRVSPfill +", paste(dummyvarnames, collapse = " + ")))
formula3_60 <- as.formula(paste("Ifvc60 ~ prevcount60 + ns(prevFVC2, 2) + ns(prevFVC, 2) + 
                                prevEFfill + prevDLCOfill + prevRVSPfill +",
                                paste(dummyvarnames, collapse = " + ")))


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
  
  fit1_70 <- glm(formula1_70, data = pdat.rest, family = "binomial") 
  fit2_70 <- glm(formula2_70, data = pdat.rest, family = "binomial") 
  fit3_70 <- glm(formula3_70, data = pdat.rest, family = "binomial")
  
  fit1_60 <- glm(formula1_60, data = pdat.rest, family = "binomial")
  fit2_60 <- glm(formula2_60, data = pdat.rest, family = "binomial")
  fit3_60 <- glm(formula3_60, data = pdat.rest, family = "binomial")
  
  newpreddat <-  data.frame(pdat.cut, 
                            predprob1_70 = predict(fit1_70, newdata = pdat.cut, type = "response"),
                            predprob2_70 = predict(fit2_70, newdata = pdat.cut, type = "response"),
                            predprob3_70 = predict(fit3_70, newdata = pdat.cut, type = "response"),
                            
                            predprob1_60 = predict(fit1_60, newdata = pdat.cut, type = "response"),
                            predprob2_60 = predict(fit2_60, newdata = pdat.cut, type = "response"),
                            predprob3_60 = predict(fit3_60, newdata = pdat.cut, type = "response"))
  
  preddat <- rbind(preddat, newpreddat)
  
}

preddat <- preddat %>% arrange(pid)

#save(preddat, file = "fvcemppred.Rdata")




##### ROC curves & AUC #####

par(mfrow = c(1, 3))
fvc70_fit1 <- roc(preddat$Ifvc70, preddat$predprob1_70,
                 smoothed = TRUE,
                 # arguments for ci
                 ci = TRUE, ci.alpha = 0.95, stratified = FALSE,
                 # arguments for plot
                 plot = TRUE, auc.polygon = T, max.auc.polygon = T, grid=F,
                 print.auc = TRUE, show.thres=TRUE, main = "LM1")

fvc70_fit2 <- roc(preddat$Ifvc70, preddat$predprob2_70,
                 smoothed = TRUE,
                 # arguments for ci
                 ci = TRUE, ci.alpha = 0.95, stratified = FALSE,
                 # arguments for plot
                 plot = TRUE, auc.polygon = T, max.auc.polygon = T, grid=F,
                 print.auc = TRUE, show.thres=TRUE, main = "LM2")


fvc70_fit3 <- roc(preddat$Ifvc70, preddat$predprob3_70,
                 smoothed = TRUE,
                 # arguments for ci
                 ci = TRUE, ci.alpha = 0.95, stratified = FALSE,
                 # arguments for plot
                 plot = TRUE, auc.polygon = T, max.auc.polygon = T, grid=F,
                 print.auc = TRUE, show.thres=TRUE, main = "LM3")


fvc60_fit1 <- roc(preddat$Ifvc60, preddat$predprob1_60,
                  smoothed = TRUE,
                  # arguments for ci
                  ci = TRUE, ci.alpha = 0.95, stratified = FALSE,
                  # arguments for plot
                  plot = TRUE, auc.polygon = T, max.auc.polygon = T, grid=F,
                  print.auc = TRUE, show.thres=TRUE, main = "LM1")

fvc60_fit2 <- roc(preddat$Ifvc60, preddat$predprob2_60,
                  smoothed = TRUE,
                  # arguments for ci
                  ci = TRUE, ci.alpha = 0.95, stratified = FALSE,
                  # arguments for plot
                  plot = TRUE, auc.polygon = T, max.auc.polygon = T, grid=F,
                  print.auc = TRUE, show.thres=TRUE, main = "LM2")


fvc60_fit3 <- roc(preddat$Ifvc60, preddat$predprob3_60,
                  smoothed = TRUE,
                  # arguments for ci
                  ci = TRUE, ci.alpha = 0.95, stratified = FALSE,
                  # arguments for plot
                  plot = TRUE, auc.polygon = T, max.auc.polygon = T, grid=F,
                  print.auc = TRUE, show.thres=TRUE, main = "LM3")
