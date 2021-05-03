
##### Packages & Read in Data #####

Packages <- c("MCMCglmm", "splines", "dplyr", "plyr", "magrittr", "lme4", 
              "ggplot2", "tidyr", "matrixcalc", "abind", "gtable", "data.table", "pROC")

lapply(Packages, library, character.only = TRUE)

# load data
load("minndat0.Rdata")

##### Filter Data #####

# assign dataset
dat <- minndat %>% filter(YearsSinceOnset >= 0 & YearsSinceOnset <= 40)
pidname <- "Patient.ID"
pidall <- unique(c(t(dat %>% select(pidname))))

variablenames <- c("pFVCt", "pDLCOt", "RVSPt") 

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
              nrvsp =  sum(!is.na(RVSPt)))

# for m = 3
th <- 3
filterdat <- npid %>% filter(nfvc >= th & ndlco >= th & nrvsp >= th)
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

# find the most recent past fvc, dlco 
d[, prevFVC := c(NA, pFVCt[-.N]), by = pid]
d[, prevDLCO := c(NA, pDLCOt[-.N]), by = pid]

d[, prevFVCfill := nafill(prevDLCO, "locf"), by = pid]
d[, prevDLCOfill := nafill(prevDLCO, "locf"), by = pid]


d <- data.table(d %>% filter(!is.na(RVSPt)))

# find the most recent past rvsp 
d[, prevrvsp := c(NA, RVSPt[-.N]), by = pid]
d[, tprev := c(NA, btime[-.N]), by = pid]
d[, gap1 := btime - tprev]
d[gap1 > 3, prevrvsp := NA]

# find the second to the most recent past rvsp
d[, prevrvsp2 := c(NA, prevrvsp[-.N]), by = pid]
d[, tprev2 := c(NA, tprev[-.N]), by = pid]
d[, gap2 := tprev - tprev2]
d[(gap1 + gap2) > 3, prevrvsp2 := NA]

# set threshold
d <- d %>% mutate(Iprevrvsp = ifelse(is.na(prevrvsp), 0, prevrvsp), 
                  Iprevrvsp2 = ifelse(is.na(prevrvsp2), 0, prevrvsp2),
                  Irvsp45 = ifelse(RVSPt <= rvsp45, 1, 0),
                  Irvsp50 = ifelse(RVSPt <= rvsp50, 1, 0))

setDT(d)
d[, count45 := cumsum(Irvsp45), by = pid]
d[, prevcount45 := c(NA, count45[-.N]), by = pid]
d[, prevcount45 := ifelse(is.na(prevcount45), 0, prevcount45)]

d[, count50 := cumsum(Irvsp50), by = pid]
d[, prevcount50 := c(NA, count50[-.N]), by = pid]
d[, prevcount50 := ifelse(is.na(prevcount50), 0, prevcount50)]

##### fit model #####

vnames <- c("pid", "Irvsp45", "Irvsp50", "prevcount45", "prevcount50", "Iprevrvsp", "Iprevrvsp2", "prevFVCfill", "prevDLCOfill", dummyvarnames)

compd <- d %>% select(vnames)
#compd <- compd[complete.cases(compd), ]

formula1 <- as.formula(paste("Irvsp45 ~ ns(Iprevrvsp, 2) +", paste(dummyvarnames, collapse = " + ")))
formula2 <- as.formula(paste("Irvsp45 ~ ns(Iprevrvsp2, 2) + ns(Iprevrvsp, 2) +", paste(dummyvarnames, collapse = " + ")))
formula3 <- as.formula(paste("Irvsp45 ~ prevcount45 + ns(Iprevrvsp2, 2) + ns(Iprevrvsp, 2) +",
                             paste(dummyvarnames, collapse = " + ")))

formula1 <- as.formula(paste("Irvsp50 ~ ns(Iprevrvsp, 2) +", paste(dummyvarnames, collapse = " + ")))
formula2 <- as.formula(paste("Irvsp50 ~ ns(Iprevrvsp2, 2) + ns(Iprevrvsp, 2) +", paste(dummyvarnames, collapse = " + ")))
formula3 <- as.formula(paste("Irvsp50 ~ prevcount50 + ns(Iprevrvsp2, 2) + ns(Iprevrvsp, 2) +",
                             paste(dummyvarnames, collapse = " + ")))

##### K fold cross validation #####

# patient id
allid <- compd$pid %>% unique

# Randomly shuffle selected patient id
set.seed(123)
pid.cv <- sample(allid) # Randomly shuffle the data

# Create 5 equally size folds
K <- 5
folds <- cut(1:length(pid.cv), breaks = K, labels = FALSE)

preddat <- NULL

for (k in 1:K){
  
  # patient ids in k folds
  pid.cut <- pid.cv[which(folds == k)]
  pid.rest <-  pid.cv[which(folds != k)]
  
  pdat.cut <- compd %>% filter(pid %in% pid.cut)
  pdat.rest <- compd %>% filter(pid %in% pid.rest)
  
  fit1 <- glm(formula1, data = pdat.rest, family = "binomial") 
  fit2 <- glm(formula2, data = pdat.rest, family = "binomial") 
  fit3 <- glm(formula3, data = pdat.rest, family = "binomial")
  
  newpreddat <-  data.frame(pdat.cut, 
                            predprob1 = predict(fit1, newdata = pdat.cut, type = "response"),
                            predprob2 = predict(fit2, newdata = pdat.cut, type = "response"),
                            predprob3 = predict(fit3, newdata = pdat.cut, type = "response")
  )
  
  preddat <- rbind(preddat, newpreddat)
  
}

preddat <- preddat %>% arrange(pid)


##### ROC curves & AUC #####

pROC_fit1 <- roc(preddat$Irvsp50, preddat$predprob1,
                 smoothed = TRUE,
                 # arguments for ci
                 ci = TRUE, ci.alpha = 0.95, stratified = FALSE,
                 # arguments for plot
                 plot = TRUE, auc.polygon = T, max.auc.polygon = T, grid=F,
                 print.auc = TRUE, show.thres=TRUE)

pROC_fit2 <- roc(preddat$Irvsp45, preddat$predprob1,
                 smoothed = TRUE,
                 # arguments for ci
                 ci = TRUE, ci.alpha = 0.95, stratified = FALSE,
                 # arguments for plot
                 plot = TRUE, auc.polygon = T, max.auc.polygon = T, grid=F,
                 print.auc = TRUE, show.thres=TRUE)


