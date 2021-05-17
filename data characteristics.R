
##### Packages & Read in Data #####

Packages <- c("MCMCglmm", "splines", "dplyr", "plyr", "magrittr", "lme4", 
              "ggplot2", "tidyr", "matrixcalc", "abind", "gtable", "expm", "pROC")

lapply(Packages, library, character.only = TRUE)

# load data
load("minndat0.Rdata")

##### Filter Data #####

# assign dataset
dat <- minndat %>% filter(YearsSinceOnset >= 0 & YearsSinceOnset <= 40)
pidname <- "Patient.ID"
pidall <- unique(c(t(dat %>% select(pidname))))

variablenames <- c("pFVCt", "pDLCOt", "EFt", "RVSPt") 

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


# create spline variables
regdat.d <- regdat.d %>% merge(ddply(regdat.d, .(pid), summarize, last = last(time)), by = "pid") %>%
  mutate(btime = time - last, 
         sp10yr = ifelse(btime + 10 > 0, btime + 10, 0), 
         sp3yr = ifelse(btime + 3 > 0, btime + 3, 0))

regdat.d <- regdat.d %>% filter(!is.na(ACApos) & !is.na(RNAPolpos))
allid <- regdat.d$pid %>% unique

# remove rows where all 4 variables are missing
regdat.d <- regdat.d %>% filter(!(is.na(pFVCt) & is.na(pDLCOt) & is.na(EFt) & is.na(RVSPt)))


regdat.d$pFVCt[!is.na(regdat.d$pFVCt)] %>% length
regdat.d$pDLCOt[!is.na(regdat.d$pDLCOt)] %>% length
regdat.d$EFt[!is.na(regdat.d$EFt)] %>% length
regdat.d$RVSPt[!is.na(regdat.d$RVSPt)] %>% length


# RVSP cutoffs .
rvsp50 <- -1.33; rvsp45 <- -1.09

# FVC cutoffs
fvc70 <- -0.41; fvc60 <- -0.91

# EF cutoffs 
ef50 <- -1.75; ef35 <- -2.43


e1 <- sum(regdat.d$EFt[!is.na(regdat.d$EFt)] < ef50)
e2 <- sum(regdat.d$EFt[!is.na(regdat.d$EFt)] < ef35)

e3 <- sum(regdat.d$pFVCt[!is.na(regdat.d$pFVCt)] <= fvc70)
e4 <- sum(regdat.d$pFVCt[!is.na(regdat.d$pFVCt)] <= fvc60)

e5 <- sum(regdat.d$RVSPt[!is.na(regdat.d$RVSPt)] <= rvsp45)
e6 <- sum(regdat.d$RVSPt[!is.na(regdat.d$RVSPt)] <= rvsp50)

c("$EF < 50$",  "$EF < 35$", "$pFVC \le 70$", "$pFVC \le 60$", "$RVSP \ge 45$", "$RVSP \ge 50$")
cbind(e1, e2, e3, e4, e5, e6)






