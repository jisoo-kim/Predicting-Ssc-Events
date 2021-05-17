
Packages<- c("plyr", "magrittr")
lapply(Packages, library, character.only = TRUE)


# Setwd to current directory
setwd("~/Desktop/Git/Data & Preprocessing")


# CVSP
load("cvdfrvsp.Rdata")
load("cvdffvc.Rdata")
load("cvdfef.Rdata")

# logistic regression models results

load("efemppred.Rdata")
efpreddat <- preddat

load("rvspemppred.Rdata")
rvsppreddat <- preddat

load("fvcemppred.Rdata")
fvcpreddat <- preddat

# function to perform chi-square goodness of fit test
quantest <- function(dat, varname, obsvar){
  
  x <- dat %>% select(!!as.name(varname)) %>% unlist %>% as.numeric
  obsx <- dat %>% select(!!as.name(obsvar)) %>% unlist %>% as.numeric
  
  qx <- quantile(x, probs = seq(0, 1, 0.2))
  ind <- cut(x, qx, include.lowest = TRUE)
  
  seldat <- dat %>% select(varname, obsvar) %>% set_colnames(c("predicted", "observed"))
  
  quintab <- ddply(seldat, .(ind), summarise, 
                   observedprob = mean(observed), prectedprob = mean(predicted),
                   nobs = length(predicted), O = nobs*observedprob, E = nobs*prectedprob,
                   chisq = (O-E)^2/E)
  
  chisq <- sum(quintab$chisq)
  gft <- pchisq(q = chisq, df = 5)
  
  return(list(quintab, gft, chisq))
  
}  

# CVSP
X2rvsp45 <- quantest(dat = cvdf.rvsp, varname = "p45", obsvar = "rvsp45obs")
X2rvsp50 <- quantest(dat = cvdf.rvsp, varname = "p50", obsvar = "rvsp50obs")
X2fvc60 <- quantest(dat = cvdf.fvc, varname = "p60", obsvar = "fvc60obs")
X2fvc70 <- quantest(dat = cvdf.fvc, varname = "p70", obsvar = "fvc70obs")
X2ef35 <- quantest(dat = cvdf.ef, varname = "p35", obsvar = "ef35obs")
X2ef50 <- quantest(dat = cvdf.ef, varname = "p50", obsvar = "ef50obs")

# LR
ef35_1 <- quantest(dat = efpreddat, varname = "predprob1_35", obsvar = "Ief35")
ef35_2 <- quantest(dat = efpreddat, varname = "predprob2_35", obsvar = "Ief35")
ef35_3 <- quantest(dat = efpreddat, varname = "predprob3_35", obsvar = "Ief35")

ef50_1 <- quantest(dat = efpreddat, varname = "predprob1_50", obsvar = "Ief50")
ef50_2 <- quantest(dat = efpreddat, varname = "predprob2_50", obsvar = "Ief50")
ef50_3 <- quantest(dat = efpreddat, varname = "predprob3_50", obsvar = "Ief50")

rvsp45_1 <- quantest(dat = rvsppreddat, varname = "predprob1_45", obsvar = "Irvsp45")
rvsp45_2 <- quantest(dat = rvsppreddat, varname = "predprob2_45", obsvar = "Irvsp45")
rvsp45_3 <- quantest(dat = rvsppreddat, varname = "predprob3_45", obsvar = "Irvsp45")

rvsp50_1 <- quantest(dat = rvsppreddat, varname = "predprob1_50", obsvar = "Irvsp50")
rvsp50_2 <- quantest(dat = rvsppreddat, varname = "predprob2_50", obsvar = "Irvsp50")
rvsp50_3 <- quantest(dat = rvsppreddat, varname = "predprob3_50", obsvar = "Irvsp50")

fvc70_1 <- quantest(dat = fvcpreddat, varname = "predprob1_70", obsvar = "Ifvc70")
fvc70_2 <- quantest(dat = fvcpreddat, varname = "predprob2_70", obsvar = "Ifvc70")
fvc70_3 <- quantest(dat = fvcpreddat, varname = "predprob3_70", obsvar = "Ifvc70")

fvc60_1 <- quantest(dat = fvcpreddat, varname = "predprob1_60", obsvar = "Ifvc60")
fvc60_2 <- quantest(dat = fvcpreddat, varname = "predprob2_60", obsvar = "Ifvc60")
fvc60_3 <- quantest(dat = fvcpreddat, varname = "predprob3_60", obsvar = "Ifvc60")




X2tab <- rbind(c(paste0(round(X2ef50[[3]], 2), " (p = ", round(X2ef50[[2]], 2), ")"), 
                 paste0(round(ef50_1[[3]], 2), " (p = ", round(ef50_1[[2]], 2), ")"), 
                 paste0(round(ef50_2[[3]], 2), " (p = ", round(ef50_2[[2]], 2), ")"), 
                 paste0(round(ef50_3[[3]], 2), " (p = ", round(ef50_3[[2]], 2), ")")),
               
               c(paste0(round(X2ef35[[3]], 2), " (p = ", round(X2ef35[[2]], 2), ")"), 
                 paste0(round(ef35_1[[3]], 2), " (p = ", round(ef35_1[[2]], 2), ")"), 
                 paste0(round(ef35_2[[3]], 2), " (p = ", round(ef35_2[[2]], 2), ")"), 
                 paste0(round(ef35_3[[3]], 2), " (p = ", round(ef35_3[[2]], 2), ")")),
               
               c(paste0(round(X2rvsp45[[3]], 2), " (p = ", round(X2rvsp45[[2]], 2), ")"), 
                 paste0(round(rvsp45_1[[3]], 2), " (p = ", round(rvsp45_1[[2]], 2), ")"), 
                 paste0(round(rvsp45_2[[3]], 2), " (p = ", round(rvsp45_2[[2]], 2), ")"), 
                 paste0(round(rvsp45_3[[3]], 2), " (p = ", round(rvsp45_3[[2]], 2), ")")),
               
               c(paste0(round(X2rvsp50[[3]], 2), " (p = ", round(X2rvsp50[[2]], 2), ")"), 
                 paste0(round(rvsp50_1[[3]], 2), " (p = ", round(rvsp50_1[[2]], 2), ")"), 
                 paste0(round(rvsp50_2[[3]], 2), " (p = ", round(rvsp50_2[[2]], 2), ")"), 
                 paste0(round(rvsp50_3[[3]], 2), " (p = ", round(rvsp50_3[[2]], 2), ")")),
               
               c(paste0(round(X2fvc70[[3]], 2), " (p = ", round(X2fvc70[[2]], 2), ")"), 
                 paste0(round(fvc70_1[[3]], 2), " (p = ", round(fvc70_1[[2]], 2), ")"), 
                 paste0(round(fvc70_2[[3]], 2), " (p = ", round(fvc70_2[[2]], 2), ")"), 
                 paste0(round(fvc70_3[[3]], 2), " (p = ", round(fvc70_3[[2]], 2), ")")),
               
               c(paste0(round(X2fvc60[[3]], 2), " (p = ", round(X2fvc60[[2]], 2), ")"), 
                 paste0(round(fvc60_1[[3]], 2), " (p = ", round(fvc60_1[[2]], 2), ")"), 
                 paste0(round(fvc60_2[[3]], 2), " (p = ", round(fvc60_2[[2]], 2), ")"), 
                 paste0(round(fvc60_3[[3]], 2), " (p = ", round(fvc60_3[[2]], 2), ")")))

colnames(X2tab) <- c("CVSP", "LM1", "LM2", "LM3")
rownames(X2tab) <- c("EF50", "EF<35", "RVSP>=45", "RVSP>=50", "FVC<=70", "FVC<=60")

print(xtable::xtable(X2tab), comment = F)





