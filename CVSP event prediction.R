
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


# data for estimating regression parameters 

regdat.d <- data.frame(cbind(ns(regdat.d$time, knots = c(10, 30), Boundary.knots= c(0, 40))),
                       regdat.d, int = 1)

colnames(regdat.d)[1:3] <- c("ns1", "ns2", "ns3")

FEcolnames <- c("int", "ns1", "ns2", "ns3", "AgeOnset", "race_AA", "Type_Diffuse", "Type_Sine", "Sex_Male",
                "ACApos", "RNAPolpos", "SCL70pos")

REcolnames <- c("int", "btime", "sp10yr", "sp3yr")

datfvc <- regdat.d %>% filter(complete.cases(!!as.symbol(variablenames[1])))
datdlco <- regdat.d %>% filter(complete.cases(!!as.symbol(variablenames[2])))
datef <- regdat.d %>% filter(complete.cases(!!as.symbol(variablenames[3])))
datrvsp <- regdat.d %>% filter(complete.cases(!!as.symbol(variablenames[4])))




##### K fold cross validation #####

# Randomly shuffle selected patient id
set.seed(123)
pid.cv <- sample(allid) # Randomly shuffle the data

# Create 5 equally size folds
K <- 5
folds <- cut(1:length(pid.cv), breaks = K, labels = FALSE)


##### model #####

# Set prior for combined model when m = 4

nmeas <- length(variablenames)

Mode <- diag(c(rep(1, 4), rep(0.005, 8), rep(0.00005, 4)))
p <- dim(Mode)[2]
nu <- 16

phi <- Mode * (nu + p + 1)
V <- phi/nu

# prior for R
p.R <- 4
Mode.R <- diag(4)

nu.R <- 4

phi.R <- Mode.R * (nu.R + p.R + 1)
V.R <- phi.R/nu.R


nitt <- 15000; burnin <- 2000

combmcmc <- MCMCglmm(as.formula(paste("cbind(", paste(variablenames, collapse = ","),
                                      ") ~ -1 + trait + trait:(ns(time, knots = c(10, 30), Boundary.knots= c(0, 40)) +
                                      AgeOnset + race_AA + Type_Diffuse + Type_Sine + Sex_Male + ACApos + RNAPolpos + SCL70pos)")),
                     random = ~ us(trait + trait:btime + trait:sp10yr + trait:sp3yr):pid, pr = TRUE,
                     rcov = ~ us(trait):units,
                     prior = list(R = list(V = V.R, nu = nu.R),
                                 G = list(G1 = list(V = V, nu = nu))),
                     family = rep("gaussian", nmeas),
                     nitt = nitt, burnin = burnin, data = regdat.d)


##### Fixed Effects estimates #####

# ordering of FE estimates

nfe <- combmcmc$Fixed$nfl

solmat <- combmcmc$Sol
sol.postmean <- apply(solmat[, 1:nfe], 2, mean)

cnames <- names(sol.postmean)
colorder <- c(cnames[startsWith(cnames, "traitpFVCt")], 
              cnames[startsWith(cnames, "traitpDLCOt")],
              cnames[startsWith(cnames, "traitEFt")],
              cnames[startsWith(cnames, "traitRVSPt")]) 

beta <- sol.postmean[order(match(names(sol.postmean), colorder))] 


##### covariance matrices #####

# Varcov matrix from mcmcglmm: R and G

covdat <- apply(combmcmc$VCV, 2, mean)

tmat <- covdat[endsWith(names(covdat), "pid")]
Gnames <- matrix(names(tmat), sqrt(length(tmat))) %>% diag

Gvec <- 1:length(Gnames)
ordervec <- c(Gvec[startsWith(Gnames, "traitpFVCt")], 
              Gvec[startsWith(Gnames, "traitpDLCOt")],
              Gvec[startsWith(Gnames, "traitEFt")],
              Gvec[startsWith(Gnames, "traitRVSPt")])

G <- matrix(tmat, sqrt(length(tmat)))
G2 <- G[ordervec, ordervec]

tmatr <- covdat[endsWith(names(covdat), "units")]
R2 <- matrix(tmatr, sqrt(length(tmatr)))



##### design matrices & observation vector for each patient #####

Zi = Xi = Yi = nXid1 = nXid2 = nXid3 = nXid4 = nZid1 = nZid2 = nZid3 = nZid4 =
  nYi1 = nYi2 = nYi3 = nYi4 <- list()

for(n in 1:length(allid)){
  
  # Z
  
  nZid1[[n]] <- nzid1 <- datfvc %>% filter(pid == allid[n]) %>% select(REcolnames) %>% as.matrix
  nZid2[[n]] <- nzid2 <- datdlco %>% filter(pid == allid[n]) %>% select(REcolnames) %>% as.matrix
  nZid3[[n]] <- nzid3 <- datef %>% filter(pid == allid[n]) %>% select(REcolnames) %>% as.matrix
  nZid4[[n]] <- nzid4 <- datrvsp %>% filter(pid == allid[n]) %>% select(REcolnames) %>% as.matrix
  
  datlist <- list(nzid1, nzid2, nzid3, nzid4)
  Zi[[n]] <- Z <- Reduce(direct.sum, datlist)
  
  
  # X
  
  nXid1[[n]] <- x1 <- datfvc %>% filter(pid == allid[n]) %>% select(FEcolnames) %>% as.matrix
  nXid2[[n]] <- x2 <- datdlco %>% filter(pid == allid[n]) %>% select(FEcolnames) %>% as.matrix
  nXid3[[n]] <- x3 <- datef %>% filter(pid == allid[n]) %>% select(FEcolnames) %>% as.matrix
  nXid4[[n]] <- x4 <- datrvsp %>% filter(pid == allid[n]) %>% select(FEcolnames) %>% as.matrix
  
  
  Xi[[n]] <-  Reduce(direct.sum, list(x1, x2, x3, x4))
  
  # Y 
  
  nYi1[[n]] <- y1 <- datfvc %>% filter(pid == allid[n]) %>% select(variablenames[1]) %>% as.matrix %>% set_colnames("y")
  nYi2[[n]] <- y2 <- datdlco %>% filter(pid == allid[n]) %>% select(variablenames[2]) %>% as.matrix %>% set_colnames("y")
  nYi3[[n]] <- y3 <- datef %>% filter(pid == allid[n]) %>% select(variablenames[3]) %>% as.matrix %>% set_colnames("y")
  nYi4[[n]] <- y4 <- datrvsp %>% filter(pid == allid[n]) %>% select(variablenames[4]) %>% as.matrix %>% set_colnames("y")
  
  Yi[[n]] <-  rbind(y1, y2, y3, y4)
  
  
}



##### Produce qqplots to check normality assumptions #####

U1 = U2 = U3 = U4 <- list()

for(n in 1:length(allid)){
  
  # data for patient n
  
  Zin <- Zi[[n]]
  Xin <- Xi[[n]] 
  Yin <- Yi[[n]]
  
  z1 <- nZid1[[n]]; z2 <- nZid2[[n]]; z3 <- nZid3[[n]]; z4 <- nZid4[[n]]
  
  datlistZ <- list(z1, z2, z3, z4)
  
  Sigma2n <- create_sigma2(datlist = datlistZ, R2 = R2)
  
  # Vinv, U
  #Vinvn <- solve(Zin %*% G2 %*% t(Zin) + Sigma2n)
  #sqrtVinv <- sqrtm(Vinvn) 
  
  V <- Zin %*% G2 %*% t(Zin) + Sigma2n
  sqrtVinv <- diag(1/sqrt(diag(V)))
  
  U <- sqrtVinv %*% (Yin - Xin %*% beta)
  
  U1[[n]] <- U[1:nrow(z1)]
  U2[[n]] <- U[(nrow(z1) + 1):(nrow(z1) + nrow(z2))]
  U3[[n]] <- U[(nrow(z1) + nrow(z2) + 1):(nrow(z1) + nrow(z2) + nrow(z3))]
  U4[[n]] <- U[(nrow(z1) + nrow(z2) + nrow(z3) + 1):length(U)]
  
}

U1vec <- unlist(U1); U2vec <- unlist(U2); U3vec <- unlist(U3); U4vec <- unlist(U4)
uvecALL <- list(U1vec, U2vec, U3vec, U4vec)

par(mfrow = c(2, 2))
qqnorm(uvecALL[[1]], pch = 20, main = "Normal Q-Q Plot for pFVC scaled residuals"); abline(0, 1, col = "darkred")
qqnorm(uvecALL[[2]], pch = 20, main = "Normal Q-Q Plot for pDLCO scaled residuals"); abline(0, 1, col = "darkred")
qqnorm(uvecALL[[3]], pch = 20, main = "Normal Q-Q Plot for EF scaled residuals"); abline(0, 1, col = "darkred")
qqnorm(uvecALL[[4]], pch = 20, main = "Normal Q-Q Plot for RVSP scaled residuals"); abline(0, 1, col = "darkred")




##### function to generate design matrix for future value & calculate mean and var #####


predictRVSPtraj <- function(G2, R2, beta, nZid1, nZid2, nZid3, nZid4,
                            nXid1, nXid2, nXid3, nXid4, nYi1, nYi2, nYi3, nYi4){
  
  dfi <- NULL
  
  for(n in 1:length(pid.cut)){
    
    z1 <- nZid1[[n]]; z2 <- nZid2[[n]]; z3 <- nZid3[[n]]; z4 <- nZid4[[n]]
    x1 <- nXid1[[n]]; x2 <- nXid2[[n]]; x3 <- nXid3[[n]]; x4 <- nXid4[[n]]
    y1 <- nYi1[[n]]; y2 <- nYi2[[n]]; y3 <- nYi3[[n]]; y4 <- nYi4[[n]]
    
    rvsptimes <- z4[, 2]
    
    dfj <- NULL
    
    for(j in 1:length(rvsptimes)){
      
      predtime <- rvsptimes[j]
      rvspobs <- y4[z4[, 2] == predtime, ]
      
      # Zi+: random effects
      Zip <- data.frame(int = 1, btime = predtime) %>%
        mutate(sp10yr = ifelse(btime + 10 > 0, btime + 10, 0), 
               sp3yr = ifelse(btime + 3 > 0, btime + 3, 0)) %>% as.matrix
      
      Zp <- Reduce(direct.sum, list(Zip, Zip, Zip, Zip))
      
      # Xip: fixed effects
      mdat <- regdat.d %>% filter(pid == pid.cut[n]) %>% head(1)
      
      
      Xip <- cbind(mdat$int, matrix(ns(mdat[, "last"] + predtime, knots = c(10, 30),
                                       Boundary.knots= c(0, 40))[1, ], nrow = 1),
                   mdat %>% select(FEcolnames[-c(1:4)])) %>% as.matrix
      
      Xp <- Reduce(direct.sum, list(Xip, Xip, Xip, Xip))
      
      
      # Z
      tm1 <- z1[z1[, 2] < predtime, ]
      tm2 <- z2[z2[, 2] < predtime, ]
      tm3 <- z3[z3[, 2] < predtime, ]
      tm4 <- z4[z4[, 2] < predtime, ]
      
      if(is.null(dim(tm1))){
        tm1 <- z1[z1[, 2] < predtime, ] %>% matrix(nrow = 1)
      }
      
      if(is.null(dim(tm2))){
        tm2 <- z2[z2[, 2] < predtime, ] %>% matrix(nrow = 1)
      }
      
      if(is.null(dim(tm3))){
        tm3 <- z3[z3[, 2] < predtime, ] %>% matrix(nrow = 1)
      }
      
      if(is.null(dim(tm4))){
        tm4 <- z4[z4[, 2] < predtime, ] %>% matrix(nrow = 1)
      }
      
      datlistZsigma = datlistZ <- list(tm1, tm2, tm3, tm4)
      
      rowlist <- sapply(datlistZ, nrow)
      
      nfvc3 <- sum(predtime - 3 <= z1[, 2] & z1[, 2] < predtime)
      ndlco3 <- sum(predtime - 3 <= z2[, 2] & z2[, 2] < predtime)
      nef3 <- sum(predtime - 3 <= z3[, 2] & z3[, 2] < predtime)
      nrvsp3 <- sum(predtime - 3 <= z4[, 2] & z4[, 2] < predtime)
      
      # no observed value for all measures
      if(all(rowlist == 0)){       
        
        # Vp
        Vp <- Zp %*% G2 %*% t(Zp) + R2
        
        # predicted values & confidence interval
        mu <- Xp %*% beta
        sigma <- Vp  
        
      }else{
        
        xm1 <- x1[1:nrow(tm1), ]
        xm2 <- x2[1:nrow(tm2), ]
        xm3 <- x3[1:nrow(tm3), ]
        xm4 <- x4[1:nrow(tm4), ]
        
        
        if(is.null(dim(xm1))){
          xm1 <- x1[1:nrow(tm1), ] %>% matrix(nrow = 1)
        }
        
        if(is.null(dim(xm2))){
          xm2 <- x2[1:nrow(tm2), ] %>% matrix(nrow = 1)
        }
        
        if(is.null(dim(xm3))){
          xm3 <- x3[1:nrow(tm3), ] %>% matrix(nrow = 1)
        }
        
        if(is.null(dim(xm4))){
          xm4 <- x4[1:nrow(tm4), ] %>% matrix(nrow = 1)
        }
        
        
        datlistX <- list(xm1, xm2, xm3, xm4)
        
        ty1 <- y1; ty2 <- y2; ty3 <- y3; ty4 <- y4
        
        # if there is no observed value for any one of the measures, remove empty dataframe
        if (any(rowlist == 0)){
          
          datlistZsigma <- datlistZ[!is.na(rowlist)]
          
          nZ <- ncol(datlistZ[[1]]); nX <- ncol(datlistX[[1]])
          
          if(rowlist[1] == 0){ty1 <- NULL
          datlistZ[[1]] <- matrix(rep(0, nZ), nrow = 1)
          datlistX[[1]] <- matrix(rep(0, nX), nrow = 1)}
          
          if(rowlist[2] == 0){ty2 <- NULL
          datlistZ[[2]] <- matrix(rep(0, nZ), nrow = 1)
          datlistX[[2]] <- matrix(rep(0, nX), nrow = 1)}
          
          if(rowlist[3] == 0){ty3 <- NULL
          datlistZ[[3]] <- matrix(rep(0, nZ), nrow = 1)
          datlistX[[3]] <- matrix(rep(0, nX), nrow = 1)}
          
          if(rowlist[4] == 0){ty4 <- NULL
          datlistZ[[4]] <- matrix(rep(0, nZ), nrow = 1)
          datlistX[[4]] <- matrix(rep(0, nX), nrow = 1)}
          
        }
        
        Zin <- Reduce(direct.sum, datlistZ) 
        Xin <- Reduce(direct.sum, datlistX)
        
        # observations from a single measure available
        if(sum(rowSums(Zin) != 0) == 1){
          
          # remove empty rows
          Zin <- Zin[rowSums(Zin) != 0, ] %>% matrix(nrow = 1)
          Xin <- Xin[rowSums(Xin) != 0, ] %>% matrix(nrow = 1)
          
        }else{
          
          # remove empty rows
          Zin <- Zin[rowSums(Zin) != 0, ]
          Xin <- Xin[rowSums(Xin) != 0, ]
          
        }
        
        Yin <- c(ty1[1:nrow(tm1)], ty2[1:nrow(tm2)], ty3[1:nrow(tm3)], ty4[1:nrow(tm4)]) %>%
          matrix(ncol = 1)
        
        Sigma2n <- create_sigma2(datlist = datlistZsigma, R2 = R2)
        
        # Vp, Cp
        Vinvn <- solve(Zin %*% G2 %*% t(Zin) + Sigma2n)
        Vp <- Zp %*% G2 %*% t(Zp) + R2
        Cp <- Zin %*% G2 %*% t(Zp)
        
        # predicted values & confidence interval
        mu <- (Xp %*% beta) + (t(Cp) %*% Vinvn %*% (Yin - Xin %*% beta))
        sigma <- Vp - (t(Cp) %*% Vinvn %*% Cp)
        
      }
      
      print(paste("n =", n, ", j =", j))
      
      newdfj <- data.frame(pid = pid.cut[n], predtime = predtime, rvspobs = rvspobs,
                           predmu = mu[4], predsigma = sigma[4, 4],
                           p45 = pnorm(q = rvsp45, mean =  mu[4], sd = sqrt(sigma[4, 4])),
                           p50 = pnorm(q = rvsp50, mean =  mu[4], sd = sqrt(sigma[4, 4])),
                           nfvc = rowlist[1], ndlco = rowlist[2], nef = rowlist[3], nrvsp = rowlist[4], 
                           nfvc3 = nfvc3, ndlco3 = ndlco3, nef3 = nef3, nrvsp3 = nrvsp3)
      
      dfj <- rbind(dfj, newdfj)
      
    }
    
    dfi <- rbind(dfi, dfj)
    
  }
  
  return(dfi)
  
}


predictFVCtraj <- function(G2, R2, beta, nZid1, nZid2, nZid3, nZid4,
                           nXid1, nXid2, nXid3, nXid4, nYi1, nYi2, nYi3, nYi4){
  
  dfi <- NULL
  
  for(n in 1:length(pid.cut)){
    
    z1 <- nZid1[[n]]; z2 <- nZid2[[n]]; z3 <- nZid3[[n]]; z4 <- nZid4[[n]]
    x1 <- nXid1[[n]]; x2 <- nXid2[[n]]; x3 <- nXid3[[n]]; x4 <- nXid4[[n]]
    y1 <- nYi1[[n]]; y2 <- nYi2[[n]]; y3 <- nYi3[[n]]; y4 <- nYi4[[n]]
    
    fvctimes <- z1[, 2]
    
    dfj <- NULL
    
    for(j in 1:length(fvctimes)){
      
      predtime <- fvctimes[j]
      fvcobs <- y1[z1[, 2] == predtime, ]
      
      # Zi+: random effects
      Zip <- data.frame(int = 1, btime = predtime) %>%
        mutate(sp10yr = ifelse(btime + 10 > 0, btime + 10, 0), 
               sp3yr = ifelse(btime + 3 > 0, btime + 3, 0)) %>% as.matrix
      
      Zp <- Reduce(direct.sum, list(Zip, Zip, Zip, Zip))
      
      # Xip: fixed effects
      mdat <- regdat.d %>% filter(pid == pid.cut[n]) %>% head(1)
      
      
      Xip <- cbind(mdat$int, matrix(ns(mdat[, "last"] + predtime, knots = c(10, 30),
                                       Boundary.knots= c(0, 40))[1, ], nrow = 1),
                   mdat %>% select(FEcolnames[-c(1:4)])) %>% as.matrix
      
      Xp <- Reduce(direct.sum, list(Xip, Xip, Xip, Xip))
      
      
      # Z
      tm1 <- z1[z1[, 2] < predtime, ]
      tm2 <- z2[z2[, 2] < predtime, ]
      tm3 <- z3[z3[, 2] < predtime, ]
      tm4 <- z4[z4[, 2] < predtime, ]
      
      if(is.null(dim(tm1))){
        tm1 <- z1[z1[, 2] < predtime, ] %>% matrix(nrow = 1)
      }
      
      if(is.null(dim(tm2))){
        tm2 <- z2[z2[, 2] < predtime, ] %>% matrix(nrow = 1)
      }
      
      if(is.null(dim(tm3))){
        tm3 <- z3[z3[, 2] < predtime, ] %>% matrix(nrow = 1)
      }
      
      if(is.null(dim(tm4))){
        tm4 <- z4[z4[, 2] < predtime, ] %>% matrix(nrow = 1)
      }
      
      datlistZsigma = datlistZ <- list(tm1, tm2, tm3, tm4)
      
      rowlist <- sapply(datlistZ, nrow)
      
      nfvc3 <- sum(predtime - 3 <= z1[, 2] & z1[, 2] < predtime)
      ndlco3 <- sum(predtime - 3 <= z2[, 2] & z2[, 2] < predtime)
      nef3 <- sum(predtime - 3 <= z3[, 2] & z3[, 2] < predtime)
      nrvsp3 <- sum(predtime - 3 <= z4[, 2] & z4[, 2] < predtime)
      
      # no observed value for all measures
      if(all(rowlist == 0)){       
        
        # Vp
        Vp <- Zp %*% G2 %*% t(Zp) + R2
        
        # predicted values & confidence interval
        mu <- Xp %*% beta
        sigma <- Vp  
        
      }else{
        
        
        xm1 <- x1[1:nrow(tm1), ]
        xm2 <- x2[1:nrow(tm2), ]
        xm3 <- x3[1:nrow(tm3), ]
        xm4 <- x4[1:nrow(tm4), ]
        
        
        if(is.null(dim(xm1))){
          xm1 <- x1[1:nrow(tm1), ] %>% matrix(nrow = 1)
        }
        
        if(is.null(dim(xm2))){
          xm2 <- x2[1:nrow(tm2), ] %>% matrix(nrow = 1)
        }
        
        if(is.null(dim(xm3))){
          xm3 <- x3[1:nrow(tm3), ] %>% matrix(nrow = 1)
        }
        
        if(is.null(dim(xm4))){
          xm4 <- x4[1:nrow(tm4), ] %>% matrix(nrow = 1)
        }
        
        datlistX <- list(xm1, xm2, xm3, xm4)
        
        ty1 <- y1; ty2 <- y2; ty3 <- y3; ty4 <- y4
        
        # if there is no observed value for any one of the measures, remove empty dataframe
        if (any(rowlist == 0)){
          
          datlistZsigma <- datlistZ[!is.na(rowlist)]
          
          nZ <- ncol(datlistZ[[1]]); nX <- ncol(datlistX[[1]])
          
          if(rowlist[1] == 0){ty1 <- NULL
          datlistZ[[1]] <- matrix(rep(0, nZ), nrow = 1)
          datlistX[[1]] <- matrix(rep(0, nX), nrow = 1)}
          
          if(rowlist[2] == 0){ty2 <- NULL
          datlistZ[[2]] <- matrix(rep(0, nZ), nrow = 1)
          datlistX[[2]] <- matrix(rep(0, nX), nrow = 1)}
          
          if(rowlist[3] == 0){ty3 <- NULL
          datlistZ[[3]] <- matrix(rep(0, nZ), nrow = 1)
          datlistX[[3]] <- matrix(rep(0, nX), nrow = 1)}
          
          if(rowlist[4] == 0){ty4 <- NULL
          datlistZ[[4]] <- matrix(rep(0, nZ), nrow = 1)
          datlistX[[4]] <- matrix(rep(0, nX), nrow = 1)}
          
        }
        
        Zin <- Reduce(direct.sum, datlistZ) 
        Xin <- Reduce(direct.sum, datlistX)
        
        # observations from a single measure available
        if(sum(rowSums(Zin) != 0) == 1){
          
          # remove empty rows
          Zin <- Zin[rowSums(Zin) != 0, ] %>% matrix(nrow = 1)
          Xin <- Xin[rowSums(Xin) != 0, ] %>% matrix(nrow = 1)
          
        }else{
          
          # remove empty rows
          Zin <- Zin[rowSums(Zin) != 0, ]
          Xin <- Xin[rowSums(Xin) != 0, ]
          
        }
        
        Yin <- c(ty1[1:nrow(tm1)], ty2[1:nrow(tm2)], ty3[1:nrow(tm3)], ty4[1:nrow(tm4)]) %>% matrix(ncol = 1)
        
        Sigma2n <- create_sigma2(datlist = datlistZsigma, R2 = R2)
        
        # Vp, Cp
        Vinvn <- solve(Zin %*% G2 %*% t(Zin) + Sigma2n)
        Vp <- Zp %*% G2 %*% t(Zp) + R2
        Cp <- Zin %*% G2 %*% t(Zp)
        
        # predicted values & confidence interval
        mu <- (Xp %*% beta) + (t(Cp) %*% Vinvn %*% (Yin - Xin %*% beta))
        sigma <- Vp - (t(Cp) %*% Vinvn %*% Cp)
        
      }
      
      print(paste("n =", n, ", j =", j))
      
      
      newdfj <- data.frame(pid = pid.cut[n], predtime = predtime, fvcobs = fvcobs,
                           predmu = mu[1], predsigma = sigma[1, 1],
                           p60 = pnorm(q = fvc60, mean =  mu[1], sd = sqrt(sigma[1, 1])),
                           p70 = pnorm(q = fvc70, mean =  mu[1], sd = sqrt(sigma[1, 1])),
                           nfvc = rowlist[1], ndlco = rowlist[2], nef = rowlist[3], nrvsp = rowlist[4],
                           nfvc3 = nfvc3, ndlco3 = ndlco3, nef3 = nef3, nrvsp3 = nrvsp3)
      
      dfj <- rbind(dfj, newdfj)
      
    }
    
    dfi <- rbind(dfi, dfj)
    
  }
  
  return(dfi)
  
}


predictEFtraj <- function(G2, R2, beta, nZid1, nZid2, nZid3, nZid4,
                          nXid1, nXid2, nXid3, nXid4, nYi1, nYi2, nYi3, nYi4){
  
  dfi <- NULL
  
  for(n in 1:length(pid.cut)){
    
    z1 <- nZid1[[n]]; z2 <- nZid2[[n]]; z3 <- nZid3[[n]]; z4 <- nZid4[[n]]
    x1 <- nXid1[[n]]; x2 <- nXid2[[n]]; x3 <- nXid3[[n]]; x4 <- nXid4[[n]]
    y1 <- nYi1[[n]]; y2 <- nYi2[[n]]; y3 <- nYi3[[n]]; y4 <- nYi4[[n]]
    
    eftimes <- z3[, 2]
    
    dfj <- NULL
    
    for(j in 1:length(eftimes)){
      
      predtime <- eftimes[j]
      efobs <- y3[z3[, 2] == predtime, ]
      
      # Zi+: random effects
      Zip <- data.frame(int = 1, btime = predtime) %>%
        mutate(sp10yr = ifelse(btime + 10 > 0, btime + 10, 0), 
               sp3yr = ifelse(btime + 3 > 0, btime + 3, 0)) %>% as.matrix
      
      Zp <- Reduce(direct.sum, list(Zip, Zip, Zip, Zip))
      
      # Xip: fixed effects
      mdat <- regdat.d %>% filter(pid == pid.cut[n]) %>% head(1)
      
      
      Xip <- cbind(mdat$int, matrix(ns(mdat[, "last"] + predtime, knots = c(10, 30),
                                       Boundary.knots= c(0, 40))[1, ], nrow = 1),
                   mdat %>% select(FEcolnames[-c(1:4)])) %>% as.matrix
      
      Xp <- Reduce(direct.sum, list(Xip, Xip, Xip, Xip))
      
      
      # Z
      tm1 <- z1[z1[, 2] < predtime, ]
      tm2 <- z2[z2[, 2] < predtime, ]
      tm3 <- z3[z3[, 2] < predtime, ]
      tm4 <- z4[z4[, 2] < predtime, ]
      
      if(is.null(dim(tm1))){
        tm1 <- z1[z1[, 2] < predtime, ] %>% matrix(nrow = 1)
      }
      
      if(is.null(dim(tm2))){
        tm2 <- z2[z2[, 2] < predtime, ] %>% matrix(nrow = 1)
      }
      
      if(is.null(dim(tm3))){
        tm3 <- z3[z3[, 2] < predtime, ] %>% matrix(nrow = 1)
      }
      
      if(is.null(dim(tm4))){
        tm4 <- z4[z4[, 2] < predtime, ] %>% matrix(nrow = 1)
      }
      
      
      datlistZsigma = datlistZ <- list(tm1, tm2, tm3, tm4)
      
      rowlist <- sapply(datlistZ, nrow)
      
      nfvc3 <- sum(predtime - 3 <= z1[, 2] & z1[, 2] < predtime)
      ndlco3 <- sum(predtime - 3 <= z2[, 2] & z2[, 2] < predtime)
      nef3 <- sum(predtime - 3 <= z3[, 2] & z3[, 2] < predtime)
      nrvsp3 <- sum(predtime - 3 <= z4[, 2] & z4[, 2] < predtime)
      
      # no observed value for all measures
      if(all(rowlist == 0)){       
        
        # Vp
        Vp <- Zp %*% G2 %*% t(Zp) + R2
        
        # predicted values & confidence interval
        mu <- Xp %*% beta
        sigma <- Vp  
        
      }else{
        
        
        xm1 <- x1[1:nrow(tm1), ]
        xm2 <- x2[1:nrow(tm2), ]
        xm3 <- x3[1:nrow(tm3), ]
        xm4 <- x4[1:nrow(tm4), ]
        
        
        if(is.null(dim(xm1))){
          xm1 <- x1[1:nrow(tm1), ] %>% matrix(nrow = 1)
        }
        
        if(is.null(dim(xm2))){
          xm2 <- x2[1:nrow(tm2), ] %>% matrix(nrow = 1)
        }
        
        if(is.null(dim(xm3))){
          xm3 <- x3[1:nrow(tm3), ] %>% matrix(nrow = 1)
        }
        
        if(is.null(dim(xm4))){
          xm4 <- x4[1:nrow(tm4), ] %>% matrix(nrow = 1)
        }
        
        datlistX <- list(xm1, xm2, xm3, xm4)
        
        ty1 <- y1; ty2 <- y2; ty3 <- y3; ty4 <- y4
        
        # if there is no observed value for any one of the measures, remove empty dataframe
        if (any(rowlist == 0)){
          
          datlistZsigma <- datlistZ[!is.na(rowlist)]
          
          nZ <- ncol(datlistZ[[1]]); nX <- ncol(datlistX[[1]])
          
          if(rowlist[1] == 0){ty1 <- NULL
          datlistZ[[1]] <- matrix(rep(0, nZ), nrow = 1)
          datlistX[[1]] <- matrix(rep(0, nX), nrow = 1)}
          
          if(rowlist[2] == 0){ty2 <- NULL
          datlistZ[[2]] <- matrix(rep(0, nZ), nrow = 1)
          datlistX[[2]] <- matrix(rep(0, nX), nrow = 1)}
          
          if(rowlist[3] == 0){ty3 <- NULL
          datlistZ[[3]] <- matrix(rep(0, nZ), nrow = 1)
          datlistX[[3]] <- matrix(rep(0, nX), nrow = 1)}
          
          if(rowlist[4] == 0){ty4 <- NULL
          datlistZ[[4]] <- matrix(rep(0, nZ), nrow = 1)
          datlistX[[4]] <- matrix(rep(0, nX), nrow = 1)}
          
        }
        
        Zin <- Reduce(direct.sum, datlistZ) 
        Xin <- Reduce(direct.sum, datlistX)
        
        # observations from a single measure available
        if(sum(rowSums(Zin) != 0) == 1){
          
          # remove empty rows
          Zin <- Zin[rowSums(Zin) != 0, ] %>% matrix(nrow = 1)
          Xin <- Xin[rowSums(Xin) != 0, ] %>% matrix(nrow = 1)
          
        }else{
          
          # remove empty rows
          Zin <- Zin[rowSums(Zin) != 0, ]
          Xin <- Xin[rowSums(Xin) != 0, ]
          
        }
        
        Yin <- c(ty1[1:nrow(tm1)], ty2[1:nrow(tm2)], ty3[1:nrow(tm3)], ty4[1:nrow(tm4)]) %>% matrix(ncol = 1)
        
        Sigma2n <- create_sigma2(datlist = datlistZsigma, R2 = R2)
        
        # Vp, Cp
        Vinvn <- solve(Zin %*% G2 %*% t(Zin) + Sigma2n)
        Vp <- Zp %*% G2 %*% t(Zp) + R2
        Cp <- Zin %*% G2 %*% t(Zp)
        
        # predicted values & confidence interval
        mu <- (Xp %*% beta) + (t(Cp) %*% Vinvn %*% (Yin - Xin %*% beta))
        sigma <- Vp - (t(Cp) %*% Vinvn %*% Cp)
        
      }
      
      print(paste("n =", n, ", j =", j))
      
      
      newdfj <- data.frame(pid = pid.cut[n], predtime = predtime, efobs = efobs,
                           predmu = mu[3], predsigma = sigma[3, 3],
                           p35 = pnorm(q = ef35, mean =  mu[3], sd = sqrt(sigma[3, 3])),
                           p50 = pnorm(q = ef50, mean =  mu[3], sd = sqrt(sigma[3, 3])),
                           nfvc = rowlist[1], ndlco = rowlist[2], nef = rowlist[3], nrvsp = rowlist[4],
                           nfvc3 = nfvc3, ndlco3 = ndlco3, nef3 = nef3, nrvsp3 = nrvsp3)
      
      dfj <- rbind(dfj, newdfj)
      
    }
    
    dfi <- rbind(dfi, dfj)
    
  }
  
  return(dfi)
  
}



##### fit model & compute probability for all folds #####

# cutoffs generated from the cohort data

# RVSP cutoffs .

rvsp50 <- -1.33; rvsp45 <- -1.09

# FVC cutoffs
fvc70 <- -0.41; fvc60 <- -0.91

# EF cutoffs 
ef50 <- -1.75; ef35 <- -2.43


cvdf.rvsp = cvdf.fvc = cvdf.ef <- NULL

for (k in 1:K){
  
  # patient ids in k folds
  
  pid.cut <- pid.cv[which(folds == k)]
  pid.rest <-  pid.cv[which(folds != k)]
  
  pdat.rest <- regdat.d %>% filter(pid %in% pid.rest)
  
  
  # fit model
  
  if (!file.exists(paste0("rvspeffvcdlco", k,".Rdata"))){
    
    combmcmc <- MCMCglmm(as.formula(paste("cbind(", paste(variablenames, collapse = ","),
                                          ") ~ -1 + trait + trait:(ns(time, knots = c(10, 30), Boundary.knots= c(0, 40)) +
                                      AgeOnset + race_AA + Type_Diffuse + Type_Sine + Sex_Male + ACApos + RNAPolpos + SCL70pos)")),
                         random = ~ us(trait + trait:btime + trait:sp10yr + trait:sp3yr):pid, pr = TRUE,
                         rcov = ~ us(trait):units,
                         prior = list(R = list(V = V.R, nu = nu.R),
                                      G = list(G1 = list(V = V, nu = nu))),
                         family = rep("gaussian", nmeas),
                         nitt = nitt, burnin = burnin, data = pdat.rest)
    
    save(combmcmc, file = paste0("rvspeffvcdlco", k,".Rdata"))
    
  }else{
    
    load(paste0("rvspeffvcdlco", k,".Rdata"))
    
  }
  
  # Varcov matrix from mcmcglmm: R and G
  
  covdat <- apply(combmcmc$VCV, 2, mean)
  
  tmat <- covdat[endsWith(names(covdat), "pid")]
  G <- matrix(tmat, sqrt(length(tmat)))
  G2 <- G[ordervec, ordervec]
  
  tmatr <- covdat[endsWith(names(covdat), "units")]
  R2 <- matrix(tmatr, sqrt(length(tmatr)))
  
  
  # Fixed effect etimates
  
  solmat <- combmcmc$Sol
  sol.postmean <- apply(solmat[, 1:nfe], 2, mean)
  
  cnames <- names(sol.postmean)
  colorder <- c(cnames[startsWith(cnames, "traitpFVCt")], 
                cnames[startsWith(cnames, "traitpDLCOt")],
                cnames[startsWith(cnames, "traitEFt")],
                cnames[startsWith(cnames, "traitRVSPt")])
  
  beta <- sol.postmean[order(match(names(sol.postmean), colorder))] 
  
  pidind <- match(pid.cut, allid)
  
  newcvdf.rvsp <- predictRVSPtraj(G2 = G2, R2 = R2, beta = beta, 
                                  nZid1 = nZid1[pidind], nZid2 = nZid2[pidind], nZid3 = nZid3[pidind], nZid4 = nZid4[pidind], 
                                  nXid1 = nXid1[pidind], nXid2 = nXid2[pidind], nXid3 = nXid3[pidind], nXid4 = nXid4[pidind], 
                                  nYi1 = nYi1[pidind], nYi2 = nYi2[pidind], nYi3 = nYi3[pidind], nYi4 = nYi4[pidind])
  
  newcvdf.fvc <- predictFVCtraj(G2 = G2, R2 = R2, beta = beta, 
                                nZid1 = nZid1[pidind], nZid2 = nZid2[pidind], nZid3 = nZid3[pidind], nZid4 = nZid4[pidind], 
                                nXid1 = nXid1[pidind], nXid2 = nXid2[pidind], nXid3 = nXid3[pidind], nXid4 = nXid4[pidind], 
                                nYi1 = nYi1[pidind], nYi2 = nYi2[pidind], nYi3 = nYi3[pidind], nYi4 = nYi4[pidind])
  
  newcvdf.ef <- predictEFtraj(G2 = G2, R2 = R2, beta = beta, 
                              nZid1 = nZid1[pidind], nZid2 = nZid2[pidind], nZid3 = nZid3[pidind], nZid4 = nZid4[pidind], 
                              nXid1 = nXid1[pidind], nXid2 = nXid2[pidind], nXid3 = nXid3[pidind], nXid4 = nXid4[pidind], 
                              nYi1 = nYi1[pidind], nYi2 = nYi2[pidind], nYi3 = nYi3[pidind], nYi4 = nYi4[pidind])
  
  cvdf.rvsp <- rbind(cvdf.rvsp, newcvdf.rvsp)
  cvdf.fvc <- rbind(cvdf.fvc, newcvdf.fvc)
  cvdf.ef <- rbind(cvdf.ef, newcvdf.ef)
  
}


# estimated probabilities
cvdf.rvsp <- cvdf.rvsp %>% mutate(rvsp50obs = ifelse(rvspobs <= rvsp50, 1, 0), rvsp45obs = ifelse(rvspobs <= rvsp45, 1, 0))
cvdf.fvc <- cvdf.fvc %>% mutate(fvc70obs = ifelse(fvcobs <= fvc70, 1, 0), fvc60obs = ifelse(fvcobs <= fvc60, 1, 0))
cvdf.ef <- cvdf.ef %>% mutate(ef50obs = ifelse(efobs < ef50, 1, 0), ef35obs = ifelse(efobs < ef35, 1, 0))

#save(cvdf.rvsp, file = "cvdfrvsp.Rdata")
#save(cvdf.fvc, file = "cvdffvc.Rdata")
#save(cvdf.ef, file = "cvdfef.Rdata")

# rvsp
pROC_50 <- roc(cvdf.rvsp$rvsp50obs, cvdf.rvsp$p50,
               smoothed = TRUE,
               # arguments for ci
               ci = TRUE, ci.alpha = 0.95, stratified = FALSE,
               # arguments for plot
               plot = TRUE, auc.polygon = T, max.auc.polygon = T, grid=F,
               print.auc = TRUE, show.thres=TRUE, main = "RVSP > 50")


pROC_45 <- roc(cvdf.rvsp$rvsp45obs, cvdf.rvsp$p45,
               smoothed = TRUE,
               # arguments for ci
               ci = TRUE, ci.alpha = 0.95, stratified = FALSE,
               # arguments for plot
               plot = TRUE, auc.polygon = T, max.auc.polygon = T, grid=F,
               print.auc = TRUE, show.thres=TRUE, main = "RVSP > 45")


# FVC
pROC_60 <- roc(cvdf.fvc$fvc60obs, cvdf.fvc$p60,
               smoothed = TRUE,
               # arguments for ci
               ci = TRUE, ci.alpha = 0.95, stratified = FALSE,
               # arguments for plot
               plot = TRUE, auc.polygon = T, max.auc.polygon = T, grid=F,
               print.auc = TRUE, show.thres=TRUE, main = "FVC <= 60")


pROC_70 <- roc(cvdf.fvc$fvc70obs, cvdf.fvc$p70,
               smoothed = TRUE,
               # arguments for ci
               ci = TRUE, ci.alpha = 0.95, stratified = FALSE,
               # arguments for plot
               plot = TRUE, auc.polygon = T, max.auc.polygon = T, grid=F,
               print.auc = TRUE, show.thres=TRUE, main = "FVC <= 70")


# EF
pROC_50 <- roc(cvdf.ef$ef50obs, cvdf.ef$p50,
               smoothed = TRUE,
               # arguments for ci
               ci = TRUE, ci.alpha = 0.95, stratified = FALSE,
               # arguments for plot
               plot = TRUE, auc.polygon = T, max.auc.polygon = T, grid=F,
               print.auc = TRUE, show.thres=TRUE, main = "EF < 50")


pROC_35 <- roc(cvdf.ef$ef35obs, cvdf.ef$p35,
               smoothed = TRUE,
               # arguments for ci
               ci = TRUE, ci.alpha = 0.95, stratified = FALSE,
               # arguments for plot
               plot = TRUE, auc.polygon = T, max.auc.polygon = T, grid=F,
               print.auc = TRUE, show.thres=TRUE , main = "EF < 35")


##### CV AUC by the number of observation #####

rvsp3yr0 <- cvdf.rvsp %>% filter(nrvsp3 == 0)
rvsp3yr1 <- cvdf.rvsp %>% filter(nrvsp3 == 1)
rvsp3yr2 <- cvdf.rvsp %>% filter(nrvsp3 == 2)
rvsp3yr3up <- cvdf.rvsp %>% filter(nrvsp3 > 2)

par(mfrow = c(2, 2))

pROC_rvspn0 <- roc(rvsp3yr0$rvsp50obs, rvsp3yr0$p50,
                   smoothed = TRUE,
                   # arguments for ci
                   ci = TRUE, ci.alpha = 0.95, stratified = FALSE,
                   # arguments for plot
                   plot = TRUE, auc.polygon = T, max.auc.polygon = T, grid=F,
                   print.auc = TRUE, show.thres=TRUE, main = "rvsp nobs = 0")

pROC_rvspn1 <- roc(rvsp3yr1$rvsp50obs, rvsp3yr1$p50,
                   smoothed = TRUE,
                   # arguments for ci
                   ci = TRUE, ci.alpha = 0.95, stratified = FALSE,
                   # arguments for plot
                   plot = TRUE, auc.polygon = T, max.auc.polygon = T, grid=F,
                   print.auc = TRUE, show.thres=TRUE, main = "rvsp nobs = 1")

pROC_rvspn2 <- roc(rvsp3yr2$rvsp50obs, rvsp3yr2$p50,
                   smoothed = TRUE,
                   # arguments for ci
                   ci = TRUE, ci.alpha = 0.95, stratified = FALSE,
                   # arguments for plot
                   plot = TRUE, auc.polygon = T, max.auc.polygon = T, grid=F,
                   print.auc = TRUE, show.thres=TRUE, main = "rvsp nobs = 2")


pROC_rvspn3up <- roc(rvsp3yr3up$rvsp50obs, rvsp3yr3up$p50,
                     smoothed = TRUE,
                     # arguments for ci
                     ci = TRUE, ci.alpha = 0.95, stratified = FALSE,
                     # arguments for plot
                     plot = TRUE, auc.polygon = T, max.auc.polygon = T, grid=F,
                     print.auc = TRUE, show.thres=TRUE, main = "rvsp nobs > 2")



par(mfrow = c(2, 2))

pROC_rvspn0 <- roc(rvsp3yr0$rvsp45obs, rvsp3yr0$p45,
                   smoothed = TRUE,
                   # arguments for ci
                   ci = TRUE, ci.alpha = 0.95, stratified = FALSE,
                   # arguments for plot
                   plot = TRUE, auc.polygon = T, max.auc.polygon = T, grid=F,
                   print.auc = TRUE, show.thres=TRUE, main = "rvsp nobs = 0")

pROC_rvspn1 <- roc(rvsp3yr1$rvsp45obs, rvsp3yr1$p45,
                   smoothed = TRUE,
                   # arguments for ci
                   ci = TRUE, ci.alpha = 0.95, stratified = FALSE,
                   # arguments for plot
                   plot = TRUE, auc.polygon = T, max.auc.polygon = T, grid=F,
                   print.auc = TRUE, show.thres=TRUE, main = "rvsp nobs = 1")

pROC_rvspn2 <- roc(rvsp3yr2$rvsp45obs, rvsp3yr2$p45,
                   smoothed = TRUE,
                   # arguments for ci
                   ci = TRUE, ci.alpha = 0.95, stratified = FALSE,
                   # arguments for plot
                   plot = TRUE, auc.polygon = T, max.auc.polygon = T, grid=F,
                   print.auc = TRUE, show.thres=TRUE, main = "rvsp nobs = 2")

pROC_rvspn3up <- roc(rvsp3yr3up$rvsp45obs, rvsp3yr3up$p45,
                     smoothed = TRUE,
                     # arguments for ci
                     ci = TRUE, ci.alpha = 0.95, stratified = FALSE,
                     # arguments for plot
                     plot = TRUE, auc.polygon = T, max.auc.polygon = T, grid=F,
                     print.auc = TRUE, show.thres=TRUE, main = "rvsp nobs > 2")




ef3yr0 <- cvdf.ef %>% filter(nef3 == 0)
ef3yr1 <- cvdf.ef %>% filter(nef3 == 1)
ef3yr2 <- cvdf.ef %>% filter(nef3 == 2)
ef3yr3up <- cvdf.ef %>% filter(nef3 > 2)

par(mfrow = c(2, 2))

pROC_efn0 <- roc(ef3yr0$ef50obs, ef3yr0$p50,
                 smoothed = TRUE,
                 # arguments for ci
                 ci = TRUE, ci.alpha = 0.95, stratified = FALSE,
                 # arguments for plot
                 plot = TRUE, auc.polygon = T, max.auc.polygon = T, grid=F,
                 print.auc = TRUE, show.thres=TRUE, main = expression('n'[EFi]*' = 0'))

pROC_efn1 <- roc(ef3yr1$ef50obs, ef3yr1$p50,
                 smoothed = TRUE,
                 # arguments for ci
                 ci = TRUE, ci.alpha = 0.95, stratified = FALSE,
                 # arguments for plot
                 plot = TRUE, auc.polygon = T, max.auc.polygon = T, grid=F,
                 print.auc = TRUE, show.thres=TRUE, main = expression('n'[EFi]*' = 1'))

pROC_efn2 <- roc(ef3yr2$ef50obs, ef3yr2$p50,
                 smoothed = TRUE,
                 # arguments for ci
                 ci = TRUE, ci.alpha = 0.95, stratified = FALSE,
                 # arguments for plot
                 plot = TRUE, auc.polygon = T, max.auc.polygon = T, grid=F,
                 print.auc = TRUE, show.thres=TRUE, main = expression('n'[EFi]*' = 2'))


pROC_efn3up <- roc(ef3yr3up$ef50obs, ef3yr3up$p50,
                   smoothed = TRUE,
                   # arguments for ci
                   ci = TRUE, ci.alpha = 0.95, stratified = FALSE,
                   # arguments for plot
                   plot = TRUE, auc.polygon = T, max.auc.polygon = T, grid=F,
                   print.auc = TRUE, show.thres=TRUE, main = expression('n'[EFi]*' > 2'))



par(mfrow = c(2, 2))

pROC_efn0 <- roc(ef3yr0$ef35obs, ef3yr0$p35,
                 smoothed = TRUE,
                 # arguments for ci
                 ci = TRUE, ci.alpha = 0.95, stratified = FALSE,
                 # arguments for plot
                 plot = TRUE, auc.polygon = T, max.auc.polygon = T, grid=F,
                 print.auc = TRUE, show.thres=TRUE, main = "EF nobs = 0")

pROC_efn1 <- roc(ef3yr1$ef35obs, ef3yr1$p35,
                 smoothed = TRUE,
                 # arguments for ci
                 ci = TRUE, ci.alpha = 0.95, stratified = FALSE,
                 # arguments for plot
                 plot = TRUE, auc.polygon = T, max.auc.polygon = T, grid=F,
                 print.auc = TRUE, show.thres=TRUE, main = "EF nobs = 1")

pROC_efn2 <- roc(ef3yr2$ef35obs, ef3yr2$p35,
                 smoothed = TRUE,
                 # arguments for ci
                 ci = TRUE, ci.alpha = 0.95, stratified = FALSE,
                 # arguments for plot
                 plot = TRUE, auc.polygon = T, max.auc.polygon = T, grid=F,
                 print.auc = TRUE, show.thres=TRUE, main = "EF nobs = 2")

pROC_efn3up <- roc(ef3yr3up$ef35obs, ef3yr3up$p35,
                   smoothed = TRUE,
                   # arguments for ci
                   ci = TRUE, ci.alpha = 0.95, stratified = FALSE,
                   # arguments for plot
                   plot = TRUE, auc.polygon = T, max.auc.polygon = T, grid=F,
                   print.auc = TRUE, show.thres=TRUE, main = "EF nobs > 2")




