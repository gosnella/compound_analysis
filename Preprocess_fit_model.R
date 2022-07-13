IMODEL <- __IMODEL__

library(compounds)

############################# Read Data ##############################

CPD <- read.csv("data/CPD_FP_PRJ01_UoB.csv", stringsAsFactors = FALSE)
CPD$PRJID <- NULL
Fingerprints <- CPD[, 2:ncol(CPD)] > 0
CPDID <- CPD$CPDID

SpeciesData <- read.csv("data/Species_UoB.csv", stringsAsFactors = FALSE)

EPS <- read.csv("data/EPS_PRJ01_UoB.csv", stringsAsFactors = FALSE)
EPS$CPDID <- factor(EPS$CPDID, levels = CPDID)
EPSDATA <- NULL
icol <- c(1, 3, 4) # Keep these columns
for (i in grep("^POST_SPC[0-9]+", names(EPS))) {
  Damage <- as.integer(EPS[, i])
  Species <- as.integer(sub("^POST_SPC", "", names(EPS)[i]))
  subdata <- c(EPS[, icol],
               list(Stage = "POST", Species = Species,
                    Damage = Damage),
               SpeciesData[Species, 2:4])
  subdata <- data.frame(subdata, stringsAsFactors = FALSE)
  subdata <- subdata[!is.na(Damage), ]
  EPSDATA <- rbind(EPSDATA, subdata)
}
for (i in grep("^PRE_SPC[0-9]+", names(EPS))) {
  Damage <- as.integer(EPS[, i])
  Species <- as.integer(sub("^PRE_SPC", "", names(EPS)[i]))
  subdata <- c(EPS[, icol],
               list(Stage = "PRE", Species = Species,
                    Damage = Damage),
               SpeciesData[Species, 2:4])
  subdata <- data.frame(subdata, stringsAsFactors = FALSE)
  subdata <- subdata[!is.na(Damage), ]
  EPSDATA <- rbind(EPSDATA, subdata)
}
EPSDATA$Stage <- factor(EPSDATA$Stage, levels = c("PRE", "POST"))
EPSDATA$Species <- factor(EPSDATA$Species, levels = SpeciesData$Species)
EPSDATA$TransfDam <- qlogis(pmin(95, pmax(5, EPSDATA$Damage))/100)
ToxicityBreaks <- c(20L, 50L, 80L)
EPSDATA$ToxicityLevel <- cut(EPSDATA$Damage, c(-10L, ToxicityBreaks, 110L), labels = 0:length(ToxicityBreaks), ordered_result = TRUE)

############################# Fit model ##############################

flink_fit <- c("logit", "probit", "GEVmod", "GEVmodNS")
fcor_fit <- c("tanimoto", "exponential", "gaussian", "independent")
models_fit <- expand.grid(flink = flink_fit, fcor = fcor_fit,
                          stringsAsFactors = FALSE)

nData <- nrow(EPSDATA)
nID <- table(EPSDATA$CPDID)
lFit <- round(nID*.8)
nFit <- sum(lFit)
set.seed(150)
iFit <- sort(do.call(c, lapply(seq_along(CPDID),
               function(i) sample(which(EPSDATA$CPDID == CPDID[i]), lFit[i]))))
iPrd <- (1:nData)[-iFit]

Fixed <- ToxicityLevel ~ log(Rate/1000) + Monocot + Warm + Stage
Random <- ~ CPDID

fit <- fit_ordinal_compound(Fixed, Random,
                            Fingerprints, EPSDATA,
                            subset = iFit, 
                            flink = models_fit$flink[IMODEL], 
                            fcor = models_fit$fcor[IMODEL],
                            returnData = FALSE)
prd <- predict_compound(fit, EPSDATA[iPrd, ])
slg <- score_model_compound(EPSDATA$ToxicityLevel[iPrd], prd$Prediction, "log")
sbr <- score_model_compound(EPSDATA$ToxicityLevel[iPrd], prd$Prediction, "Bri")
ssp <- score_model_compound(EPSDATA$ToxicityLevel[iPrd], prd$Prediction, "sph")

RDATAFILE <- Sys.getenv("RDATAFILE")
RFILENAME <- Sys.getenv("RFILENAME")
if (RDATAFILE == "") {
  if (RFILENAME == "") RFILENAME <- tempfile("run_", ".")
  RDATAFILE <- paste0(RFILENAME,
                      "_",format(Sys.Date(),"%y%m%d"),".RData")
}
save(IMODEL, fit, prd, slg, sbr, ssp, file = RDATAFILE)
