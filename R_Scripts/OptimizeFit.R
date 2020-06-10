
require(Matrix) 	
require(data.table) 
require(deSolve)

if (!exists("OutputStringPrefix")) OutputStringPrefix = ""

# Extract params from ParamVec
Tau         		= ParamVec$Tau        
SwitchFreq			= ParamVec$SwitchFreq	
Pi          		= ParamVec$Pi         
delta				= ParamVec$delta		
H           		= ParamVec$H     
K           		= ParamVec$K     
Duration			= ParamVec$Duration     
TimeStep 			= ParamVec$TimeStep		
Equilib_Diversity 	= ParamVec$Equilib_Diversity		

Timepoints 				= seq(0, Duration + TimeStep, by = TimeStep)
MaxStochasticDuration	= max(Timepoints) + 1 

#### import CloneDataFrame_Initial and PromotionDataFrame_Initial data frames
CloneDataFrame_Initial 		= as.data.table(read.table(file = file.path(OutputDir, paste0(OutputStringPrefix, "CloneLookUpTable_", dta, ".txt"))	, header = T, sep = "\t"))
PromotionDataFrame_Initial 	= read.table(file = file.path(OutputDir, paste0(OutputStringPrefix, "PromotionDFrame_", dta, ".txt"))	, header = T, sep = "\t")

#### import functions for hybrid
source(file.path(R_Script_Dir, "HybridFunctions.R"))

# import grid search results
ParamGridFileName	= file.path(OutputDir, paste0(OutputStringPrefix, "ParamGridFit_", dta, ".txt"))
ResultsDFrame 		= read.table(file = ParamGridFileName, header = T, sep = "\t")

#### ==== #### ==== #### ==== #### ==== #### ==== #### ==== #### ==== #### ==== #### ==== #### ==== #### ==== #### ==== #### ==== #### ==== #### ==== #### ==== 
#### ==== #### ==== #### ==== #### ==== #### ==== #### ==== #### ==== #### ==== #### ==== #### ==== #### ==== #### ==== #### ==== #### ==== #### ==== #### ==== 

#### Various pre-calculated quantities to reduce runtime.  
MaxPlausibleNoDeterministicBins = 15  ### helps for calculating WhichClonesAreDeterministic quickly - no need to search entire data frame. 
NoRows_CloneDF 					= dim(CloneDataFrame_Initial)[1]
CloneDFIndices_List = list() 
for (Freq in SwitchFreq:1)	CloneDFIndices_List[[Freq]] = which(CloneDataFrame_Initial$Zi == Freq & CloneDataFrame_Initial$Age == 0 ) + 0:PromotionDataFrame_Initial[Freq, "HowManyTimeSteps"]

WhichClonesAreDeterministic		= which(CloneDataFrame_Initial[1:MaxPlausibleNoDeterministicBins, DorS] == "D")
WhichClonesAreStochastic		= which(CloneDataFrame_Initial[, DorS] == "S")

#### ==== #### ==== #### ==== #### ==== #### ==== #### ==== #### ==== #### ==== #### ==== #### ==== #### ==== #### ==== #### ==== #### ==== #### ==== #### ==== 
#### ==== #### ==== #### ==== #### ==== #### ==== #### ==== #### ==== #### ==== #### ==== #### ==== #### ==== #### ==== #### ==== #### ==== #### ==== #### ==== 

XaxisValues 		= Duration	 
YaxisValues 		= rep(Equilib_Diversity, length(XaxisValues))
dat 				= cbind(XaxisValues, YaxisValues)
colnames(dat) 		= c("time", "S")

BestToWorstFitIndices 			= order(ResultsDFrame$ModCost)
HowManyGuessesToInformBounds 	= min(c(dim(ResultsDFrame)[1], 10))

ResultsDFrame[BestToWorstFitIndices[1], ]
BestGuess		= ResultsDFrame$r_I[BestToWorstFitIndices[1]]
LowerBounds 	= min(ResultsDFrame$r_I  [BestToWorstFitIndices  [1:HowManyGuessesToInformBounds]], na.rm = TRUE)
UpperBounds 	= max(ResultsDFrame$r_I  [BestToWorstFitIndices  [1:HowManyGuessesToInformBounds]], na.rm = TRUE)

S_t_LowerBounds 	= ResultsDFrame$S_t[which(ResultsDFrame$r_I == LowerBounds)]
S_t_UpperBounds 	= ResultsDFrame$S_t[which(ResultsDFrame$r_I == UpperBounds)]

log10_BestGuess		= log10(BestGuess	)
log10_LowerBounds 	= log10(LowerBounds )
log10_UpperBounds 	= log10(UpperBounds	) 

##### to avoid numerical errors with logging and exponentiating r_I
LowerBounds = 0.98 * LowerBounds
UpperBounds = 1.02 * UpperBounds

RootHybrid 	<- function(Equilib_Diversity)  
{
	HybridCost <- function(log10_p)
	{
		real_rI = 10^log10_p
		
		iteration <<- iteration + 1
		cat(paste0("i", iteration, " log10_rI = ", signif(log10_p, 3), ", rI = " , signif(real_rI,3), " ", Sys.time(), "\n"))
		
		if(any(real_rI < LowerBounds) | any(real_rI > UpperBounds))   Ypred = rep(Inf, length(XaxisValues)) else
		{
			CloneDataFrame 		<<- copy(CloneDataFrame_Initial) 		## want CloneDataFrame as a global
			PromotionDataFrame 	<<- copy(PromotionDataFrame_Initial) 	## want PromotionDataFrame as a global
			
			log10_p 	<<- log10_p
			real_rI 	<<- real_rI
			HybridOut	<<- HybridFun(r_I = real_rI)
			Ypred 		= HybridOut$DeterS_t[paste0("t=", Duration)]
			if (is.na(Ypred)) Ypred = HybridOut$StochS_t[paste0("t=", Duration)]
			if(any(Ypred == "NaN")) 	Ypred = rep(Inf, length(XaxisValues)) 
		}
		Discrep 	<<- Ypred - Equilib_Diversity
		
		if(abs(Discrep) < abs(CurrentBest_Discrep)) 
		{
			CurrentBestR_I_log10 	<<- log10_p
			CurrentBestR_I 			<<- real_rI
			CurrentBest_Discrep		<<- Discrep
			HybridOut_best			<<- HybridOut
		} 
		cat(paste0("i", iteration, " finished: Discrep ", signif(Discrep, 3), " ", Sys.time(), "\n"))
		
		### Save fit information 
		BestFitInfo = c(CurrentBestR_I = CurrentBestR_I, CurrentBestR_I_log10 = CurrentBestR_I_log10, CurrentBest_Discrep = CurrentBest_Discrep)
		write.table(t(BestFitInfo), file = file.path(OutputDir, paste0(OutputStringPrefix, "Pt_", dta, "_BestFitInfo.txt")), sep = "\t", col.names = T, row.names = FALSE)
		
		return(Discrep) 
	}
}
f <- RootHybrid(Equilib_Diversity)

CurrentBestR_I 			= BestGuess
CurrentBestR_I_log10	= log10_BestGuess
CurrentBest_Discrep		= Inf
iteration 				= 0

ans <- uniroot(f, c(log10_LowerBounds, log10_UpperBounds))
cat("uniroot DONE\n")
HybridOut			<<- HybridOut_best


## Extract and write time-series from Hybrid
DeterTimeSeries = data.frame(Time = seq(0, Duration + TimeStep/2 + 2	, by = TimeStep/2), 
		DeterS_t = HybridOut$DeterS_t, DeterN_t = HybridOut$DeterN_t, DeterInc = HybridOut$DeterInc)
write.table(DeterTimeSeries, file = file.path(OutputDir, paste0(OutputStringPrefix, "Pt_", dta, "_DeterHybridOut.txt")), sep = "\t", col.names = T, row.names = FALSE)

StochTimeSeries = data.frame(Time = seq(0, Duration + TimeStep + 1	, by = TimeStep), 
		StochS_t = HybridOut$StochS_t, StochN_t = HybridOut$StochN_t, StochInc = HybridOut$StochInc)
write.table(StochTimeSeries, file = file.path(OutputDir, paste0(OutputStringPrefix, "Pt_", dta, "_StochHybridOut.txt")), sep = "\t", col.names = T, row.names = FALSE)

warnings()














