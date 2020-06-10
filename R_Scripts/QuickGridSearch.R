
require(Matrix) 	
require(data.table) 
require(deSolve)

if (!exists("OutputStringPrefix")) OutputStringPrefix = ""

ParamVec 		= read.table(file = file.path(OutputDir, paste0(OutputStringPrefix, "Pt_", dta, "_ParamVec.txt")), header = T, sep = "\t")

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

#### pre-calculate quantities to reduce runtime.  
MaxPlausibleNoDeterministicBins = 15  ### for calculating WhichClonesAreDeterministic quickly (no need to search entire data frame) 
NoRows_CloneDF 					= dim(CloneDataFrame_Initial)[1]
CloneDFIndices_List = list()  
for (Freq in SwitchFreq:1) CloneDFIndices_List[[Freq]] = which(CloneDataFrame_Initial$Zi == Freq & CloneDataFrame_Initial$Age == 0 ) + 0:PromotionDataFrame_Initial[Freq, "HowManyTimeSteps"]

WhichClonesAreDeterministic	<<- which(CloneDataFrame_Initial[1:MaxPlausibleNoDeterministicBins, DorS] == "D")
WhichClonesAreStochastic	<<- which(CloneDataFrame_Initial[, DorS] == "S")

#### ==== #### ==== #### ==== #### ==== #### ==== #### ==== #### ==== #### ==== #### ==== #### ==== #### ==== #### ==== #### ==== #### ==== #### ==== #### ==== 
#### ==== #### ==== #### ==== #### ==== #### ==== #### ==== #### ==== #### ==== #### ==== #### ==== #### ==== #### ==== #### ==== #### ==== #### ==== #### ==== 

source(file.path(R_Script_Dir, "HybridFunctions.R"))

if (!exists("NumberOfGuesses"	)) NumberOfGuesses 	= 25
if (!exists("BestGuess"			)) BestGuess 		= 1e-9
if (!exists("LowerBound"		)) LowerBound 		= 1e-11
if (!exists("UpperBound"		)) UpperBound 		= 3e-9	

Ri_ValuesToTry 		= 10^(seq(log10(UpperBound), log10(LowerBound), length.out = NumberOfGuesses))
CurrentBestR_I 		= BestGuess
CurrentBestSOS 		= Inf
ParamGridFileName	= file.path(OutputDir, paste0(OutputStringPrefix, "ParamGridFit_", dta, ".txt"))

param_index = 1
for (param_index in 1:NumberOfGuesses)
{
	rIguess = Ri_ValuesToTry[param_index]
	cat(paste0(Sys.time(), " guess ", param_index, ": rIguess = ", signif(rIguess, 3), "\n"))
	
	## Run hybrid
	CloneDataFrame 		<<- copy(CloneDataFrame_Initial) 
	PromotionDataFrame 	<<- copy(PromotionDataFrame_Initial) 	 
	HybridOut			= HybridFun(r_I = rIguess)
	
	Ypred = as.numeric(HybridOut$DeterS_t[paste0("t=", Duration)])
	if (is.na(Ypred)) Ypred = HybridOut$StochS_t[paste0("t=", Duration)]
	SumOfSquares 		= 2 * ((Ypred - Equilib_Diversity)^2)
	Discrep 			= Ypred - Equilib_Diversity
	
	## store results of fit
	temp_dframe = data.frame(r_I = rIguess, S_t = Ypred, Discrep, ModCost = SumOfSquares)
	if (param_index == 1) ResultsDFrame = temp_dframe	else 	ResultsDFrame = rbind(ResultsDFrame, temp_dframe)
	write.table(ResultsDFrame, file = ParamGridFileName, sep = "\t", col.names = T, row.names = FALSE)
	
	cat(paste0(Sys.time(), " guess ", param_index, " finished. Ypred = ", round(Ypred), " Discrep = ", round(Discrep), "\n\n"))
}

warnings()







