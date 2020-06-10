## For given patient data set, this script calculates look up tables required for the hybrid model

require(Matrix) 	
require(data.table) 

if (!exists("OutputStringPrefix")) OutputStringPrefix = ""

FreqChar				= "Zi"
Timepoints 				= seq(0, Duration + TimeStep, by = TimeStep)
MaxStochasticDuration	= max(Timepoints) + 1 

## import functions
source(file.path(R_Script_Dir, "CloneBirthDeath.R"))

### === ### === ### === LOAD & PROCESS ESTIMATED DISTRIBUTION
EstDistFileName = paste0("EstCloneFreqBody_Sample_", dta, ".txt")
DATA = read.table(file = file.path(RootDir, "Inputs", "EstDists", EstDistFileName), header = T, sep = "\t")
DATA = DATA[DATA$Frequency <= FreqCutOff, ] ### set FreqCutOff to Inf if not cutting. 
str(DATA)
DATA$CloneNumber	= NULL 
rownames(DATA) 		= paste0("C_", 1:dim(DATA)[1])
colnames(DATA) 		= "Zi"
str(DATA)

MagBin = floor(log10(DATA[, "Zi"])) + SwitchFreq + 1	
for (Freq in 1:SwitchFreq) MagBin[which(DATA[, "Zi"] == Freq)] = Freq 

DATA = cbind(DATA, MagBin)
if (!exists("Equilib_Diversity")) Equilib_Diversity = dim(DATA)[1]

NumberMagBins	 	= max(MagBin)	
NoClonesInBin		= c()
AverageFreqInBin 	= c()

for (Freq in 1:SwitchFreq)
{
	NoClonesInBin	[Freq] = length(which(MagBin == (Freq)))
	AverageFreqInBin[Freq] = Freq
}
if (NumberMagBins > (SwitchFreq+1))
	for (magnitude in (SwitchFreq+1):NumberMagBins)
	{
		NoClonesInBin	[magnitude] = length(which(MagBin == magnitude))
		AverageFreqInBin[magnitude] = mean(DATA[which(MagBin == magnitude), "Zi"])
	}
AverageFreqInBin 	= rev(AverageFreqInBin)
NoClonesInBin		= rev(NoClonesInBin)
SummaryCloneDataFrame = data.frame(Bin = 1:NumberMagBins, NoClonesInBin, AverageFreqInBin)
names(SummaryCloneDataFrame)= c("Bin", "NoClonesInBin", FreqChar)
rm(Freq, magnitude, MagBin, NumberMagBins, NoClonesInBin, AverageFreqInBin)

### === ### === ### === DEFINE (person-specific) PROLIFERATION AND DEATH PARAMETERS

N_0 		= sum(DATA[, "Zi"])	
Pi_dash2 	= 1
TotLymphoBody	= 2e12
CD3overLympho	= 0.75
CD4overCD3		= 0.7
ROposoverCD4	= 0.5
TotROposBody	= TotLymphoBody * CD3overLympho * CD4overCD3 * ROposoverCD4 

## Equations to solve are: 
# 1) Pi / (K + U + I + TotROposBody) = Pi_dash1 = 0.0316, where U = number of uninfected cells (person-specific) and K is the carrying capacity (not person specific). 
# 2) Pi / (K + U + I) 				 = Pi_dash2 = 1
## Therefore K + U = Pi from 2)
Pi 		= (Pi_dash1 * TotROposBody)/(1 - Pi_dash1)
K 		= Pi
K 		= KplusU = K + TotROposBody - N_0
H 		= K + N_0
ParamVec = c(Tau = Tau, SwitchFreq = SwitchFreq, Pi = Pi, delta = delta, H = H, K = K, TimeStep = TimeStep, Duration = Duration, Equilib_Diversity = Equilib_Diversity)

### === ### === ### === Make precalculated clone trajectories
DefineCloneDataFrame = function()
{
	BaseDFrame = SummaryCloneDataFrame
	CloneDataFrame = cbind( BaseDFrame, 
			Age 			= rep(0 , dim(BaseDFrame)[1]),
			ExtinctProb 	= rep(0 , dim(BaseDFrame)[1]),
			DorS 			= rep(NA, dim(BaseDFrame)[1]), 
			Var_Zi 			= rep(0 , dim(BaseDFrame)[1]))
	
	DetermIndices				= which(CloneDataFrame[,FreqChar] > SwitchFreq)
	DetermCloneDataFrame 		= CloneDataFrame[DetermIndices, ]
	if (length(DetermIndices) > 0) DetermCloneDataFrame$DorS	= "D"
	StochastCloneDataFrame 		= CloneDataFrame[which(CloneDataFrame[,FreqChar] <= SwitchFreq), ]
	
	CloneDataFrame 				= DetermCloneDataFrame 
	if (length(DetermIndices) > 0) Bin = tail(DetermCloneDataFrame$Bin, 1) else Bin = 0 
	for (Freq in SwitchFreq:1) ## reverse ordered
	{
		CloneProgressionList 	= ProgressionsList[[paste("CloneEvo_StartOff_", Freq, sep ="")]]
		HowManyAges 			= length(CloneProgressionList$Means)	
		HowManyAlreadyExist 	= round(StochastCloneDataFrame[which(StochastCloneDataFrame$Zi == Freq), "NoClonesInBin"] ) 
		
		AddOn 			= matrix(nrow = HowManyAges + 1, ncol = dim(CloneDataFrame)[2])		
		AddOn 			= as.data.frame(AddOn)
		colnames(AddOn) = colnames(CloneDataFrame)
		
		## Populate Stochastic Summary
		AddOn[, "Bin"			]	= Bin + 1:(HowManyAges + 1) ## want an age zero for new clones
		AddOn[, "NoClonesInBin"	]	= c(HowManyAlreadyExist	, rep(0, HowManyAges))
		AddOn[, "Age"			]	= c(0, 1:HowManyAges * TimeStep)
		AddOn[, "Zi"			]	= c(Freq, CloneProgressionList$Means)
		AddOn[, "DorS"			]	= "S"
		AddOn[, "ExtinctProb"	]	= c(0, CloneProgressionList$PDFs[1,])	
		AddOn[, "Var_Zi" 		]	= c(0, CloneProgressionList$Vars)
		
		Bin = tail(AddOn$Bin, 1)
		CloneDataFrame = rbind(CloneDataFrame, AddOn)
	}
	CloneDataFrame = as.data.table(CloneDataFrame)
	return(CloneDataFrame)
}
DefinePromotionDataFrame = function()
{
	PromotionDataFrame = data.frame(	
			X_0					= 1:SwitchFreq						,
			HowManyTimeSteps	= rep(NA		, SwitchFreq)		,
			Promoted			= rep(NA		, SwitchFreq)		,
			LegallyDead			= rep(NA		, SwitchFreq)		,
			NoPromoted			= rep(0			, SwitchFreq)		,
			ExtProb				= rep(NA		, SwitchFreq)		,
			NZi					= rep(NA		, SwitchFreq)		,
			NVar_Zi 			= rep(NA		, SwitchFreq)		,  
			OldVNewZiDiscrep	= rep(NA		, SwitchFreq)		)
	
	for (Freq in SwitchFreq:1) 
	{
		CloneProgressionList = ProgressionsList[[paste("CloneEvo_StartOff_", Freq, sep ="")]]
		HowManyAges 		 = length(CloneProgressionList$Means)	## as in, "how many before you promoted"?
		
		## Populate Promotion Data Frame
		Index = which(PromotionDataFrame[, "X_0"] == Freq)
		PromotionDataFrame[Index, "Promoted"			]   = CloneProgressionList$PromotionConditions["Promoted"]
		PromotionDataFrame[Index, "LegallyDead"			]   = CloneProgressionList$PromotionConditions["LegallyDead"]
		PromotionDataFrame[Index, "HowManyTimeSteps"	]   = HowManyAges
		NormDummyProg = NormalizeMomentsAfterExtinction(CloneProgressionList, HowManyAges)
		
		PromotionDataFrame[Index, "ExtProb"				]   = NormDummyProg$ExtProb
		PromotionDataFrame[Index, "NZi"					]   = NormDummyProg$Mean
		PromotionDataFrame[Index, "NVar_Zi"				]	= NormDummyProg$Var
		PromotionDataFrame[Index, "OldVNewZiDiscrep"	]   = NormDummyProg$Mean - CloneProgressionList$Means[HowManyAges]
	}
	return(PromotionDataFrame)
}

Sys.time()
RunStartTime 		= proc.time()
ProgressionsList	= CalculatePDFEvolutions()		
RunRunTime 			= proc.time() - RunStartTime
RunRunTime
Sys.time()

Sys.time()
CloneDataFrame_Initial 		= DefineCloneDataFrame()
PromotionDataFrame_Initial	= DefinePromotionDataFrame() 	
Sys.time()

write.table(t(ParamVec)					, file = file.path(OutputDir, paste0(OutputStringPrefix, "Pt_", dta, "_ParamVec.txt"))		, sep = "\t", col.names = T, row.names = FALSE)
write.table(CloneDataFrame_Initial		, file = file.path(OutputDir, paste0(OutputStringPrefix, "CloneLookUpTable_", dta, ".txt"))	, sep = "\t", col.names = T, row.names = FALSE)
write.table(PromotionDataFrame_Initial	, file = file.path(OutputDir, paste0(OutputStringPrefix, "PromotionDFrame_", dta, ".txt"))	, sep = "\t", col.names = T, row.names = FALSE)

warnings()









