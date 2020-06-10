## Author: dlaydon 

## Script contains functions for hybrid model of chronic HTLV-1 infection/within host persistence.
## Model divides HTLV-1 proviral load (number of infected cells) into clones.
## Hybrid model is comprised of two systems; i) deterministic system modelled by series of ODES (for large clones); ii) stochastic system modelled by multiple birth-death processes (for smaller clones).
## Clones proliferate via mitotic spread. Clones are created by infectious spread.
## Model works by updating two data tables created outside this script: CloneDataFrame and PromotionDataFrame (global variables). These are look up tables correspnding to "clone trajectories"...
## ... i.e. summary statistics (mean, variance and extinction probabilities) of multiple probability distributions, which are solutions of the master equation governing a birth-death process.
## For example, for a clone that started with frequency X, what are its mean, variance, and extinction probability after 1 day, 2 days etc.? 
## This information is necessary to calculate the expected number of clones and infected cells at time t.


`%MultiAssign%` <- function(x, y)
{
	mapply(assign, x, y, MoreArgs = list(envir = parent.frame()))
	invisible()
}

BinnedCloneODEs = function(t, state, parameters, NoStochasticCells, TotalNormalisedMeansDiscrepancy) ## "state" refers to "deterministic state", not state of entire deterministic and stochastic systems. 
{
	with(as.list(c(state, parameters)),
	{
		SumOfAllCellsPlusK = sum(as.numeric(mget(paste0("NoInBin_", 1:length(state)))) *  
						as.numeric(mget(paste0("EZ_", 1:length(state))	))) + NoStochasticCells + K - TotalNormalisedMeansDiscrepancy 
		## Create clone derivatives using multi assign. 
		paste0("dEZ_", 1:length(state))  %MultiAssign% 	(((Pi/(SumOfAllCellsPlusK)) - delta) * as.numeric(mget(paste0("EZ_", 1:length(state)) )	))
		return(list(as.numeric(mget(paste0("dEZ_", 1:length(state)) 	))))
	})
}
## Counting functions for cells and clones. 
CalculateExpectedNoStochasticCells 	= function()	CloneDataFrame[ WhichClonesAreStochastic, sum(Zi * NoClonesInBin)]
CalculateExpectedNoClones 			= function()
{
	HowManyPromotedAreExtinct 			= sum(PromotionDataFrame$NoPromoted * PromotionDataFrame$ExtProb)
	HowManyStochasticsAreExtinct 		= CloneDataFrame[, sum(NoClonesInBin * ExtinctProb)]
	
	return(CloneDataFrame[, sum(NoClonesInBin)] - HowManyPromotedAreExtinct - HowManyStochasticsAreExtinct)
}
CalculateTotalNormalisedMeansDiscrepancy 	= function() sum(PromotionDataFrame$NoPromoted * PromotionDataFrame$OldVNewZiiDiscrep)
CalculateHowManyInfectedCells 				= function() CloneDataFrame[ , sum(Zi * NoClonesInBin)] - CalculateTotalNormalisedMeansDiscrepancy()		

DetermHalfTimeStep 		= function()
{
	if (length(WhichClonesAreDeterministic) > 0)	
	{
		ODE_In 			= CloneDataFrame[WhichClonesAreDeterministic,Zi]
		names(ODE_In) 	= paste0("EZ_" , WhichClonesAreDeterministic)
		
		## ODEparameters must be updated between timesteps to include number of clones in each bin
		ODEparameters 			= c(Pi, delta, K, CloneDataFrame[WhichClonesAreDeterministic,NoClonesInBin]	)	
		names(ODEparameters)	= c("Pi", "delta", "K",  paste0("NoInBin_", 1:length(WhichClonesAreDeterministic)))
		TotalNormalisedMeansDiscrepancy = sum(PromotionDataFrame$NoPromoted		 * 		PromotionDataFrame$OldVNewZiDiscrep)
	
		ODE_Pars_AddOn = c(NoStochasticCells = CalculateExpectedNoStochasticCells(), TotalNormalisedMeansDiscrepancy = TotalNormalisedMeansDiscrepancy)
		ODE_Out = ode(y = ODE_In, times = c(0, (TimeStep/2)), func = BinnedCloneODEs, parms = c(ODEparameters, ODE_Pars_AddOn))	
		
		## Amend CloneDataFrame
		set (CloneDataFrame, i = WhichClonesAreDeterministic, j = "Age"	, CloneDataFrame[WhichClonesAreDeterministic, Age] + 0.5	)		
		set (CloneDataFrame, i = WhichClonesAreDeterministic, j = "Zi"	, ODE_Out[-1, grep("EZ_" , colnames(ODE_Out))])		
	}
}
DemotionFunction		= function()
{
	if (length(WhichClonesAreDeterministic) > 0)	
	{
		WhichClonesSankBelow = WhichClonesAreDeterministic[which(CloneDataFrame[WhichClonesAreDeterministic, Zi	] < SwitchFreq)]	
		if(length(WhichClonesSankBelow) > 0)
			for (CloneCategory in WhichClonesSankBelow)
			{
				HowManyInCategory = CloneDataFrame$NoClonesInBin[CloneCategory]
				
				CloneCategoryFreq = round(CloneDataFrame[CloneCategory, Zi])
				Vec = which(CloneDataFrame[, Zi] == CloneCategoryFreq 	& 	CloneDataFrame$Age == 0)
				set (CloneDataFrame, i = as.integer(Vec), j = "NoClonesInBin",  CloneDataFrame[Vec, NoClonesInBin] +  HowManyInCategory    )		
	
				## Subract out of original bin
				set(CloneDataFrame, i = as.integer(CloneCategory), j = "NoClonesInBin", 0)
			}
	}
}
CloneCounterFun 		= function(ContinousCount, DiscreteCount, Increment)
{
	a = c(	ContinousCount + Increment							,	## Running Total continuous
			floor(ContinousCount + Increment) 					, 	## Running Total discrete
			floor(ContinousCount - DiscreteCount + Increment) 	)
	names(a) = c("ContinousCount", "DiscreteCount", "NumNewClones")
	return(a)
}
InfectiousSpreadFun		= function(Determ_Or_Stochastic_Step, ContinousCount, DiscreteCount, r_I)
{
	if(Determ_Or_Stochastic_Step == "Deterministic")
	{
		HowManyInfectingCells	= CloneDataFrame[WhichClonesAreDeterministic, sum(Zi * NoClonesInBin)] - CalculateTotalNormalisedMeansDiscrepancy()	
		LengthOfTimeStep 		= TimeStep/2		
	
	} else ## i.e. Determ_Or_Stochastic_Step == "Stochastic"
	{
		HowManyInfectingCells 	= CalculateExpectedNoStochasticCells()			
		LengthOfTimeStep 		= TimeStep			
	}
	
	IncrementToCloneCount = r_I * HowManyInfectingCells * LengthOfTimeStep
	CloneCount 		= CloneCounterFun(ContinousCount, DiscreteCount, Increment = IncrementToCloneCount)
	ContinousCount 	= CloneCount["ContinousCount"	]
	DiscreteCount 	= CloneCount["DiscreteCount"	]						
	NumNewClones	= CloneCount["NumNewClones"		]						
	
	if (NumNewClones > 0)
	{
		SingletonVec = which(CloneDataFrame$Zi == 1 & CloneDataFrame$Age == 0) # age = 0 : Zi == 1 because they are singletons 
		set (CloneDataFrame, i = SingletonVec, j = "NoClonesInBin", CloneDataFrame[SingletonVec, NoClonesInBin] + NumNewClones)		
	}
	
	List = list()
	List[["ContinousCount"			]] = ContinousCount
	List[["DiscreteCount"			]] = DiscreteCount
	List[["IncrementToCloneCount"	]] = IncrementToCloneCount
	return(List)
}
StochasticWholeTimeStep	= function()
{
	if (length(WhichClonesAreStochastic) > 0)
	{
		## cycle through frequencies, and move clones from one age group to the next
		for (Freq in SwitchFreq:1) ## reverse ordered
		{
			NoClonesInEachAgeGroupForThisCategory 	= CloneDataFrame$NoClonesInBin[CloneDFIndices_List[[Freq]]]
			NoClonesPromotedInThisTimeStep 			= 0
			
			if(PromotionDataFrame[Freq, "Promoted"] == TRUE)		
			{	
				## Promote oldest clone
				NoClonesPromotedInThisTimeStep 			= tail(NoClonesInEachAgeGroupForThisCategory, 1) 
				## Age all other clones by one year, and there are now zero age 0 clones (at least until next deterministic half time step). 
				## redefine NoClonesInEachAgeGroupForThisCategory
				NoClonesInEachAgeGroupForThisCategory = c(0, NoClonesInEachAgeGroupForThisCategory[-length(NoClonesInEachAgeGroupForThisCategory)])
			}
			if((PromotionDataFrame[Freq, "LegallyDead"] == TRUE ) || (PromotionDataFrame[Freq, "HowManyTimeSteps"] >= length(Timepoints)))		
			{
				NoClonesPromotedInThisTimeStep = 0
				## Age all other clones by one year, and there are now zero age 0 clones (at least until next deterministic half time step). 
				NoClonesInEachAgeGroupForThisCategory = c(0, NoClonesInEachAgeGroupForThisCategory[-length(NoClonesInEachAgeGroupForThisCategory)])
			}
			set (CloneDataFrame, i = CloneDFIndices_List[[Freq]], j = "NoClonesInBin",  NoClonesInEachAgeGroupForThisCategory)
			
			### Promote 
			if (NoClonesPromotedInThisTimeStep > 0)
			{
				NZMean 		= PromotionDataFrame[Freq, "NZi"] ## Normalized Mean of this
				WhichBin 	= max(which( (CloneDataFrame$DorS == "D") &  (CloneDataFrame$Zi > NZMean)))	## Which Deterministic bin will you amend?
				
				# what is the average freqeuncy, variance etc. (for zi) of the bin, with new clones (weighted average, where weights are the number of clones already in the bin, and those being added to the bin). 
				CloneDataFrame[WhichBin, "Zi"] 	<<- (((CloneDataFrame[, "Zi"] * CloneDataFrame[, "NoClonesInBin"])[WhichBin]) +	(PromotionDataFrame[Freq, "NZi"] * NoClonesPromotedInThisTimeStep)) / (NoClonesPromotedInThisTimeStep + CloneDataFrame[WhichBin, "NoClonesInBin"])
				
				## add to number of clones in this bin. (recursive defn)
				set(CloneDataFrame, i = as.integer(WhichBin), j = "NoClonesInBin",  CloneDataFrame[WhichBin, NoClonesInBin] + NoClonesPromotedInThisTimeStep)
				
				## Record number of clones have been promoted over all time. 
				PromotionDataFrame[Freq, "NoPromoted"] <<- PromotionDataFrame[Freq, "NoPromoted"] + NoClonesPromotedInThisTimeStep
				
				rm(NZMean, WhichBin)
			}
		}		
	}
}
HybridFun				= function(r_I)
{
	###### Initialize Clone Count Quantities
	ContinousCount = DiscreteCount = 0		
	
	Names_DeterS_t = paste0("t=",as.character(seq(0, Duration + TimeStep/2 + 2	, by = TimeStep/2))) 
	Names_DeterN_t = paste0("t=",as.character(seq(0, Duration + TimeStep/2 + 2	, by = TimeStep/2))) 
	Names_DeterInc = paste0("t=",as.character(seq(0, Duration + TimeStep/2 + 2	, by = TimeStep/2))) 
	Names_StochS_t = paste0("t=",as.character(seq(0, Duration + TimeStep + 1	, by = TimeStep)))   
	Names_StochN_t = paste0("t=",as.character(seq(0, Duration + TimeStep + 1	, by = TimeStep)))   
	Names_StochInc = paste0("t=",as.character(seq(0, Duration + TimeStep + 1	, by = TimeStep)))   
	
	### Expected no. clones over time
	DeterS_t 	= rep(NA, length(Names_DeterS_t)) # records expected number of clones at deterministic half timesteps
	StochS_t 	= rep(NA, length(Names_StochS_t)) # records expected number of clones at stochastic full timesteps 
	
	### Expected no. cells over time
	DeterN_t 	= rep(NA, length(Names_DeterN_t)) # records expected number of infected cells at deterministic half timesteps
	StochN_t 	= rep(NA, length(Names_StochN_t)) # records expected number of infected cells at stochastic full timesteps
	
	### Increments to Clone Count - how many new clones born at half timesteps (deterministic) and whole timesteps (stochastic)
	DeterInc 	= rep(NA, length(Names_DeterInc)) #  how many new clones born at half timesteps (deterministic) 
	StochInc 	= rep(NA, length(Names_StochInc)) #  how many new clones born at whole timesteps (stochastic)
	
	###### Begin Loop
	for (Step in 0:length(Timepoints))
	{
		###### ====== ###### ====== ###### ====== ###### ====== ###### ====== ###### ====== ###### ====== ###### ====== ###### ====== ###### ====== ###### ====== ###### ====== ###### ====== ###### ====== ###### ====== 
		###### ## first half time step - deterministic 
		###### ====== ###### ====== ###### ====== ###### ====== ###### ====== ###### ====== ###### ====== ###### ====== ###### ====== ###### ====== ###### ====== ###### ====== ###### ====== ###### ====== ###### ====== 
		
		if (length(WhichClonesAreDeterministic) > 0)
		{
			DeterS_t[(Step*2) + 1] 	= CalculateExpectedNoClones()
			DeterN_t[(Step*2) + 1] 	= CalculateHowManyInfectedCells()
			
			DetermHalfTimeStep()
			DemotionFunction() 
			Dummy = InfectiousSpreadFun(Determ_Or_Stochastic_Step = "Deterministic", ContinousCount, DiscreteCount, r_I = r_I)
			
			ContinousCount 			= Dummy[["ContinousCount"		]]
			DiscreteCount 			= Dummy[["DiscreteCount"		]]
			DeterInc[(Step*2) + 1] 	= Dummy[["IncrementToCloneCount"]]
			rm(Dummy)
			
			DeterS_t[(Step*2) + 2] 	= CalculateExpectedNoClones()
			DeterN_t[(Step*2) + 2] 	= CalculateHowManyInfectedCells()
		}
		
		###### ====== ###### ====== ###### ====== ###### ====== ###### ====== ###### ====== ###### ====== ###### ====== ###### ====== ###### ====== ###### ====== ###### ====== ###### ====== ###### ====== ###### ====== 
		###### ## Move stochastic part of system forward by a whole timestep - keep deterministic constant. 

		if (length(WhichClonesAreStochastic) > 0)
		{
			StochasticWholeTimeStep()
			StochS_t[Step+1] 	= CalculateExpectedNoClones()
			StochN_t[Step+1] 	= CalculateHowManyInfectedCells()
			
			### Infectious Spread - here function considers only stochastically modelled cells
			Dummy = InfectiousSpreadFun(Determ_Or_Stochastic_Step = "Stochastic", ContinousCount, DiscreteCount, r_I = r_I)
			
			ContinousCount 		= Dummy[["ContinousCount"		]]
			DiscreteCount 		= Dummy[["DiscreteCount"		]]
			StochInc[Step+1] 	= Dummy[["IncrementToCloneCount"]]
			rm(Dummy)
		}
		
		###### ====== ###### ====== ###### ====== ###### ====== ###### ====== ###### ====== ###### ====== ###### ====== ###### ====== ###### ====== ###### ====== ###### ====== ###### ====== ###### ====== ###### ====== 
		## second deterministic half time step - keep stochastic constant
		###### ====== ###### ====== ###### ====== ###### ====== ###### ====== ###### ====== ###### ====== ###### ====== ###### ====== ###### ====== ###### ====== ###### ====== ###### ====== ###### ====== ###### ====== 
		
		if (length(WhichClonesAreDeterministic) > 0)
		{
			DetermHalfTimeStep()
			DemotionFunction()
			Dummy = InfectiousSpreadFun(Determ_Or_Stochastic_Step = "Deterministic", ContinousCount, DiscreteCount, r_I = r_I)
			
			ContinousCount 			= Dummy[["ContinousCount"		]]
			DiscreteCount 			= Dummy[["DiscreteCount"		]]
			DeterInc[(Step*2) + 2] 	= Dummy[["IncrementToCloneCount"]]
			rm(Dummy)
		}
	}	
	names(DeterS_t) = Names_DeterS_t 
	names(DeterN_t) = Names_DeterN_t 
	names(DeterInc) = Names_DeterInc 
	names(StochS_t) = Names_StochS_t 
	names(StochN_t) = Names_StochN_t 
	names(StochInc) = Names_StochInc 
	
	List = list()
	List[["r_I"]]						= r_I
	List[["ContinousCount"]]			= ContinousCount
	List[["DiscreteCount"]]				= DiscreteCount
	List[["DeterS_t"]]					= DeterS_t
	List[["DeterN_t"]]					= DeterN_t
	List[["DeterInc"]]					= DeterInc
	List[["StochS_t"]]					= StochS_t
	List[["StochN_t"]]					= StochN_t
	List[["StochInc"]]					= StochInc
	
	return(List)
}


