## Contains functions for birth-death process. 
### MakeTransitionMatrix: Make Transition matrix
### MakeP_0: Make initial probability distribution P_0
### CalculateMatrixExponential: Matrix Exponential (requires package "Matrix") of transition matrix
### MakeP_tstep: Make single probability distribution P_t given timestep
### PDFEvolution: Recursively calculate P_t for multiple values of t, from given starting frequency.
### NormalizeMomentsAfterExtinction: calculates conditional distribution (on assumpmtion that clone/species has "survived") what are means and variances?
### CalculatePDFEvolutions: Calls PDFEvolution repeatedly for multiple clones from multiple starting frequencies. 

MakeTransitionMatrix 			= function(Tau, Pi, delta, H, ArbitraryColumnErrorTolerance = 10^-10)
{
	## Create empty sparse matrix
	TransitionMatrix = Matrix(0 , nrow = Tau + 1, ncol = Tau + 1, sparse = TRUE)
	
	# Populate diagaonal, super-diagonal and sub=diagonal
	diag(TransitionMatrix) 					= c(0, ((1:(Tau-1)) * -((Pi/ (H + (1:(Tau-1)))	) + delta)), -delta * Tau)
	diag(TransitionMatrix[-(Tau+1),-1])  	= 1:Tau * delta
	diag(TransitionMatrix[-1,-(Tau+1)]) 	= c(0, (1:(Tau-1) * Pi)  / (H + (1:(Tau-1)))	)
	
	if (any(abs(colSums(TransitionMatrix)) > ArbitraryColumnErrorTolerance)) 	 	
		stop("MakeTransitionMatrix error: columns do not sum to 0")
	
	return(TransitionMatrix)
}
MakeP_0 						= function(Tau, StartOff) 
{
	if(StartOff > Tau) stop(paste0("MakeP_0 error: StartOff > Tau: StartOff = ", StartOff, ", Tau = ", Tau))
	P_0 	= matrix(0, ncol = 1, nrow = Tau + 1)		
	P_0[StartOff + 1] = 1
	return(P_0)
}
CalculateMatrixExponential 		= function(TransitionMatrix, TimeStep)
{
	eAt 					= expm(TransitionMatrix * TimeStep)
	List 					= list()
	List[["eAt"]] 			= eAt
	List[["TimeStep"]] 		= TimeStep
	List[["call"]]			= match.call()
	return (List)
}
MakeP_tstep 					= function(Tau, Pi, delta, H, Time, TimeStep = 1 , ArbitraryColumnErrorTolerance = 10^-10, StartOff, eAt = NULL, P_0 = NULL, TransitionMatrix = NULL)
{
	if (eAt$TimeStep != TimeStep) stop ("MakeP_tstep error: Timesteps do not match")
	if (is.null(TransitionMatrix)) 	TransitionMatrix 	= MakeTransitionMatrix(Tau = Tau, Pi = Pi, delta = delta, H = H, ArbitraryColumnErrorTolerance = ArbitraryColumnErrorTolerance)
	if (is.null(P_0))				P_0 				= MakeP_0(Tau = Tau, StartOff = StartOff)
	if (is.null(eAt))				eAt					= CalculateMatrixExponential(TransitionMatrix = TransitionMatrix, TimeStep = TimeStep)
	
	PDF 		= eAt$eAt %*% P_0
	Mean  		= sum(0:(Tau) * PDF	)								
	Var			= sum(PDF * ((0:(Tau)) - Mean)^2)
	Mode		= which(PDF == max(PDF)) - 1
	ParamVec 	= c(Pi = Pi, delta = delta, H = H, Tau = Tau, TimeStep = TimeStep, StartOff = StartOff)
	
	PDFlist 							= list()
	PDFlist[["PDF"	]] 					= PDF
	PDFlist[["Mean"	]]					= Mean
	PDFlist[["Var"	]]					= Var	
	PDFlist[["Mode"	]]					= Mode
	PDFlist[["ExtProb"	]]				= PDF[1]
	PDFlist[["TimeStep"	]] 				= TimeStep
	PDFlist[["ParamVec"	]] 				= ParamVec
	PDFlist[["MakeP_tstep_call"	]] 		= match.call()
	
	return(PDFlist)
}

PDFEvolution = function(Pi, delta, H, Tau, TimeStep = 1, Duration, StartOff, ArbitraryColumnErrorTolerance = 10^-10,	 
		eAt = NULL, P_0 = NULL, TransitionMatrix = NULL, P_tstep = NULL, LegallyDeadThreshold = 1, RubberNeckThreshold = 0.01)			

{
	if (eAt$TimeStep != TimeStep) 	stop("PDFEvolution error: Timesteps do not match")
	if (is.null(TransitionMatrix)) 	TransitionMatrix 	= MakeTransitionMatrix(Tau = Tau, Pi = Pi, delta = delta, H = H, ArbitraryColumnErrorTolerance = ArbitraryColumnErrorTolerance)
	if (is.null(P_0))				P_0 				= MakeP_0(Tau = Tau, StartOff = StartOff)
	if (is.null(eAt))				eAt					= CalculateMatrixExponential(TransitionMatrix = TransitionMatrix, TimeStep = TimeStep )
	if (is.null(P_tstep))			P_tstep 			= MakeP_tstep(Pi = Pi, delta = delta, H = H, Tau = Tau, TimeStep = TimeStep, StartOff = StartOff, eAt = eAt, P_0 = P_0, TransitionMatrix = TransitionMatrix )
	
	Timepoints 	= seq(TimeStep, Duration + TimeStep, by = TimeStep)
	
	PDFs 			= matrix(nrow = Tau + 1, ncol = length(Timepoints))
	rownames(PDFs) 	= rownames(P_tstep$PDF)
	colnames(PDFs)	= paste("t =", Timepoints)
	Means			= rep(NA, length(Timepoints))
	Vars			= rep(NA, length(Timepoints))
	Modes			= rep(NA, length(Timepoints))
	
	for (Time in 1:length(Timepoints))
	{
		if (Time == 1)	
			SingleTimeStepPDF 	= P_tstep 	else	
			SingleTimeStepPDF	= MakeP_tstep(Pi = Pi, delta = delta, H = H, Tau = Tau, TimeStep = TimeStep, StartOff = StartOff, 
					eAt = eAt, P_0 = SingleTimeStepPDF$PDF, TransitionMatrix = TransitionMatrix )
		
		PDFs	[, 	Time]	= as.numeric(SingleTimeStepPDF$PDF)
		Means	[	Time] 	= SingleTimeStepPDF$Mean	
		Vars 	[	Time] 	= SingleTimeStepPDF$Var	
		Modes	[	Time] 	= SingleTimeStepPDF$Mode	
		
		## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 
		## Define Promotion Conditions
		## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 
		
		DidCloneRubberneck 	= (PDFs[Tau + 1, Time] - PDFs[Tau, Time]) > RubberNeckThreshold
		## Is the clone "legally dead"? i.e. if the expected value of the clone is less than 1, and it's extinction prob is giant, declare it dead
		IsCloneDead 		= (Means[Time] < 1) & (PDFs[1, Time] > LegallyDeadThreshold)
		
		if (DidCloneRubberneck | IsCloneDead) break
	}
	
	## remove NAs if clone died or was promoted
	PDFs 	= PDFs	[ , 1:Time]
	Means	= Means	[	1:Time]
	Vars	= Vars	[	1:Time]
	Modes	= Modes	[	1:Time]
	
	ParamVec = c(	Pi = Pi, delta = delta, H = H, 	Tau = Tau, TimeStep = TimeStep, 
			StartOff 				= StartOff				, 
			LegallyDeadThreshold 	= LegallyDeadThreshold	,
			RubberNeckThreshold		= RubberNeckThreshold	)
	
	## DidCloneRubberneck = Promoted? IsCloneDead = Legally Dead
	PromotionConditions 		= c(DidCloneRubberneck, IsCloneDead)
	names(PromotionConditions) 	= c("Promoted", "LegallyDead")
	
	List = list()
	List[[ "PDFs"				]]		= PDFs
	List[[ "Means"				]]		= Means 
	List[[ "Vars"				]]		= Vars	 
	List[[ "Modes"				]]		= Modes
	List[[ "TimeStep"			]]		= TimeStep 
	List[[ "StartOff"			]]		= StartOff 
	List[[ "Tau"				]]		= Tau 
	List[[ "ParamVec"			]]		= ParamVec 
	List[[ "PromotionConditions"]]		= PromotionConditions 
	List[[ "PDFEvolution_call"	]]		= match.call()
	
	return(List)
}	
NormalizeMomentsAfterExtinction = function(PDF_PROG, WhichTimePoint)	
{
	## Extract Important Quantities (Extinction Prob, and PDFs and Means 
	ParamVec	= PDF_PROG$ParamVec
	Mode		= PDF_PROG$Mode
	Tau 		= PDF_PROG$Tau
	ExtProb 	= PDF_PROG$PDFs	[1, WhichTimePoint]
	Multiplier	= (1 / (1 - ExtProb))
	
	## Remove Extinction Probability and Normalise PDFs 
	PDF	= PDF_PROG$PDFs	[, WhichTimePoint]
	PDF	= PDF[-1] * Multiplier	
	
	Mean 	= PDF_PROG$Means[WhichTimePoint] * Multiplier
	Var		= sum(PDF *	((1:Tau) - Mean)^2)
	
	PDFlist 							= list()
	PDFlist[["PDF"				]]		= PDF
	PDFlist[["ExtProb"			]] 		= ExtProb
	PDFlist[["Mean"				]]		= Mean
	PDFlist[["Var"				]]		= Var	
	PDFlist[["Mode"				]]		= Mode
	PDFlist[["TimeStep"			]]		= TimeStep
	PDFlist[["ParamVec"			]]		= ParamVec
	PDFlist[["MakeP_tstep_call"	]] 		= match.call()
	
	return(PDFlist)
}


CalculatePDFEvolutions = function()
{
	cat(paste0("Calculate transition matrix\n"))
	TransitionMatrix 	= MakeTransitionMatrix(Tau = Tau, Pi = Pi, delta = delta, H = H)
	cat(paste0("Calculate matrix exponential\n"))
	eAt					= CalculateMatrixExponential(TransitionMatrix = TransitionMatrix, TimeStep = TimeStep )
	cat(paste0("Prob Dist Evolutions\n"))
	for (Freq in 1:SwitchFreq)
	{
		cat(paste("F", Freq, Sys.time(), "\n"))
		assign(paste0("CloneEvo_StartOff_", Freq), 
				PDFEvolution(Tau = Tau, Pi = Pi, delta = delta, H = H, 
						LegallyDeadThreshold = 1, RubberNeckThreshold = 1/100	, 
						TimeStep = TimeStep, Duration = MaxStochasticDuration,
						eAt = eAt, TransitionMatrix = TransitionMatrix,
						StartOff = Freq)
		)
	}
	List = mget(paste("CloneEvo_StartOff_", 1:SwitchFreq, sep =""))
	List[["eAt"]] = eAt
	List[["TransitionMatrix"]] = TransitionMatrix
	
	return(List)		
}

















