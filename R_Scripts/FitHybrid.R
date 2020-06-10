## Author: dlaydon 

## script combines: i) CloneTrajectories.R; ii) QuickGridSearch.R; iii) OptimmizeFit.R. 
## Script returns an estimate of the rate of infectious spread (de novo infection) for a particular HTLV-1 infected person data set, given the rate of mitotic spread (infected cell proliferation).
# i) 	CloneTrajectories.R calculates look-up tables for stochastically modelled clones
# ii) 	QuickGridSearch.R performs crude parameter grid search, running model using look-up tables calculated in i)
# iii) 	OptimmizeFit.R performs 1-dimentional optimization using parameter range and initial values calculated in ii), using look up tables from i). 

RootDir 		= "." ### Change as appropriate
R_Script_Dir 	= file.path(RootDir, "R_Scripts")
OutputDir 		= file.path(RootDir, "Output")

## script combines: i) CloneTrajectories; ii) QuickGridSearch; iii) OptimmizeFit.R

dta 				= 8 		# patient data set number to fit.
OutputStringPrefix 	= "Test"

Pi_dash1 	= 0.0316	# rate of mitotic spread
delta		= Pi_dash1	# rate of cell death.
SwitchFreq	= 460		# Threshold frequency F in manuscript: Frequency above and below which clones are respectively modelled deterministically and stochastically.
Tau 		= 1500		# maximum frequency in individual clone state space for clone birth-death process. 
TimeStep 	= 1			# units in days
Duration	= 3133		# units in days. For how long will hybrid run?
FreqCutOff 	= Inf		# Clones of freqeuncy greater than this will be excluded (set to Inf if modelling all clones).

# For quick grid search = range of values to be considered
NumberOfGuesses = 25
BestGuess 		= 1e-9
LowerBound 		= 1e-11
UpperBound 		= 3e-9	

source(file.path(R_Script_Dir, "CloneTrajectories.R"))
source(file.path(R_Script_Dir, "QuickGridSearch.R"))
source(file.path(R_Script_Dir, "OptimizeFit.R"))
























