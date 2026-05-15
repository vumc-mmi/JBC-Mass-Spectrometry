# JBC-Mass-Spectrometry

Launch the App: Open the R script Vulcanoplot v16.R file in RStudio and click 'Run App', or run via the console:

R
shiny::runApp('path_to_your_script/R script Vulcanoplot v16.R')


Upload Data: Upload a CSV file containing your protein groups.

Requirement: Columns for LFQ intensities must contain the string "LFQ intensity".

Requirement (Optional): For tryptic filtering, columns named "Amino acid before" and "Last amino acid" are required.

Define Groups: Select the columns corresponding to your experimental and control groups.

Configure Parameters: Adjust the normalization, imputation settings, and significance thresholds.

Generate & Export: Click "Generate Plot". The results can be exported as a high-resolution PDF and the full statistical table as a CSV.

Groups present in the data set:
						
A	Wt	 Mycobacterium marinum	Surface fraction (Genapol)	FASP + OASIS HLB (10mg), ISD	repl1	(A) Wt_Repl1  
B	D PecA	 Mycobacterium marinum PecA Mutant	Surface fraction (Genapol)	FASP + OASIS HLB (10mg), ISD	repl1	(B) PecA_Repl1  
C	D PecA/B/C	 Mycobacterium marinum PecA/PecB/PecC Triple Mutant	Surface fraction (Genapol)	FASP + OASIS HLB (10mg), ISD	repl1	(C) PecABC_Repl1  
D	Wt	 Mycobacterium marinum	Surface fraction (Genapol)	FASP + OASIS HLB (10mg), ISD	repl2	(D) Wt_Repl2  
E	D PecA	 Mycobacterium marinum PecA Mutant	Surface fraction (Genapol)	FASP + OASIS HLB (10mg), ISD	repl2	(E) PecA_Repl2  
F	D PecA/B/C	 Mycobacterium marinum PecA/PecB/PecC Triple Mutant	Surface fraction (Genapol)	FASP + OASIS HLB (10mg), ISD	repl2	(F) PecABC_Repl2  
G	Wt	 Mycobacterium marinum	Surface fraction (Genapol)	FASP + OASIS HLB (10mg), ISD	repl3	(G) Wt_Repl3  
H	D PecA	 Mycobacterium marinum PecA Mutant	Surface fraction (Genapol)	FASP + OASIS HLB (10mg), ISD	repl3	(H) PecA_Repl3  
I	D PecA/B/C	 Mycobacterium marinum PecA/PecB/PecC Triple Mutant	Surface fraction (Genapol)	FASP + OASIS HLB (10mg), ISD	repl3	(I) PecABC_Repl3  


[![DOI](https://zenodo.org/badge/1239865762.svg)](https://doi.org/10.5281/zenodo.20209201)
