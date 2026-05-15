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

DOI: 
