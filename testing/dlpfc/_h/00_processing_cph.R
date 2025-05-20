## Load the package at the top of your script
library("sessioninfo")

## Reproducibility information
print('Reproducibility information:')
Sys.time()
proc.time()
options(width = 120)
session_info()