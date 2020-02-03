# KIPIK
Fractional inhibition data and inhibitor library size analysis for Watson et al. (2020)

To run the prediction analysis:

*Checkout this repository from gitHub
*Open R
*Set your working directory to the Watson_2019 directory
*Run the "Predictions.R" script

This script will use all available CPUs on your machine (less one).  Typically takes several hours to complete on a 24 core machine.  To run faster, but with fewer samples, edit the text file "Predictions.R" so that the number of repeats is a low value, 10 for example:

```
nreps = 10
```
