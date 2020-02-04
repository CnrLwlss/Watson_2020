# KIPIK
Fractional inhibition data and inhibitor library size analysis for Watson et al. (2020)

To run the prediction analysis, please install the latest version of R (3.6.0 or greater):

* Checkout this repository from gitHub
* Open R
* Set your working directory to the Watson_2020 directory
* Run the "Predictions.R" script

This script will use all available CPUs on your machine (less one).  Typically takes several hours to complete on a 24 core machine.  To run a small, demo version of analysis that is faster, but with fewer samples, before carrying out the final step, edit the text file "Predictions.R" so that the number of repeats is a low value, 100 for example:

```
nreps = 10
```

Once either the demo or the full analysis have completed, you should find that the R code has run without errors and two output files ("Sampling_fixed.pdf" & "Sampling_varying.pdf") generated in the root directory.
