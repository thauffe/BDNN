# Proboscidean BDNN analyses using PyRate

This directory contains all data to reproduce the proboscidean analyses. These include:
* `proboscideans.py`, ten replicates of dated fossil occurrences to account for their age uncertainty. Sourced from [Cantalapiedra et al. (2021)](https://www.nature.com/articles/s41559-021-01498-w).
* `Traits_*.txt`, ten trait files containing ecomorphological traits summarized on two NMDS axes (NMDS1 and 2), the phylogenetic reletedness summarized by two phylogenetic eigenvectors (PVR1 and 2), and the geographic distribution (America, Africa, Eurasia, and Island).
* `FilterTaxa.txt`, specifying which species of the fossil occurrences should be used. Ten species were excluded due to missing traits and singleton occurrences before 40 Ma.
* `sampling_epochs.txt`, allowing for episodic changing fossil smapling rates in the analysis.
* `FixedTimeShifts.txt`, setting time windows to run the BDNN not in bins with a duration of 1 million years but with geological time scale to reduce run time.
* `TimeVar.txt`, the time-variable diversification predictors of paleotemperature, presence of humans, and the emergence of open grasslands. All these are used only to initiate the BDNN analysis, while the python script in the [Make_species_specific_timevars](https://github.com/thauffe/BDNN/tree/main/Proboscidea/Make_species_specific_timevars) directory adds the species-specific trajectory of paleotemperature and the continent-specific appearence of grasslands and spatial-temporal overlap with humans.
* `Cont_stats_*.txt`, ten files containing the mean and standard deviation of continuous diversification predictors, which are used convert z-standardized predictors back to their original scale for plotting the predictor's effect.


Instructions on how to setup python are available at the [PyRate website](https://github.com/dsilvestro/PyRate/blob/master/tutorials/pyrate_tutorial_0.md)


