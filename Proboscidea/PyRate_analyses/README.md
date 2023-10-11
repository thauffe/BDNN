# Proboscidean BDNN analyses using PyRate

## The data

This directory contains all data to reproduce the proboscidean analyses with the steps shown below. These include:
* `proboscideans.py`, ten replicates of dated fossil occurrences to account for their age uncertainty. Sourced from [Cantalapiedra et al. (2021)](https://www.nature.com/articles/s41559-021-01498-w).
* `Traits_*.txt`, ten trait files containing ecomorphological traits summarized on two NMDS axes (NMDS1 and 2), the phylogenetic reletedness summarized by two phylogenetic eigenvectors (PVR1 and 2), and the geographic distribution (America, Africa, Eurasia, and Island).
* `FilterTaxa.txt`, specifying which species of the fossil occurrences should be used. Ten species were excluded due to missing traits and singleton occurrences before 40 Ma.
* `sampling_epochs.txt`, allowing for episodic changing fossil smapling rates in the analysis.
* `FixedTimeShifts.txt`, setting time windows to run the BDNN not in bins with a duration of 1 million years but with geological time scale to reduce run time.
* `TimeVar.txt`, the time-variable diversification predictors of paleotemperature, presence of humans, and the emergence of open grasslands. All these are used only to initiate the BDNN analysis, while the python script in the [Make_species_specific_timevars](https://github.com/thauffe/BDNN/tree/main/Proboscidea/Make_species_specific_timevars) directory adds the species-specific trajectory of paleotemperature and the continent-specific appearence of grasslands and spatial-temporal overlap with humans.
* `pyrate_mcmc_logs/proboscideans_*_G_BDS_BDNN_16_8TVcb.pkl`, python objects where species-specific time series of paleotemperature and the spatial-temporal overlap with grassland and humans is already included. These files can only be viewed in python because they serialized.
* `Cont_stats_*.txt`, ten files containing the mean and standard deviation of continuous diversification predictors, which are used convert z-standardized predictors back to their original scale for plotting the predictor's effect.


For the first replicate of the proboscidean dataset, the sub-directory [pyrate_mcmc_logs](https://github.com/thauffe/BDNN/tree/main/Proboscidea/PyRate_analyses/pyrate_mcmc_logs) contains the MCMC log-files of the BDNN inference and the output of the post-processing steps. These use approaches from explainable artificial intelligence to obtain the importance of diversification predictors.


## BDNN inference

Instructions on how to setup python for PyRate are available at the [PyRate website](https://github.com/dsilvestro/PyRate/blob/master/tutorials/pyrate_tutorial_0.md)


For replicate 2, the BDNN model can be inferred with the following command:
```
python /Path to PyRate/PyRate.py \
/.../BDNN/Proboscidea/PyRate_analyses/proboscideans.py \
-j 2 \
-filter_taxa /.../BDNN/Proboscidea/PyRate_analyses/FilterTaxa.txt \
-trait_file /.../BDNN/Proboscidea/PyRate_analyses/Traits_2.txt \
-A 0 \
-mG \
-qShift /.../BDNN/Proboscidea/PyRate_analyses/sampling_epochs.txt \
-BDNNmodel 1 \
-BDNNpklfile /.../BDNN/Proboscidea/PyRate_analyses/pyrate_mcmc_logs/proboscideans_2_G_BDS_BDNN_16_8TVcb.pkl \
-BDNNtimevar /.../BDNN/Proboscidea/PyRate_analyses/TimeVar.txt \
-fixShift /.../BDNN/Proboscidea/PyRate_analyses/FixedTimeShifts.txt \
-n 50000000 -s 10000 -p 1000000
```

The */.../* should be replace with the path to this GitHub repository. This Bayesian analysis lasts ca. two weeks. It can be be made faster by specifiying less MCMC generation in the `-n` flag. On Linux systems, the flag `-thread 4 0` allows parallel calculation on four CPUs.


## Post-processing

### Visualizing the effect of diversification predictors on rates

```
python3 /home/thauffe/BDNN/PyRate/PyRate.py \
-plotBDNN_effects /data/users/thauffe/BDNN/Proboscideans/PyRateAnalyses40Ma/Humans_Island_SpTemp_Grass_Juan/pyrate_mcmc_logs/proboscideans_${SLURM_ARRAY_TASK_ID}_G_BDS_BDNN_16_8TVcb_mcmc.log \
-plotBDNN_transf_features /data/users/thauffe/BDNN/Proboscideans/PyRateAnalyses40Ma/Humans_Island_SpTemp_Grass_Juan/Cont_stats_${SLURM_ARRAY_TASK_ID}.txt \
-BDNN_groups '{"Geography": ["Africa", "America", "Eurasia", "Island"]}' \
-b 0.25 -resample 100 -thread 10 0
```
