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

The */.../* should be replace with the path to this GitHub repository. This Bayesian analysis lasts ca. two weeks. It can be be made faster by specifiying less MCMC generation in the `-n` flag. On Linux systems, the flag `-thread 4 0` allows parallel calculation on, for instance, four CPUs.


## Post-processing

### Visualizing the effect of diversification predictors on rates

We can visualize the effect of a single predictor or the interaction between two predictors on speciation and extinction rates using partial dependence plots. These marginalize over the effects of the remaining predictors. The following command will produce a PDF file called *proboscideans_2_G_BDS_BDNN_16_8TVcb_PDP.pdf* in the pyrate_mcmc_logs sub-directory.
```
python /Path to PyRate/PyRate.py \
-plotBDNN_effects /.../BDNN/Proboscidea/PyRate_analyses/pyrate_mcmc_logs/proboscideans_2_G_BDS_BDNN_16_8TVcb_mcmc.log \
-plotBDNN_transf_features /.../BDNN/Proboscidea/PyRate_analyses/Cont_stats_2.txt \
-BDNN_groups "{\"Geography\": [\"Africa\", \"America\", \"Eurasia\", \"Island\"]}" \
-b 0.25 -resample 1000
```

This can be speed-up on Linux and Mac systems with the flag `-thread 4 0` or by using less MCMC samples, which is specified by the `-resample` argument.

### Ranking the predictors' importance

Finally, we rank the importance of the predictors on the inferred speciation and extinction rates.

```
python /Path to PyRate/PyRate.py \
-BDNN_pred_importance /.../BDNN/Proboscidea/PyRate_analyses/pyrate_mcmc_logs/proboscideans_2_G_BDS_BDNN_16_8TVcb_mcmc.log \
-BDNN_groups "{\"Geography\": [\"Africa\", \"America\", \"Eurasia\", \"Island\"]}" \
-BDNN_pred_importance_window_size 0.1 \
-b 0.25 -resample 1000
```

This step is slow but it can be speed-up on Linux and Mac systems with the flag `-thread 4 0`, by using less MCMC samples, specified by the `-resample` argument, and with `-BDNN_pred_importance_only_main`, which calculates only the importance of individual predictors but not their interactions.

This step will produce five text files and one PDF with plots in the pyrate_mcmc_logs sub-directory. The coefficient of rate variation is safed in the file *proboscideans_2_G_BDS_BDNN_16_8TVcb_coefficient_of_rate_variation.csv*. The text file *proboscideans_2_G_BDS_BDNN_16_8TVcb_ex_predictor_influence.csv* (and *sp* for speciation) summarizes the results of three artificial explainable intelligence metrics and show their consensus rank. The text file *proboscideans_2_G_BDS_BDNN_16_8TVcb_ex_shap_per_species.csv* (and *sp* for speciation) shows with Shapley additive value (SHAP) for each species how much their species-time-specific rates deviate from the average acroos all species due to the effect of each predictor. The latter is used to produce the *proboscideans_2_G_BDS_BDNN_16_8TVcb_contribution_per_species_rates.pdf* plot, which shows the rates and SHAP values per species.

![SHAP values extinction](https://github.com/thauffe/BDNN/blob/main/Proboscidea/Plots/Contribution_extinction.png)

